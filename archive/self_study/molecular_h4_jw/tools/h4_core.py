"""Standalone H4 teaching runtime; no file access, optimization, or accounts.

Blocked spin order alpha0..3,beta0..3; little endian; BPDE ancilla q8.
All circuit builders return unmeasured unitary circuits for exact verification.
"""
import numpy as np
from scipy.linalg import eigh
from scipy.sparse import csr_matrix
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import UCRYGate, UnitaryGate, PauliEvolutionGate
from qiskit.quantum_info import Operator, Statevector, SparsePauliOp

DT = 0.1
TAU = 2*np.pi


def fock_permutation():
    """Interleaved -> blocked: signed fermionic permutation, not bit swaps."""
    mapping = [0,4,1,5,2,6,3,7]
    indices, signs = [], []
    for n in range(256):
        occupied = [mapping[p] for p in range(8) if n >> p & 1]
        indices.append(sum(1 << p for p in occupied))
        inversions = sum(a > b for j,a in enumerate(occupied) for b in occupied[j+1:])
        signs.append((-1)**inversions)
    return np.asarray(indices), np.asarray(signs)


def to_blocked(matrix):
    indices, signs = fock_permutation()
    out = np.empty_like(matrix)
    out[np.ix_(indices,indices)] = signs[:,None]*matrix*signs[None,:]
    return out


def fermion_operators():
    operators = []
    for p in range(8):
        columns = [n for n in range(256) if n >> p & 1]
        rows = [n ^ (1 << p) for n in columns]
        values = [(-1)**((n & ((1 << p)-1)).bit_count()) for n in columns]
        operators.append(csr_matrix((values,(rows,columns)),shape=(256,256)))
    return operators


def hamiltonian_from_integrals(h1, eri):
    """Chemist ERI (pq|rs); independent explicit fermionic construction."""
    a = fermion_operators()
    h = csr_matrix((256,256),dtype=float)
    for p in range(4):
        for q in range(4):
            for sigma in range(2):
                h += h1[p,q]*(a[p+4*sigma].T@a[q+4*sigma])
            for r in range(4):
                for s in range(4):
                    if abs(eri[p,q,r,s]) < 1e-15:
                        continue
                    for sigma in range(2):
                        for tau in range(2):
                            h += .5*eri[p,q,r,s]*(a[p+4*sigma].T@a[r+4*tau].T@a[s+4*tau]@a[q+4*sigma])
    return h.toarray()


def spin_operators():
    a = fermion_operators()
    na = sum(x.T@x for x in a[:4]).toarray()
    nb = sum(x.T@x for x in a[4:]).toarray()
    sz = (na-nb)/2
    plus = sum(a[p].T@a[p+4] for p in range(4)).toarray()
    s2 = sz@sz + (plus@plus.T+plus.T@plus)/2
    return na, nb, s2


def solve_states(h):
    """Resolve spin first, then energy. T denotes lowest S=1, Ms=0 state."""
    sector = np.array([n for n in range(256)
                       if (n & 15).bit_count()==2 and (n >> 4).bit_count()==2])
    na,nb,s2 = spin_operators()
    spins, vectors = eigh(s2[np.ix_(sector,sector)])
    states, energies = {}, {}
    for prefix, spin, count in [('S',0,2),('T',2,1)]:
        basis = vectors[:,abs(spins-spin)<1e-8]
        e,v = eigh(basis.T@h[np.ix_(sector,sector)]@basis)
        for j in range(count):
            name = prefix+str(j) if prefix=='S' else 'T'
            psi = np.zeros(256)
            psi[sector] = basis@v[:,j]
            psi *= np.sign(psi[np.argmax(abs(psi))])
            states[name], energies[name] = psi,float(e[j])
    hf = np.zeros(256); hf[51] = 1  # q0,q1,q4,q5 occupied
    states['HF'] = hf
    return states, energies, (na,nb,s2)


def rotation_angles(vector, target):
    """Binary-tree real-amplitude preparation; signs retained at leaf level."""
    groups = np.asarray(vector,float).reshape(-1,2,2**target)
    if target == 0:
        left,right = groups[:,0,0], groups[:,1,0]
    else:
        left,right = np.linalg.norm(groups[:,0],axis=1),np.linalg.norm(groups[:,1],axis=1)
    return (2*np.arctan2(right,left)).tolist()


def prepare_real(vector, name='prep'):
    vector = np.asarray(vector)
    if np.max(abs(vector.imag)) > 1e-12 or not np.isclose(np.linalg.norm(vector),1):
        raise ValueError('Expected normalized real amplitudes')
    n = int(np.log2(len(vector)))
    qc = QuantumCircuit(n,name=name)
    for target in reversed(range(n)):
        angles = rotation_angles(vector.real,target)
        if target==n-1:
            qc.ry(angles[0],target)
        else:
            qc.append(UCRYGate(angles),[target]+list(range(target+1,n)))
    return qc


def prepare_pair(ground, excited):
    """A = diag(U0,Uj) (H_a tensor I); coherent branch phase is fixed."""
    qc = QuantumCircuit(9,name='A_0j')
    qc.h(8)
    for target in reversed(range(8)):
        angles = rotation_angles(ground,target)+rotation_angles(excited,target)
        qc.append(UCRYGate(angles),[target]+list(range(target+1,8))+[8])
    return qc


class Propagator:
    """Spectral acceleration evaluates the repeated circuit, not its gate cost."""
    def __init__(self, matrix):
        self.matrix = np.asarray(matrix,complex)
        if unitary_defect(self.matrix)>1e-8:
            raise ValueError('Nonunitary reference MPO cannot be a circuit propagator')
        # All one-step operators here are near I: cos(phase)>0. On this arc
        # sin(phase) is one-to-one. Hermitian diagonalization keeps a unitary
        # eigenbasis even for the large, exactly degenerate spin multiplets.
        cosine=(self.matrix+self.matrix.conj().T)/2
        if np.min(np.linalg.eigvalsh(cosine))<=0:
            raise ValueError('One-step phase must lie strictly within (-pi/2,pi/2)')
        sine=(self.matrix-self.matrix.conj().T)/(2j)
        _,self.vectors=eigh(sine)
        diagonal=self.vectors.conj().T@self.matrix@self.vectors
        self.phases=np.angle(np.diag(diagonal))
        if relative_error(self.power(1),self.matrix)>1e-8:
            raise ValueError('Spectral reconstruction failed')

    def power(self, steps):
        return (self.vectors*np.exp(1j*self.phases*int(steps)))@self.vectors.conj().T

    def evolve(self, vector, steps):
        return self.vectors@(np.exp(1j*self.phases*int(steps))*(self.vectors.conj().T@vector))

    def survival(self, vector, steps):
        weights = abs(self.vectors.conj().T@vector)**2
        return np.sum(weights*np.exp(1j*self.phases*int(steps)))


def exact_unitary(h,dt=DT):
    e,v = eigh(h)
    return (v*np.exp(-1j*e*dt))@v.T.conj()


def relative_error(a,b):
    return float(np.linalg.norm(a-b)/np.linalg.norm(b))


def unitary_defect(u):
    return float(np.linalg.norm(u.conj().T@u-np.eye(len(u)))/np.sqrt(len(u)))


def trotter_step(labels,coefficients,slices=1,dt=DT):
    """One full dt, symmetric Pauli splitting, including the identity term."""
    size=256; indices=np.arange(size); matrix=np.eye(size,dtype=complex)
    terms=list(zip(labels,coefficients))
    for label,coefficient in terms+terms[::-1]:
        mask=0; phase=np.ones(size,complex)
        for q,axis in enumerate(reversed(str(label))):
            bit=(indices>>q)&1
            if axis in 'XY': mask |= 1 << q
            if axis=='Y': phase *= 1j*(1-2*bit)
            elif axis=='Z': phase *= 1-2*bit
        applied=np.empty_like(matrix); applied[indices^mask]=phase[:,None]*matrix
        angle=float(np.real(coefficient))*dt/(2*slices)
        matrix=np.cos(angle)*matrix-1j*np.sin(angle)*applied
    return np.linalg.matrix_power(matrix,slices)


def trotter_microcircuit(labels,coefficients,slices=1,dt=DT):
    qc=QuantumCircuit(8,name=f'T2micro_s{slices}')
    terms=list(zip(labels,coefficients))
    for label,coefficient in terms+terms[::-1]:
        op=SparsePauliOp(str(label),coeffs=[complex(coefficient)])
        qc.append(PauliEvolutionGate(op,time=dt/(2*slices)),range(8))
    return qc


def author_circuit(data):
    """Unmodified exported 4x4 gates, in original interleaved orbital order."""
    qc=QuantumCircuit(8,name='author_d8_interleaved')
    for layer in range(1,9):
        for unit,start in enumerate(range((layer-1)%2,7,2),1):
            gate=data[f'author_gate_{layer}_{unit}']
            np.testing.assert_allclose(gate.conj().T@gate,np.eye(4),atol=1e-10)
            qc.unitary(gate,[start,start+1],label=f'G{layer}.{unit}')
    return qc


def orbital_bridge():
    """Signed adjacent fermionic SWAPs: interleaved -> blocked."""
    desired=[0,2,4,6,1,3,5,7]; current=list(range(8))
    qc=QuantumCircuit(8,name='P_fock')
    for target,mode in enumerate(desired):
        position=current.index(mode)
        while position>target:
            qc.swap(position-1,position);qc.cz(position-1,position)
            current[position-1],current[position]=current[position],current[position-1]
            position-=1
    return qc


def bpde_circuit(prep, matrix, beta=0):
    qc=prep.copy(); qc.unitary(matrix,range(8),label='U(t)')
    qc.p(beta,8); qc.compose(prep.inverse(),inplace=True)
    return qc


def bpde_probabilities(prop, ground, excited, steps, beta=0, branch_matrices=None):
    v0=prop.evolve(ground,steps); v1=prop.evolve(excited,steps)
    a0=np.vdot(ground,v0); a1=np.vdot(excited,v1)
    pc=abs(a0+np.exp(1j*beta)*a1)**2/4
    if branch_matrices is None:
        return float(np.clip(pc,0,1))
    u0,uj=branch_matrices
    b0=u0.conj().T@v0; b1=uj.conj().T@v1
    pb=(1+np.real(np.exp(1j*beta)*np.vdot(b0,b1)))/2
    return float(np.clip(pb,0,1)),float(np.clip(pc,0,1))


def gap_probability(gap,steps,beta):
    return (1+np.cos(beta-np.asarray(gap)*DT*steps))/2


def schedule_comparison():
    stages=[2**r for r in range(15)]
    return {name:[(L,b) for L in lengths for b in [0,np.pi/2]]
            for name,lengths in [('short',[1]*15),('staged',stages),('long',[16384]*15)]}


def bayes_grid(grid, settings, truth, shots=256, seed=1, expected=False, energy=False):
    rng=np.random.default_rng(seed)
    # energy and gap have identical physical minus sign; different grid units only.
    model=np.array([np.clip(gap_probability(grid,L,b),1e-12,1-1e-12) for L,b in settings])
    probabilities=np.clip([truth(L,b) for L,b in settings],0,1)
    counts=shots*probabilities if expected else rng.binomial(shots,probabilities)
    increments=counts[:,None]*np.log(model)+(shots-counts[:,None])*np.log1p(-model)
    logs=np.cumsum(increments,axis=0); logs-=logs.max(axis=1,keepdims=True)
    history=np.exp(logs); history/=history.sum(axis=1,keepdims=True)
    records=np.array([(L,b,n,shots,p) for (L,b),n,p in zip(settings,counts,probabilities)])
    return history, records


def posterior_summary(grid,p):
    lo,hi=np.interp([.025,.975],np.cumsum(p),grid)
    return {'map':float(grid[np.argmax(p)]),'lo':float(lo),'hi':float(hi),'width':float(hi-lo)}


def repeated_rmse(grid,settings,truth,target,trials=64,shots=256,seed=1):
    model=np.array([np.clip(gap_probability(grid,L,b),1e-12,1-1e-12) for L,b in settings])
    probabilities=np.clip([truth(L,b) for L,b in settings],0,1)
    counts=np.random.default_rng(seed).binomial(shots,probabilities,size=(trials,len(settings)))
    logs=counts@np.log(model)+(shots-counts)@np.log1p(-model)
    estimates=grid[np.argmax(logs,axis=1)]
    return float(np.sqrt(np.mean((estimates-target)**2))),estimates


def inverse_qft(m):
    qc=QuantumCircuit(m,name='QFTdg')
    for j in range(m//2): qc.swap(j,m-1-j)
    for j in range(m):
        for k in range(j): qc.cp(-np.pi/2**(j-k),j,k)
        qc.h(j)
    return qc


def controlled_matrix(matrix):
    controlled=np.eye(512,dtype=complex)
    controlled[1::2,1::2]=matrix
    return controlled


def qpe_circuit(prep,prop,m=4):
    qc=QuantumCircuit(m+8)
    qc.compose(prep,range(m,m+8),inplace=True)
    for r in range(m):
        qc.h(r)
        qc.append(UnitaryGate(controlled_matrix(prop.power(2**r)),label=f'cU^{2**r}'),[r]+list(range(m,m+8)))
    qc.compose(inverse_qft(m),range(m),inplace=True)
    return qc


def qpe_formula(psi,h,m=4):
    e,v=eigh(h); weights=abs(v.conj().T@psi)**2
    phi=(-e*DT/TAU)%1; M=2**m
    delta=phi[:,None]-np.arange(M)[None,:]/M
    # Dirichlet kernel; sinc handles exact grid phases by continuity.
    amplitudes=np.sinc(M*delta)/np.sinc(delta)
    return weights@abs(amplitudes)**2


def iqpe_distribution(psi,prop,m=4,reprepare=False):
    branches=[(0,1.,psi)]
    for ell,r in enumerate(reversed(range(m))):
        matrix=prop.power(2**r); new=[]
        for y,weight,state in branches:
            incoming=psi if reprepare else state
            beta=-TAU*2**r*y/2**m
            for bit in [0,1]:
                out=(incoming+(-1)**bit*np.exp(1j*beta)*matrix@incoming)/2
                p=float(np.vdot(out,out).real)
                if p>1e-15: new.append((y+bit*2**ell,weight*p,out/np.sqrt(p)))
        branches=new
    result=np.zeros(2**m)
    for y,p,_ in branches: result[y]+=p
    return result


def bpe_circuit(prep,matrix,beta=0):
    qc=QuantumCircuit(9);qc.compose(prep,range(8),inplace=True);qc.h(8)
    qc.append(UnitaryGate(controlled_matrix(matrix),label='cU(t)'),[8]+list(range(8)))
    qc.p(beta,8);qc.h(8)
    return qc


def compile_counts(qc):
    compiled=transpile(qc,basis_gates=['ry','rz','h','x','cx'],optimization_level=1,seed_transpiler=7)
    return compiled,{'cx':int(compiled.count_ops().get('cx',0)),
                     'depth':int(compiled.depth()),'qubits':compiled.num_qubits}


def regenerate_chemistry():
    """Optional independent PySCF regeneration; molecular orbital phases may differ."""
    from pyscf import gto,scf,ao2mo,fci
    mol=gto.M(atom=[('H',(0,y,0)) for y in [3,1,-1,-3]],basis='sto-3g',unit='Angstrom',spin=0,verbose=0)
    mf=scf.RHF(mol).run(conv_tol=1e-12,max_cycle=150)
    if not mf.converged: raise RuntimeError('RHF did not converge')
    c=mf.mo_coeff; h1=c.T@mf.get_hcore()@c
    eri=ao2mo.restore(1,ao2mo.kernel(mol,c),4)
    h=hamiltonian_from_integrals(h1,eri)
    states,energies,_=solve_states(h)
    singlet_solver=fci.direct_spin0.FCI(mol)
    singlet_solver.conv_tol=1e-12
    singlets,_=singlet_solver.kernel(h1,eri,4,(2,2),nroots=2)
    triplet,_=fci.direct_spin1.kernel(h1,eri,4,(3,1),tol=1e-12)
    np.testing.assert_allclose([energies['S0'],energies['S1'],energies['T']],
                               [*singlets,triplet],atol=1e-9,rtol=0)
    return energies,float(mol.energy_nuc())
