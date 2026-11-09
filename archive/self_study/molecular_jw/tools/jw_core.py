"""Small dense verification helpers; optimization itself is the upstream Julia code.
Qiskit little endian, system q0..q3; BPDE ancilla q4.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import UnitaryGate
from qiskit.quantum_info import Operator, Statevector, SparsePauliOp

TAU = 2*np.pi

def state_prep(name, theta):
    qc = QuantumCircuit(4, name=f'prep_{name}')
    if name in ('T+', 'T-'):
        for q in ([0,1] if name == 'T+' else [2,3]):
            qc.ry(np.pi, q)
        return qc
    base, pivot, angle = {
        'S0': (5,1,theta), 'S2': (5,1,theta+np.pi),
        'S1': (6,0,np.pi/2), 'T': (6,0,-np.pi/2),
        'HF': (5,1,0.),
    }[name]
    for q in range(4):
        if base >> q & 1:
            qc.ry(np.pi, q)
    qc.ry(angle, pivot)
    for q in range(4):
        if q != pivot:
            qc.cx(pivot, q)
    return qc

def joint_prep(target, theta):
    """Exactly (|0> S0 + |1> target)/sqrt(2), including relative phase."""
    c = int(target in ('S1','T'))
    angle = {'S2':theta+np.pi, 'S1':np.pi/2, 'T':-np.pi/2}[target]
    qc = QuantumCircuit(5, name=f'joint_S0_{target}')
    # The minus compensates the open-shell encoder's branch phase.
    qc.ry(-np.pi/2 if c else np.pi/2, 4)
    if c:
        qc.cx(4,1)
    # Exact controlled angle selection with two Ry and two CX.
    qc.ry((theta+angle)/2,0); qc.cx(4,0)
    qc.ry((theta-angle)/2,0); qc.cx(4,0)
    # Common reversible encoder: b=q0, c=q1 -> configurations 5/10 or 6/9.
    qc.cx(0,1); qc.cx(0,2); qc.cx(0,3); qc.cx(1,0); qc.cx(2,0)
    qc.ry(np.pi,0); qc.ry(np.pi,2)
    return qc

def spin_operators():
    a=[]
    for q in range(4):
        op=np.zeros((16,16))
        for n in range(16):
            if n>>q&1:
                op[n^(1<<q),n]=(-1)**((n&((1<<q)-1)).bit_count())
        a.append(op)
    number=sum(x.T@x for x in a)
    sz=(a[0].T@a[0]+a[1].T@a[1]-a[2].T@a[2]-a[3].T@a[3])/2
    plus=a[0].T@a[2]+a[1].T@a[3]
    s2=sz@sz+(plus@plus.T+plus.T@plus)/2
    return number, sz, s2

def brickwall(arrays, depth):
    qc=QuantumCircuit(4, name=f'author_d{depth}')
    for layer in range(1,depth+1):
        starts=range(0,3,2) if layer%2 else range(1,3,2)
        for unit,start in enumerate(starts,1):
            gate=arrays[f'h2_d{depth}_MPO_depth{layer}_{unit}']
            np.testing.assert_allclose(gate.conj().T@gate,np.eye(4),atol=1e-10)
            qc.unitary(gate,[start,start+1],label=f'G{layer}.{unit}')
    return qc

def contract_author_mpo(arrays):
    """Actual exported cores; physical axes right,left; site 1 is q0."""
    cores=[arrays[f'reference_core_{j}'] for j in range(1,5)]
    tensor=cores[0]
    for core in cores[1:]:
        tensor=np.tensordot(tensor,core,axes=(-1,0))
    # Julia reshape is column-major: q0 is fastest within each physical group.
    return tensor.transpose(0,2,4,6,1,3,5,7).reshape(16,16,order='F')

def dense_trotter(h,dt,slices=100,order=2):
    p=SparsePauliOp.from_operator(Operator(h),atol=1e-12)
    terms=[(float(c.real), SparsePauliOp(label).to_matrix()) for label,c in zip(p.paulis.to_labels(),p.coeffs)]
    step=np.eye(16,dtype=complex)
    sequence=terms+terms[::-1] if order==2 else terms
    factor=2 if order==2 else 1
    for c,pauli in sequence:
        angle=c*dt/slices/factor
        step=(np.cos(angle)*np.eye(16)-1j*np.sin(angle)*pauli)@step
    return np.linalg.matrix_power(step,slices)

def rel_error(a,b):
    return float(np.linalg.norm(a-b)/np.linalg.norm(b))

def unitary_defect(u):
    return float(np.linalg.norm(u.conj().T@u-np.eye(len(u)))/np.sqrt(len(u)))

def controlled_matrix(u):
    # The first qarg is the low bit: control=0 -> I, control=1 -> U.
    cu=np.eye(2*len(u),dtype=complex)
    cu[1::2,1::2]=u
    return cu

def inverse_qft(m):
    qc=QuantumCircuit(m,name='QFT_dg')
    for j in range(m//2): qc.swap(j,m-1-j)
    for j in range(m):
        for k in range(j): qc.cp(-np.pi/2**(j-k),j,k)
        qc.h(j)
    return qc

def qpe_circuit(prep,u,m=3):
    qc=QuantumCircuit(m+4)
    qc.compose(prep,range(m,m+4),inplace=True)
    for r in range(m):
        qc.h(r)
        qc.append(UnitaryGate(controlled_matrix(np.linalg.matrix_power(u,2**r)),label=f'c-U^{2**r}'),[r]+list(range(m,m+4)))
    qc.compose(inverse_qft(m),range(m),inplace=True)
    return qc

def qpe_distribution(prep,u,m=3):
    return Statevector.from_instruction(qpe_circuit(prep,u,m)).probabilities(range(m))

def qpe_formula(phi,m):
    M=2**m
    return abs(np.mean(np.exp(1j*TAU*np.outer(np.arange(M),phi-np.arange(M)/M)),axis=0))**2

def iqpe_tree(psi,u,m=3,reprepare=False):
    branches=[(0,1.,np.asarray(psi,complex))]
    for ell,r in enumerate(reversed(range(m))):
        up=np.linalg.matrix_power(u,2**r);new=[]
        for y,weight,state in branches:
            incoming=psi if reprepare else state
            beta=-TAU*2**r*y/2**m
            for bit in (0,1):
                out=(incoming+(-1)**bit*np.exp(1j*beta)*up@incoming)/2
                prob=float(np.vdot(out,out).real)
                if prob>1e-14: new.append((y+bit*2**ell,weight*prob,out/np.sqrt(prob)))
        branches=new
    p=np.zeros(2**m)
    for y,w,_ in branches: p[y]+=w
    return p

def iqpe_dynamic(prep,u,m=3):
    from qiskit import QuantumRegister,ClassicalRegister
    a=QuantumRegister(1,'a');s=QuantumRegister(4,'s');c=ClassicalRegister(m,'bits')
    qc=QuantumCircuit(a,s,c);qc.compose(prep,s,inplace=True)
    for ell,r in enumerate(reversed(range(m))):
        if ell: qc.reset(a[0])
        qc.h(a[0]);qc.append(UnitaryGate(controlled_matrix(np.linalg.matrix_power(u,2**r)),label=f'c-U^{2**r}'),[a[0]]+list(s))
        for old in range(ell):
            with qc.if_test((c[old],1)): qc.p(-TAU*2**r*2**old/2**m,a[0])
        qc.h(a[0]);qc.measure(a[0],c[ell])
    return qc

def bpe_circuit(prep,u,L=1,beta=0):
    qc=QuantumCircuit(5);qc.compose(prep,range(4),inplace=True);qc.h(4)
    qc.append(UnitaryGate(controlled_matrix(np.linalg.matrix_power(u,L)),label=f'c-U^{L}'),[4,0,1,2,3])
    qc.p(beta,4);qc.h(4)
    return qc

def bpe_probability(psi,u,L,beta):
    amp=np.vdot(psi,np.linalg.matrix_power(u,L)@psi)
    return float((1+np.real(np.exp(1j*beta)*amp))/2)

def bpde_circuit(target,theta,u,L=1,beta=0):
    prep=joint_prep(target,theta)
    qc=prep.copy()
    qc.unitary(np.linalg.matrix_power(u,L),range(4),label=f'U^{L}')
    qc.p(beta,4);qc.compose(prep.inverse(),inplace=True)
    return qc

def bpde_probability(target,theta,u,L,beta,all_zero=True):
    vector=Statevector.from_instruction(bpde_circuit(target,theta,u,L,beta)).data
    return float(abs(vector[0])**2 if all_zero else sum(abs(vector[:16])**2))

def single_phase(phi,L,beta):
    return .5*(1+np.cos(TAU*L*np.asarray(phi)+beta))

def bayes(grid,settings,probability,shots=256,seed=20261108,expected=False):
    rng=np.random.default_rng(seed);logp=np.zeros_like(grid);history=[];records=[]
    for L,beta in settings:
        truth=float(np.clip(probability(L,beta),0,1))
        n0=shots*truth if expected else int(rng.binomial(shots,truth))
        p=np.clip(single_phase(grid,L,beta),1e-12,1-1e-12)
        logp+=n0*np.log(p)+(shots-n0)*np.log1p(-p);logp-=max(logp)
        post=np.exp(logp);post/=post.sum();history.append(post)
        records.append((L,beta,n0,shots,truth))
    return np.asarray(history),records

def summary(grid,p):
    lo,hi=np.interp([.025,.975],np.cumsum(p),grid)
    return float(grid[np.argmax(p)]),float(lo),float(hi)

def time_comparison(probability,target,label,seed=42):
    grid=(np.arange(32768)+.5)/32768
    plans={'L=1':[1]*6,'Staged':[1,2,4,8,16,32],'L=32':[32]*6}
    fig,axes=plt.subplots(2,3,figsize=(12,6))
    results={}
    for k,(name,lengths) in enumerate(plans.items()):
        settings=[(L,b) for L in lengths for b in [np.pi/5,np.pi/5-np.pi/2]]
        hist,records=bayes(grid,settings,probability,seed=seed+k)
        est,lo,hi=summary(grid,hist[-1]);widths=[summary(grid,p)[2]-summary(grid,p)[1] for p in hist[1::2]]
        for row in (0,1):
            axes[row,k].plot(grid,hist[-1]*len(grid));axes[row,k].axvline(target,c='C3',ls='--')
            axes[row,k].set(xlabel='Phase',ylabel='Posterior density',title=name)
        axes[1,k].set_xlim(target-.025,target+.025)
        print(name,'shots=',256*len(settings),'sum(R L)=',256*sum(lengths)*2,'MAP=',est,'95% width=',hi-lo)
        results[name]=(settings,hist,widths)
    fig.suptitle(label);plt.tight_layout();plt.show()
    fig,ax=plt.subplots(figsize=(7,3))
    for name,(_,_,widths) in results.items(): ax.semilogy(np.arange(1,7)*512,widths,'o-',label=name)
    ax.set(xlabel='Cumulative shots',ylabel='Global 95% interval width',title=label);ax.legend();plt.show()
    # L=32 aliases are exactly unresolved for the uniform prior.
    np.testing.assert_allclose(results['L=32'][1][-1].reshape(32,-1).sum(axis=1),1/32,atol=1e-9)
    for k,name in enumerate(['L=1','Staged']):
        settings=results[name][0];rng=np.random.default_rng(seed+100+k)
        probs=[probability(L,b) for L,b in settings]
        counts=rng.binomial(256,np.clip(probs,0,1),size=(64,len(settings)))
        model=np.array([np.clip(single_phase(grid,L,b),1e-12,1-1e-12) for L,b in settings])
        logp=counts@np.log(model)+(256-counts)@np.log1p(-model);logp-=logp.max(axis=1,keepdims=True)
        post=np.exp(logp);post/=post.sum(axis=1,keepdims=True)
        estimates=grid[np.argmax(post,axis=1)]
        err=(estimates-target+.5)%1-.5
        print(name,'64-trial circular RMSE=',np.sqrt(np.mean(err**2)))
    return results

def exact_two_qubit_synthesis(matrix):
    """Avoid default KAK specialization tolerances; preserve the input global phase.
    Try 0/1/2/3 CX decompositions of the unspecialized Weyl target, accept only
    a matrix-verified representation. This changes compilation, not optimization.
    """
    from qiskit.circuit.library import CXGate
    from qiskit.synthesis import TwoQubitBasisDecomposer, TwoQubitWeylDecomposition, OneQubitEulerDecomposer
    target=TwoQubitWeylDecomposition(matrix,fidelity=None)
    decomposer=TwoQubitBasisDecomposer(CXGate(),basis_fidelity=1.,euler_basis='ZSX',pulse_optimize=False)
    euler=OneQubitEulerDecomposer('ZSX')
    methods=[decomposer.decomp0,decomposer.decomp1,decomposer.decomp2_supercontrolled,decomposer.decomp3_supercontrolled]
    for count,method in enumerate(methods):
        locals_=method(target._inner_decomposition)
        circuit=QuantumCircuit(2)
        for pair in range(count+1):
            circuit.compose(euler(locals_[2*pair],simplify=False),[0],inplace=True)
            circuit.compose(euler(locals_[2*pair+1],simplify=False),[1],inplace=True)
            if pair<count:circuit.cx(0,1)
        actual=Operator(circuit).data
        # Correct the circuit's phase itself; the error metric never discards phase.
        circuit.global_phase+=float(np.angle(np.trace(actual.conj().T@matrix)))
        if np.max(abs(Operator(circuit).data-matrix))<2e-11:
            return circuit
    raise AssertionError('Exact two-qubit synthesis did not preserve its matrix')

def compile_exact(qc):
    expanded=QuantumCircuit(qc.num_qubits)
    expanded.global_phase=qc.global_phase
    for instruction in qc.data:
        qargs=[qc.find_bit(q).index for q in instruction.qubits]
        if instruction.operation.name=='unitary' and len(qargs)==2:
            expanded.compose(exact_two_qubit_synthesis(instruction.operation.to_matrix()),qargs,inplace=True)
        else:
            expanded.append(instruction.operation,qargs)
    compiled=transpile(expanded,basis_gates=['rz','sx','x','cx'],optimization_level=1,seed_transpiler=7)
    np.testing.assert_allclose(Operator(compiled).data,Operator(qc).data,atol=1e-10)
    return compiled
