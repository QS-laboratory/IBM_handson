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


