"""Small, dense reference implementations for teaching MPS/MPO and circuits.
Tensor order follows the circuit drawing from top to bottom. No cloud jobs.
"""
from functools import reduce
import numpy as np
from scipy.linalg import expm
from scipy.optimize import least_squares
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator, Statevector
from qiskit.circuit.library import StatePreparation, UnitaryGate

I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)


def kron(*ops):
    return reduce(np.kron, ops)


def circuit_matrix(qc):
    """Convert Qiskit's little endian ordering to top-to-bottom site order."""
    return Operator(qc.reverse_bits()).data


def circuit_state(qc):
    return Statevector.from_instruction(qc.reverse_bits()).data


def tt_svd(tensor, max_bond=None, tol=1e-12):
    """Tensor train SVD; each core has (left bond, physical index, right bond)."""
    dims = tensor.shape
    remainder = np.asarray(tensor, dtype=complex)
    left = 1
    cores, spectra = [], []
    for physical in dims[:-1]:
        u, s, vh = np.linalg.svd(remainder.reshape(left * physical, -1), full_matrices=False)
        spectra.append(s.copy())
        rank = max(1, int(np.count_nonzero(s > tol)))
        if max_bond is not None:
            rank = min(rank, max_bond)
        cores.append(u[:, :rank].reshape(left, physical, rank))
        remainder = s[:rank, None] * vh[:rank]
        left = rank
    cores.append(remainder.reshape(left, dims[-1], 1))
    return cores, spectra


def contract_tt(cores):
    result = cores[0][0]
    for core in cores[1:]:
        result = np.tensordot(result, core, axes=(-1, 0))
    return result[..., 0]


def state_to_mps(psi, n, max_bond=None):
    return tt_svd(np.asarray(psi).reshape((2,) * n), max_bond)


def mps_to_state(cores):
    return contract_tt(cores).reshape(-1)


def operator_to_mpo(op, n, max_bond=None):
    # Matrix axes: out0,out1,...,in0,in1,... -> (out0,in0),(out1,in1),...
    axes = [axis for site in range(n) for axis in (site, n + site)]
    paired = op.reshape((2,) * (2 * n)).transpose(axes).reshape((4,) * n)
    cores, spectra = tt_svd(paired, max_bond)
    return [c.reshape(c.shape[0], 2, 2, c.shape[-1]) for c in cores], spectra


def mpo_to_operator(cores):
    n = len(cores)
    flat = [c.reshape(c.shape[0], 4, c.shape[-1]) for c in cores]
    paired = contract_tt(flat).reshape((2,) * (2 * n))
    axes = list(range(0, 2*n, 2)) + list(range(1, 2*n, 2))
    return paired.transpose(axes).reshape(2**n, 2**n)


def ghz_mps(theta):
    """Explicit bond-2 MPS of cos(theta)|000> + sin(theta)|111>."""
    a = np.zeros((1, 2, 2), complex)
    b = np.zeros((2, 2, 2), complex)
    c = np.zeros((2, 2, 1), complex)
    a[0, 0, 0], a[0, 1, 1] = np.cos(theta), np.sin(theta)
    b[0, 0, 0], b[1, 1, 1] = 1, 1
    c[0, 0, 0], c[1, 1, 0] = 1, 1
    return [a, b, c]


def ghz_circuit(theta):
    qc = QuantumCircuit(3)
    qc.ry(2 * theta, 0)
    qc.cx(0, 1)
    qc.cx(1, 2)
    return qc


def zz_mpo(angle):
    """exp(-i angle Z⊗Z) = cos(angle) I⊗I - i sin(angle) Z⊗Z."""
    a = np.zeros((1, 2, 2, 2), complex)
    b = np.zeros((2, 2, 2, 1), complex)
    a[0, :, :, 0] = np.cos(angle) * I
    a[0, :, :, 1] = -1j * np.sin(angle) * Z
    b[0, :, :, 0], b[1, :, :, 0] = I, Z
    return [a, b]


def zz_circuit(angle):
    qc = QuantumCircuit(2)
    qc.cx(0, 1)
    qc.rz(2 * angle, 1)
    qc.cx(0, 1)
    return qc


def ising_hamiltonian(n=3, coupling=1.0, field=0.7):
    h = np.zeros((2**n, 2**n), complex)
    for site in range(n):
        ops = [I] * n
        ops[site] = X
        h += field * kron(*ops)
    for site in range(n-1):
        ops = [I] * n
        ops[site] = ops[site+1] = Z
        h += coupling * kron(*ops)
    return h


def evolution_ansatz(parameters, n=3, layers=1):
    """Local RX gates and nearest-neighbor RZZ gates; unitary for all parameters."""
    parameters = np.asarray(parameters).reshape(layers, 2*n-1)
    qc = QuantumCircuit(n)
    for layer in parameters:
        for site in range(n):
            qc.rx(layer[site], site)
        for site in range(n-1):
            qc.rzz(layer[n+site], site, site+1)
    return qc


def trotter_parameters(t, steps, n=3, coupling=1.0, field=0.7):
    return np.tile([2*field*t/steps]*n + [2*coupling*t/steps]*(n-1), (steps, 1)).ravel()


def fit_evolution(target, t, layers, n=3, seed=7):
    """Dense educational optimization, NOT the paper's tensor-network sweeps."""
    initial = trotter_parameters(t, layers, n)
    def residual(p):
        delta = (circuit_matrix(evolution_ansatz(p, n, layers))-target).ravel()
        return np.concatenate((delta.real, delta.imag))
    rng = np.random.default_rng(seed)
    # Keep the best of the baseline and two deterministic local searches.
    candidates = [initial]
    for start in (initial, initial + rng.normal(0, 0.1, len(initial))):
        fit = least_squares(residual, start, max_nfev=250, ftol=1e-10, xtol=1e-10, gtol=1e-10)
        candidates.append(fit.x)
    best = min(candidates, key=lambda p: np.linalg.norm(residual(p)))
    return evolution_ansatz(best, n, layers)


def relative_error(approx, reference):
    return np.linalg.norm(approx-reference) / np.linalg.norm(reference)


def unitarity_error(op):
    return np.linalg.norm(op.conj().T @ op - np.eye(len(op))) / np.sqrt(len(op))


def qpde_circuit(ground, excited, evolution, t, epsilon):
    """Ancilla at the top; dense exact preparation only for this tiny example."""
    n = int(round(np.log2(len(ground))))
    joint = np.concatenate((ground, excited)) / np.sqrt(2)
    prep = StatePreparation(joint)
    qc = QuantumCircuit(n+1)
    # Reversed qargs make amplitudes correspond to top-to-bottom tensor axes.
    qc.append(prep, list(reversed(qc.qubits)))
    qc.append(UnitaryGate(evolution, label='U(t)'), list(reversed(qc.qubits[1:])))
    qc.p(epsilon*t, 0)
    qc.append(prep.inverse(), list(reversed(qc.qubits)))
    return qc


def verify():
    for theta in (0, 0.3, np.pi/4):
        psi = mps_to_state(ghz_mps(theta))
        np.testing.assert_allclose(psi, circuit_state(ghz_circuit(theta)), atol=1e-12)
        mps, _ = state_to_mps(psi, 3)
        np.testing.assert_allclose(mps_to_state(mps), psi, atol=1e-12)
    for angle in (0, 0.23, 0.7):
        u = expm(-1j*angle*kron(Z, Z))
        np.testing.assert_allclose(mpo_to_operator(zz_mpo(angle)), u, atol=1e-12)
        np.testing.assert_allclose(circuit_matrix(zz_circuit(angle)), u, atol=1e-12)
    h = ising_hamiltonian()
    t = 0.4
    u = expm(-1j*h*t)
    mpo, _ = operator_to_mpo(u, 3)
    np.testing.assert_allclose(mpo_to_operator(mpo), u, atol=1e-12)
    energies, states = np.linalg.eigh(h)
    gap = energies[1]-energies[0]
    for eps in (gap, gap+0.5, gap-0.8):
        qc = qpde_circuit(states[:, 0], states[:, 1], u, t, eps)
        np.testing.assert_allclose(abs(circuit_state(qc)[0])**2,
                                   (1+np.cos((gap-eps)*t))/2, atol=1e-10)
    rng = np.random.default_rng(19)
    for n in (2, 3, 4):
        psi = rng.normal(size=2**n)+1j*rng.normal(size=2**n)
        psi /= np.linalg.norm(psi)
        mps, _ = state_to_mps(psi, n)
        np.testing.assert_allclose(mps_to_state(mps), psi, atol=1e-12)
        op = expm(-1j*ising_hamiltonian(n)*0.2)
        mpo, _ = operator_to_mpo(op, n)
        np.testing.assert_allclose(mpo_to_operator(mpo), op, atol=1e-12)
    print('Verified: MPS/MPO reconstruction, circuit ordering, ZZ implementation, QPDE likelihood.')


def preparation_ansatz(parameters, n=3, layers=3):
    """RY and CNOT circuit for real states in this teaching example."""
    qc = QuantumCircuit(n)
    for layer in np.asarray(parameters).reshape(layers, n):
        for site in range(n):
            qc.ry(layer[site], site)
        for site in range(n-1):
            qc.cx(site, site+1)
    return qc


def fit_preparation(target, layers=3, seed=11):
    n = int(round(np.log2(len(target))))
    rng = np.random.default_rng(seed)
    def residual(p):
        delta = circuit_state(preparation_ansatz(p, n, layers))-target
        return np.concatenate((delta.real, delta.imag))
    fits = [least_squares(residual, rng.normal(0, 1, n*layers), max_nfev=400)
            for _ in range(2)]
    best = min(fits, key=lambda fit: np.linalg.norm(fit.fun))
    return preparation_ansatz(best.x, n, layers)


def apply_mpo_to_mps(mpo, mps):
    """Contract local input legs; output bond dimension is product of input bonds."""
    if len(mpo) != len(mps):
        raise ValueError('MPO and MPS must have the same number of sites')
    result = []
    for w, a in zip(mpo, mps):
        joined = np.einsum('loir,aib->laorb', w, a)
        result.append(joined.reshape(w.shape[0]*a.shape[0], 2, w.shape[-1]*a.shape[-1]))
    return result


if __name__ == '__main__':
    verify()
