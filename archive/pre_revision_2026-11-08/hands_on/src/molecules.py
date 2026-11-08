"""Molecular reference generation. Only instructors need PySCF.
Arrays and circuits use Qiskit's native little-endian convention throughout v2.
"""
from pathlib import Path
import json
import numpy as np
from scipy.linalg import eigh

HARTREE_TO_EV = 27.211386245988


def geometry(name):
    """Planar teaching geometries in Angstrom, not optimized geometries."""
    if name == 'ethylene':
        carbons = np.array([[0., 0., 0.], [1.34, 0., 0.]])
        directions = [(0, 120), (0, 240), (1, 60), (1, -60)]
        bonds = [(0, 1, 2)]
    elif name == 'butadiene':
        carbons = np.array([[0., 0., 0.], [1.34, 0., 0.],
                            [2.07, 1.46*np.sqrt(3)/2, 0.], [3.41, 1.46*np.sqrt(3)/2, 0.]])
        directions = [(0, 120), (0, 240), (1, -60), (2, 120), (3, 60), (3, -60)]
        bonds = [(0, 1, 2), (1, 2, 1), (2, 3, 2)]
    else:
        raise ValueError(name)
    atoms = [('C', xyz.tolist()) for xyz in carbons]
    for index, angle in directions:
        a = np.deg2rad(angle)
        atoms.append(('H', (carbons[index]+1.09*np.array([np.cos(a), np.sin(a), 0.])).tolist()))
    bonds += [(i, len(carbons)+j, 1) for j, (i, _) in enumerate(directions)]
    return atoms, bonds


def fermionic_matrix(terms, nqubits):
    """Explicit Jordan-Wigner occupation-basis action with fermionic parity signs."""
    dim=1<<nqubits
    columns=np.arange(dim)
    parity=np.array([[(-1)**int(i & ((1<<q)-1)).bit_count() for i in columns] for q in range(nqubits)])
    matrix=np.zeros((dim,dim),complex)
    for label, coefficient in terms.items():
        rows=columns.copy(); valid=np.ones(dim,dtype=bool); phase=np.ones(dim)
        for token in reversed(label.split()):
            creation=token[0]=='+'
            q=int(token[2:])
            occupied=(rows>>q)&1
            valid &= occupied==int(not creation)
            phase *= parity[q,rows]
            rows ^= 1<<q
        matrix[rows[valid],columns[valid]] += coefficient*phase[valid]
    return matrix


def build_reference(name, output):
    from pyscf import gto, scf, mcscf, ao2mo, fci, lib
    from qiskit.quantum_info import SparsePauliOp, Operator
    lib.num_threads(2)
    atoms, bonds = geometry(name)
    mol = gto.M(atom=atoms, basis='sto-3g', unit='Angstrom', spin=0, charge=0, verbose=0)
    mf = scf.RHF(mol).run(conv_tol=1e-11)
    if not mf.converged:
        raise RuntimeError('RHF failed')
    norb = sum(atom[0]=='C' for atom in atoms)
    nelec = norb
    # All carbon 2pz orbitals span the isolated pi space for this planar minimal basis.
    pz = [i for i, label in enumerate(mol.ao_labels()) if 'C 2pz' in label]
    overlap = mf.get_ovlp()
    pi_projector = overlap[:, pz] @ np.linalg.inv(overlap[np.ix_(pz, pz)]) @ overlap[pz, :]
    weights = np.einsum('pi,pq,qi->i', mf.mo_coeff, pi_projector, mf.mo_coeff)
    active = sorted(np.argsort(weights)[-norb:].tolist())
    assert np.min(weights[active]) > .99, weights[active]
    mc = mcscf.CASCI(mf, norb, nelec)
    mo = mc.sort_mo(active, base=0)
    h1, ecore = mc.get_h1eff(mo)
    eri = ao2mo.restore(1, mc.get_h2eff(mo), norb)
    terms = {}
    for p in range(norb):
        for q in range(norb):
            for spin in (0, 1):
                if abs(h1[p,q]) > 1e-12:
                    terms[f'+_{p+spin*norb} -_{q+spin*norb}'] = float(h1[p,q])
            for r in range(norb):
                for s in range(norb):
                    value = .5*eri[p,q,r,s]
                    if abs(value) < 1e-12: continue
                    for spin in (0, 1):
                        for other in (0, 1):
                            label=f'+_{p+spin*norb} +_{r+other*norb} -_{s+other*norb} -_{q+spin*norb}'
                            terms[label] = terms.get(label, 0.)+float(value)
    full_h = fermionic_matrix(terms, 2*norb).real
    pauli = SparsePauliOp.from_operator(Operator(full_h), atol=1e-10).simplify(atol=1e-10)
    assert np.allclose(full_h, full_h.T, atol=1e-10)
    # Resolve spin explicitly in the fixed Nalpha=Nbeta sector.
    plus = fermionic_matrix({f'+_{p} -_{p+norb}':1. for p in range(norb)}, 2*norb)
    sz = np.diag([.5*((i & ((1<<norb)-1)).bit_count()-(i>>norb).bit_count()) for i in range(1<<(2*norb))])
    s2 = sz@sz + .5*(plus@plus.conj().T + plus.conj().T@plus)
    sector = [i for i in range(len(full_h)) if (i & ((1<<norb)-1)).bit_count()==nelec//2 and (i>>norb).bit_count()==nelec//2]
    hs = full_h[np.ix_(sector,sector)]
    ss = s2[np.ix_(sector,sector)].real
    spin_values, spin_vectors = eigh(ss)
    states, electronic_energies, spin_checks = [], [], []
    for spin in (0, 1):
        basis = spin_vectors[:,np.isclose(spin_values,spin*(spin+1))]
        e,v = eigh(basis.T@hs@basis)
        vec = np.zeros(len(full_h))
        vec[sector] = basis@v[:,0]
        # Stable global sign for reproducible real-state compilation.
        if vec[np.argmax(abs(vec))]<0: vec=-vec
        states.append(vec)
        electronic_energies.append(float(e[0]))
        spin_checks.append(float(vec@s2.real@vec))
        assert np.linalg.norm(full_h@vec-e[0]*vec)<1e-8
    singlet_pyscf, _ = fci.direct_spin0.kernel(h1,eri,norb,(nelec//2,nelec//2),ecore=ecore)
    triplet_pyscf, trip_ci = fci.direct_spin1.kernel(h1,eri,norb,(nelec//2+1,nelec//2-1),ecore=ecore)
    trip_spin = fci.spin_op.spin_square(trip_ci,norb,(nelec//2+1,nelec//2-1))[0]
    np.testing.assert_allclose(np.array(electronic_energies)+ecore,[singlet_pyscf,triplet_pyscf],atol=1e-8)
    assert abs(trip_spin-2)<1e-7
    constant = float(np.trace(full_h)/len(full_h))
    h = full_h-constant*np.eye(len(full_h))
    labels,coeffs = [],[]
    for label, coeff in pauli.to_list():
        if set(label)=={'I'}: continue
        assert abs(complex(coeff).imag)<1e-10
        labels.append(label); coeffs.append(float(complex(coeff).real))
    directory=Path(output)
    directory.mkdir(parents=True,exist_ok=True)
    xs=np.linspace(-1.5, max(a[1][0] for a in atoms)+1.5, 100)
    ys=np.linspace(-1.8, max(a[1][1] for a in atoms)+1.8, 75)
    xx,yy=np.meshgrid(xs,ys)
    # AO evaluation expects Bohr; plot plane 0.6 Angstrom above the molecular plane.
    coords=np.column_stack([xx.ravel(),yy.ravel(),np.full(xx.size,.6)])/lib.param.BOHR
    orbital_values=(mol.eval_gto('GTOval_sph',coords)@mo[:,mc.ncore:mc.ncore+norb]).reshape(len(ys),len(xs),norb)
    np.savez_compressed(directory/f'{name}.npz', hamiltonian=h, singlet=states[0], triplet=states[1],
                        orbital_x=xs,orbital_y=ys,orbital_values=orbital_values,
                        h1=h1, eri=eri, orbital_coefficients=mo[:,mc.ncore:mc.ncore+norb],
                        mo_energies=mf.mo_energy, pi_weights=weights, active_indices=np.array(active),
                        pauli_labels=np.array(labels), pauli_coefficients=np.array(coeffs))
    meta={'name':name,'formula': 'C2H4' if name=='ethylene' else 'C4H6',
          'basis':'STO-3G','geometry_source':'fixed idealized planar teaching geometry; not optimized',
          'atoms_angstrom':atoms,'bonds':bonds,'charge':0,'active_electrons':nelec,'active_spatial_orbitals':norb,
          'frozen_doubly_occupied_orbitals':int(mc.ncore),'active_mo_indices_zero_based':active,
          'qubit_mapping':'Jordan-Wigner; q0..q(norb-1)=alpha, q(norb)..q(2*norb-1)=beta',
          'array_order':'Qiskit native little endian; bit q is (index >> q) & 1',
          'system_qubits':2*norb,'total_qpde_qubits':2*norb+1,
          'singlet_energy_hartree':float(singlet_pyscf),'triplet_energy_hartree':float(triplet_pyscf),
          'gap_hartree':float(triplet_pyscf-singlet_pyscf),'gap_ev':float((triplet_pyscf-singlet_pyscf)*HARTREE_TO_EV),
          's_squared':spin_checks,'nuclear_and_frozen_core_energy_hartree':float(ecore),
          'removed_active_identity_shift_hartree':constant,'pauli_terms_without_identity':len(labels),
          'reference_validation':'JW fixed-particle S=0/1 eigenvalues match independent PySCF FCI solvers within 1e-8 Ha',
          'scope':'Vertical S0-T1 gap of this fixed-geometry CASCI model; not optical absorption energy or experimental prediction'}
    (directory/f'{name}.json').write_text(json.dumps(meta,indent=2,ensure_ascii=False)+'\n')
    print(name, 'gap/eV',meta['gap_ev'],'spin',spin_checks, 'terms',len(labels),flush=True)
    return meta
