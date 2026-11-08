"""Spin-resolved full-CI potential curves for neutral H4/STO-3G.
All four spatial orbitals/electrons are retained. Nuclear repulsion is included.
Use symmetric orthogonalized AOs: no SCF or orbital-root-following approximation.
"""
from pathlib import Path
import csv,json
import numpy as np
from scipy.linalg import eigh
from pyscf import gto,ao2mo,fci,lib
import pyscf

ROOT=Path(__file__).resolve().parent
EV=27.211386245988
NORB=4
NELEC=(2,2)
DIM=36  # choose(4,2)^2; complete M_S=0 determinant space


def spin_bases():
    columns=[]
    for j in range(DIM):
        vector=np.eye(DIM)[:,j].reshape(6,6)
        columns.append(fci.spin_op.contract_ss(vector,NORB,NELEC).reshape(-1))
    s2=np.column_stack(columns)
    values,vectors=eigh(s2)
    singlet=vectors[:,np.isclose(values,0.,atol=1e-10)]
    triplet=vectors[:,np.isclose(values,2.,atol=1e-10)]
    assert singlet.shape==(36,20)
    assert triplet.shape==(36,15)
    return s2,{0:singlet,1:triplet}


def solve(coordinates,s2,bases):
    mol=gto.M(atom=[('H',xyz.tolist()) for xyz in coordinates],basis='sto-3g',unit='Angstrom',spin=0,charge=0,verbose=0)
    overlap=mol.intor('int1e_ovlp')
    eigenvalues,eigenvectors=eigh(overlap)
    coeff=(eigenvectors/np.sqrt(eigenvalues))@eigenvectors.T
    assert np.linalg.norm(coeff.T@overlap@coeff-np.eye(4))<1e-10
    hcore=mol.intor('int1e_kin')+mol.intor('int1e_nuc')
    h1=coeff.T@hcore@coeff
    eri=ao2mo.restore(1,ao2mo.kernel(mol,coeff),NORB)
    # Exact dense FCI Hamiltonian, independent of the quantum mapping's term thresholds.
    h2eff=fci.direct_spin1.absorb_h1e(h1,eri,NORB,NELEC,fac=.5)
    columns=[]
    for j in range(DIM):
        vector=np.eye(DIM)[:,j].reshape(6,6)
        columns.append(fci.direct_spin1.contract_2e(h2eff,vector,NORB,NELEC).reshape(-1))
    h=np.column_stack(columns)
    assert np.linalg.norm(h-h.T)<1e-10
    assert np.linalg.norm(h@s2-s2@h)<1e-10
    result={'nuclear_repulsion_hartree':float(mol.energy_nuc()),'minimum_ao_overlap_eigenvalue':float(eigenvalues.min())}
    for spin,label in [(0,'singlet'),(1,'triplet')]:
        basis=bases[spin]
        e,v=eigh(basis.T@h@basis)
        state=basis@v[:,0]
        residual=np.linalg.norm(h@state-e[0]*state)
        measured=float(state@s2@state)
        assert residual<1e-9
        assert abs(measured-spin*(spin+1))<1e-9
        result[f'{label}_electronic_hartree']=float(e[0])
        result[f'{label}_total_hartree']=float(e[0]+mol.energy_nuc())
        result[f'{label}_s_squared']=measured
        result[f'{label}_fci_residual']=float(residual)
    result['gap_hartree']=result['triplet_electronic_hartree']-result['singlet_electronic_hartree']
    result['gap_ev']=result['gap_hartree']*EV
    # Independent standard PySCF solvers, without the explicit spin-projector construction.
    e0,ci0=fci.direct_spin0.kernel(h1,eri,NORB,NELEC,tol=1e-12)
    et,cit=fci.direct_spin1.kernel(h1,eri,NORB,(3,1),tol=1e-12)
    np.testing.assert_allclose(e0,result['singlet_electronic_hartree'],atol=1e-9,rtol=0)
    np.testing.assert_allclose(et,result['triplet_electronic_hartree'],atol=1e-9,rtol=0)
    return result


def linear(distance):
    return np.array([[0.,1.5*distance,0.],[0.,.5*distance,0.],
                     [0.,-.5*distance,0.],[0.,-1.5*distance,0.]])


def save_csv(path,rows):
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)


def main():
    lib.num_threads(2)
    s2,bases=spin_bases()
    scans={}
    for scan,values,coordinates,parameter in [
        ('linear_stretch',np.round(np.linspace(.6,4.,69),8),linear,'distance_angstrom')]:
        rows=[]
        for value in values:
            row={parameter:float(value),**solve(coordinates(value),s2,bases)}
            rows.append(row)
        scans[scan]=rows
        save_csv(ROOT/f'{scan}.csv',rows)
        print(scan,len(rows),'points',flush=True)
    # Same R=2 A geometry used in the author-method exercise; compare independent FCI references.
    reference=json.loads((ROOT.parent/'data/geometry_0/canonical.json').read_text())
    at_two=next(r for r in scans['linear_stretch'] if r['distance_angstrom']==2.)
    for label in ['singlet','triplet']:
        np.testing.assert_allclose(at_two[f'{label}_electronic_hartree'],reference[f'{label}_electronic_energy_hartree'],atol=1e-9,rtol=0)
    hydrogen=gto.M(atom='H 0 0 0',basis='sto-3g',unit='Angstrom',spin=1,verbose=0)
    e_h=float((hydrogen.intor('int1e_kin')+hydrogen.intor('int1e_nuc'))[0,0]/hydrogen.intor('int1e_ovlp')[0,0])
    metadata={'method':'Full-CI in the complete 4-electron/4-spatial-orbital STO-3G space; no frozen core, no active-space truncation',
              'basis':'STO-3G','states':['lowest S=0 singlet','lowest S=1 triplet'],
              'energies_include_nuclear_repulsion':True,'geometry_optimized':False,
              'linear_geometry':'positions (0,1.5R,0),(0,0.5R,0),(0,-0.5R,0),(0,-1.5R,0), Angstrom',
              'orbital_basis':'symmetric orthogonalized complete AO basis; FCI is invariant to orbital rotations',
              'spin_resolution':'diagonalize S^2 in full Ms=0 determinant space (dimension 36), then Hamiltonian in S=0 (20) and S=1 (15) subspaces',
              'independent_validation':'standard PySCF direct_spin0 and direct_spin1(Ms=1) energies at every point; R=2 compared to the stored author-method exercise reference',
              'pyscf_version':pyscf.__version__,'hydrogen_atom_energy_hartree':e_h,'four_hydrogen_dissociation_limit_hartree':4*e_h,
              'linear_at_2_angstrom':at_two,
              'sampled_linear_minima':{label:min(scans['linear_stretch'],key=lambda r:r[f'{label}_total_hartree']) for label in ['singlet','triplet']}}
    (ROOT/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print('R=2 A:',at_two,flush=True)
    for label,row in metadata['sampled_linear_minima'].items():print('sampled minimum',label,row['distance_angstrom'],row[f'{label}_total_hartree'],flush=True)

if __name__=='__main__':main()
