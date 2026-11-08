"""Full electronic SU(2) DMRG PEC, following pyqmd2's block2 benchmark policy.
No FCI Hamiltonian, determinant enumeration, FCI seed or FCI energy is used by DMRG.
"""
from pathlib import Path
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import json,csv,tempfile,shutil,time,hashlib,sys,importlib.metadata
import numpy as np
from scipy.linalg import eigh
from pyscf import gto,ao2mo,lib
from pyblock2.driver.core import DMRGDriver,SymmetryTypes,MPOAlgorithmTypes
import block2

ROOT=Path(__file__).resolve().parent
REFERENCE_SOURCE=Path('/Users/sato/project/pyqmd2/src/pyqmd2/experimental/dmrg_benchmark.py')
SEEDS=[20261024,20261025]
BOND=16  # 4^(4//2), a conservative full-space Schmidt bound for four Electron sites
ENERGY_TOL=1e-11
DW_TOL=1e-12


def integrals(distance):
    coordinates=np.array([[0.,1.5*distance,0.],[0.,.5*distance,0.],[0.,-.5*distance,0.],[0.,-1.5*distance,0.]])
    mol=gto.M(atom=[('H',xyz.tolist()) for xyz in coordinates],basis='sto-3g',unit='Angstrom',spin=0,verbose=0)
    values,vectors=eigh(mol.intor('int1e_ovlp'))
    coeff=(vectors/np.sqrt(values))@vectors.T
    h=coeff.T@(mol.intor('int1e_kin')+mol.intor('int1e_nuc'))@coeff
    eri=ao2mo.restore(1,ao2mo.kernel(mol,coeff),4)
    return np.ascontiguousarray(h),np.ascontiguousarray(eri),float(mol.energy_nuc())


def run_state(h,v,ecore,spin,seed):
    scratch=Path(tempfile.mkdtemp(prefix='h4-pec-block2-'))
    start=time.perf_counter()
    driver=DMRGDriver(scratch=str(scratch),stack_mem=256*1024**2,symm_type=SymmetryTypes.SU2,
                      n_threads=1,n_mkl_threads=1,clean_scratch=True)
    try:
        block2.Random.rand_seed(seed)
        # block2 spin is 2*S. This is physical electronic DMRG, not seniority-zero/DOCI.
        driver.initialize_system(n_sites=4,n_elec=4,spin=2*spin,singlet_embedding=False)
        mpo=driver.get_qc_mpo(h,v,ecore=ecore,reorder=None,algo_type=MPOAlgorithmTypes.Conventional,iprint=0)
        ket=driver.get_random_mps(tag='KET',bond_dim=BOND,nroots=1,dot=2)
        driver.dmrg(mpo,ket,n_sweeps=4,tol=0.,bond_dims=[BOND],
                    noises=[1e-5,1e-5,1e-6,1e-6],thrds=[1e-12]*4,dav_max_iter=1000,iprint=0)
        energies=[float(x[0]) for x in driver._dmrg.energies]
        weights=[float(x) for x in driver._dmrg.discarded_weights]
        zero_sweeps=0;converged=False
        while zero_sweeps<12:
            energy=float(driver.dmrg(mpo,ket,n_sweeps=2,tol=0.,bond_dims=[BOND],
                                    noises=[0.,0.],thrds=[1e-14]*2,dav_max_iter=1000,iprint=0))
            step_energies=[float(x[0]) for x in driver._dmrg.energies]
            step_weights=[float(x) for x in driver._dmrg.discarded_weights]
            energies.extend(step_energies);weights.extend(step_weights);zero_sweeps+=len(step_energies)
            if not step_energies:raise RuntimeError('No DMRG sweeps returned')
            if zero_sweeps>=2 and abs(energies[-1]-energies[-2])<ENERGY_TOL and max(weights[-2:])<DW_TOL:
                converged=True;break
        spin_mpo=driver.get_spin_square_mpo(iprint=0)
        measured=float(driver.expectation(ket,spin_mpo,ket,iprint=0))
        result={'spin':spin,'seed':seed,'total_energy_hartree':energy,'s_squared':measured,
                'bond_dimension':BOND,'full_schmidt_bound':16,'sweeps':len(energies),
                'zero_noise_sweeps':zero_sweeps,'last_energy_change':abs(energies[-1]-energies[-2]),
                'last_two_max_discarded_weight':max(weights[-2:]),'converged':converged,
                'elapsed_seconds':time.perf_counter()-start,'sweep_energies_hartree':energies,
                'sweep_discarded_weights':weights}
        if not converged:raise RuntimeError(f'DMRG not converged: {result}')
        if abs(measured-spin*(spin+1))>1e-8:raise RuntimeError(f'Wrong total spin: {result}')
        return result
    finally:
        driver.finalize();shutil.rmtree(scratch,ignore_errors=True)


def main():
    lib.num_threads(1)
    reference=np.genfromtxt(ROOT/'linear_stretch.csv',delimiter=',',names=True)
    start=time.perf_counter();rows=[];reports=[]
    for point in reference:
        distance=float(point['distance_angstrom'])
        h,v,ecore=integrals(distance)
        row={'distance_angstrom':distance};details={'distance_angstrom':distance,'states':{}}
        for spin,label in [(0,'singlet'),(1,'triplet')]:
            runs=[run_state(h,v,ecore,spin,seed) for seed in SEEDS]
            energies=[r['total_energy_hartree'] for r in runs]
            spread=max(energies)-min(energies)
            if spread>ENERGY_TOL:raise RuntimeError(f'Independent seeds disagree: R={distance}, {label}, {energies}')
            energy=min(energies)
            row[f'{label}_dmrg_total_hartree']=energy
            row[f'{label}_fci_total_hartree']=float(point[f'{label}_total_hartree'])
            row[f'{label}_dmrg_minus_fci_hartree']=energy-row[f'{label}_fci_total_hartree']
            row[f'{label}_seed_spread_hartree']=spread
            row[f'{label}_s_squared']=runs[0]['s_squared']
            row[f'{label}_max_sweeps']=max(r['sweeps'] for r in runs)
            details['states'][label]=runs
            if abs(row[f'{label}_dmrg_minus_fci_hartree'])>1e-9:
                raise RuntimeError(f'Full-rank DMRG/FCI validation failed: {row}')
        rows.append(row);reports.append(details)
        print(f"R={distance:.2f} S0 error={row['singlet_dmrg_minus_fci_hartree']:.2e} T1 error={row['triplet_dmrg_minus_fci_hartree']:.2e}",flush=True)
    with (ROOT/'linear_dmrg.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    metadata={'method':'full electronic SU(2) two-site DMRG with block2; no seniority restriction',
              'reference_implementation':str(REFERENCE_SOURCE),
              'reference_sha256':hashlib.sha256(REFERENCE_SOURCE.read_bytes()).hexdigest() if REFERENCE_SOURCE.exists() else None,
              'block2_version':importlib.metadata.version('block2'),'basis':'STO-3G','electrons':4,'spatial_orbitals':4,
              'spin_targets':[0,1],'block2_spin_arguments':[0,2],'singlet_embedding':False,
              'bond_dimensions':[BOND],'schmidt_bound':16,'seeds':SEEDS,'threads':1,'memory_stack_bytes':256*1024**2,
              'warmup_noises':[1e-5,1e-5,1e-6,1e-6],'max_total_sweeps':16,'energy_tolerance':ENERGY_TOL,'discarded_weight_tolerance':DW_TOL,
              'no_fci_arrays_used_by_dmrg':True,'no_fci_initial_state_or_energy_used':True,
              'fci_comparison_is_posthoc_only':True,'points':len(rows),'dmrg_runs':len(rows)*4,
              'elapsed_seconds':time.perf_counter()-start,
              'max_abs_error_vs_fci_hartree':{label:max(abs(r[f'{label}_dmrg_minus_fci_hartree']) for r in rows) for label in ['singlet','triplet']},
              'max_seed_spread_hartree':max(r[f'{label}_seed_spread_hartree'] for r in rows for label in ['singlet','triplet']),
              'reports':reports}
    (ROOT/'dmrg_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print('SUMMARY',json.dumps({k:v for k,v in metadata.items() if k!='reports'}),flush=True)

if __name__=='__main__':
    if '--pilot' in sys.argv:
        lib.num_threads(1);h,v,e=integrals(2.)
        for spin in [0,1]:print(json.dumps(run_state(h,v,e,spin,SEEDS[0])),flush=True)
    else:main()
