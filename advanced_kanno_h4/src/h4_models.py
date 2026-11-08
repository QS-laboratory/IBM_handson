"""H4 chemistry following Kanno's molecular workflow: RHF, Boys, GA, interleaved JW.
No dependency on or modification to the existing teaching material.
"""
from pathlib import Path
import sys,json,pickle,random,itertools
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'.python_deps'))
import numpy as np
from pyscf import gto,scf,lo,ao2mo,fci,lib
from deap import base,creator,tools,algorithms
from qiskit_nature.second_q.hamiltonians import ElectronicEnergy
from qiskit_nature.second_q.mappers import InterleavedQubitMapper,JordanWignerMapper


def structure(index=0):
    """Reconstruction of Fig S6 interpolation; endpoint spacings/edges 2 Angstrom."""
    linear=np.array([[0.,3.,0.],[0.,1.,0.],[0.,-1.,0.],[0.,-3.,0.]])
    square=np.array([[-1.,1.,0.],[1.,1.,0.],[1.,-1.,0.],[-1.,-1.,0.]])
    a=index/4
    return (1-a)*linear+a*square


def exchange_cost(order,eri):
    positions=np.argsort(order)
    return float(sum(abs(eri[i,j,j,i])*(positions[i]-positions[j])**2 for i in range(4) for j in range(i)))


def ga_order(eri,seed=1):
    # Same DEAP operators, population and probabilities as the author molecular script.
    random.seed(seed);np.random.seed(seed)
    if not hasattr(creator,'KannoH4FitnessMin'):creator.create('KannoH4FitnessMin',base.Fitness,weights=(-1.,))
    if not hasattr(creator,'KannoH4Individual'):creator.create('KannoH4Individual',list,fitness=creator.KannoH4FitnessMin)
    toolbox=base.Toolbox()
    toolbox.register('indices',np.random.permutation,4)
    toolbox.register('individual',tools.initIterate,creator.KannoH4Individual,toolbox.indices)
    toolbox.register('population',tools.initRepeat,list,toolbox.individual)
    toolbox.register('mate',tools.cxOrdered)
    toolbox.register('mutate',tools.mutShuffleIndexes,indpb=.05)
    toolbox.register('select',tools.selTournament,tournsize=3)
    toolbox.register('evaluate',lambda individual:(exchange_cost(individual,eri),))
    population=toolbox.population(n=50)
    algorithms.eaSimple(population,toolbox,cxpb=.7,mutpb=.2,ngen=100,verbose=False)
    best=list(map(int,tools.selBest(population,1)[0]))
    exact=min(exchange_cost(p,eri) for p in itertools.permutations(range(4)))
    return best,exact


def build(index=0,seed=1):
    lib.num_threads(2)
    coords=structure(index)
    mol=gto.M(atom=[('H',xyz) for xyz in coords],basis='sto-3g',unit='Angstrom',spin=0,charge=0,verbose=0)
    mf=scf.RHF(mol).run(conv_tol=1e-12,max_cycle=150)
    if not mf.converged:raise RuntimeError('RHF did not converge')
    canonical=mf.mo_coeff
    boys=lo.Boys(mol,canonical).kernel()
    eri_boys=ao2mo.restore(1,ao2mo.kernel(mol,boys),4)
    order,exact_min=ga_order(eri_boys,seed)
    variants={'canonical':canonical,'boys':boys,'boys_ordered':boys[:,order]}
    output=ROOT/'data'/f'geometry_{index}'
    output.mkdir(parents=True,exist_ok=True)
    mapper=InterleavedQubitMapper(JordanWignerMapper())
    energies=[]
    for name,coeff in variants.items():
        h1=coeff.T@mf.get_hcore()@coeff
        eri=ao2mo.restore(1,ao2mo.kernel(mol,coeff),4)
        pauli=mapper.map(ElectronicEnergy.from_raw_integrals(h1,eri).second_q_op()).simplify(atol=1e-12)
        h=pauli.to_matrix().real
        eigvals,eigvecs=np.linalg.eigh(h)
        s0,_=fci.direct_spin0.kernel(h1,eri,4,(2,2),tol=1e-12)
        t1,trip_ci=fci.direct_spin1.kernel(h1,eri,4,(3,1),tol=1e-12)
        spin=fci.spin_op.spin_square(trip_ci,4,(3,1))[0]
        # The original DMRG is unconstrained: verify that its global first gap is the intended neutral gap.
        n=np.array([i.bit_count() for i in range(256)])
        neutral=[float(np.vdot(eigvecs[:,i],n*eigvecs[:,i]).real) for i in range(2)]
        if not np.allclose(neutral,[4,4],atol=1e-7):raise RuntimeError(f'Lowest full-Fock states are not neutral: {neutral}')
        if not np.isclose(eigvals[0],s0,atol=1e-7,rtol=0) or not np.isclose(eigvals[1],t1,atol=1e-7,rtol=0):
            raise RuntimeError(f'First excited state is not the desired triplet: {eigvals[:2]}, {s0}, {t1}')
        meta={'geometry_id':index,'orbital_representation':name,'basis':'STO-3G','active_electrons':4,'spatial_orbitals':4,
              'system_qubits':8,'qpde_qubits':9,'coordinates_angstrom':coords.tolist(),'geometry_source':'reconstructed linear interpolation from SI Fig S6; XYZ not published with author repository',
              'mapping':'InterleavedQubitMapper(JordanWignerMapper), alpha/beta alternating; author molecular convention',
              'orbital_order':order if name=='boys_ordered' else list(range(4)),
              'exchange_cost':exchange_cost(list(range(4)),eri),
              'ga_exhaustive_min_cost':exact_min,'ga_seed':seed,'ga_settings':{'population':50,'generations':100,'crossover':.7,'mutation':.2,'shuffle_indpb':.05},
              'singlet_electronic_energy_hartree':float(s0),'triplet_electronic_energy_hartree':float(t1),
              'nuclear_energy_hartree':float(mol.energy_nuc()),'gap_exact_hartree':float(t1-s0),'triplet_s_squared':float(spin),
              'global_lowest_states_particle_numbers':neutral,'pauli_terms':len(pauli),
              'mapped_qubit_gap_hartree':float(eigvals[1]-eigvals[0]),
              'mapping_energy_difference_vs_fci_hartree':[float(eigvals[0]-s0),float(eigvals[1]-t1)],
              'mapping_threshold_note':'Qiskit Nature mapper internally simplifies intermediate Pauli operators at its default tolerance; tiny ~1e-8 Ha differences may remain.',
              'overlap_orthonormality_error':float(np.linalg.norm(coeff.T@mf.get_ovlp()@coeff-np.eye(4))),
              'geometry_optimized':False,'frozen_core':False}
        (output/f'{name}.json').write_text(json.dumps(meta,indent=2)+'\n')
        with (output/f'{name}_pso.pickle').open('wb') as stream:pickle.dump(pauli,stream)
        (output/f'{name}.pickle').write_bytes((output/f'{name}_pso.pickle').read_bytes())
        np.savez_compressed(output/f'{name}.npz',hamiltonian=h,h1=h1,eri=eri,orbital_coefficients=coeff,
                            ao_overlap=mf.get_ovlp(),full_eigenvalues=eigvals,lowest_state=eigvecs[:,0],
                            pauli_labels=np.array(pauli.paulis.to_labels()),pauli_coefficients=pauli.coeffs,
                            rhf_orbital_energies=mf.mo_energy)
        energies.append((s0,t1))
        print('H4',index,name,'gap',t1-s0,'cost',meta['exchange_cost'],'terms',len(pauli),flush=True)
    np.testing.assert_allclose(np.array(energies),np.tile(energies[0],(3,1)),atol=1e-9,rtol=0)
    return variants

if __name__=='__main__':
    for index in map(int,sys.argv[1:] or [0,4]):build(index)
