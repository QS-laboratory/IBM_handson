from pathlib import Path
import sys,json
import numpy as np
from scipy.linalg import expm
from qiskit import qpy
from qiskit.circuit.library import StatePreparation
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector,Operator
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from circuit_building import fit_preparation,fit_evolution,evolution_circuit,compile_basic,metrics


def build(name):
    model=np.load(ROOT/f'data/reference/{name}.npz',allow_pickle=False)
    n=int(np.log2(len(model['singlet'])))
    # Ancilla q=n (highest bit); system q=0..n-1.
    target=np.concatenate([model['singlet'],model['triplet']])/np.sqrt(2)
    exact=QuantumCircuit(n+1)
    exact.append(StatePreparation(target),range(n+1))
    exact=compile_basic(exact)
    print(name,'fit preparation...',flush=True)
    layers=5 if name=='ethylene' else 12
    approx,parameters,prep_info=fit_preparation(target,layers,maxiter=2200)
    approx=compile_basic(approx)
    np.save(ROOT/f'data/circuits/{name}_prep_parameters.npy',parameters)
    print(name,'prep fidelity',prep_info['fidelity'],flush=True)
    dt=.5 if name=='ethylene' else .15
    labels=model['pauli_labels'];coeffs=model['pauli_coefficients']
    reference=evolution_circuit(labels,coeffs,np.full((6,2),dt/6))
    print(name,'fit evolution...',flush=True)
    evolved,times,evol_info=fit_evolution(model['hamiltonian'],labels,coeffs,dt,layers=1)
    reference=compile_basic(reference);evolved=compile_basic(evolved)
    target_u=expm(-1j*model['hamiltonian']*dt)
    # Independently compare Qiskit operators with the matrix reference.
    evol_info['relative_operator_error']=float(np.linalg.norm(Operator(evolved).data-target_u)/np.sqrt(len(target_u)))
    evol_info['reference_relative_operator_error']=float(np.linalg.norm(Operator(reference).data-target_u)/np.sqrt(len(target_u)))
    for qc in (reference,evolved):
        u=Operator(qc).data
        assert np.linalg.norm(u.conj().T@u-np.eye(len(u)))<1e-8
    fidelity=abs(np.vdot(target,Statevector.from_instruction(approx).data))**2
    assert abs(fidelity-prep_info['fidelity'])<1e-8
    circuits={'prep_reference':exact,'prep_compact':approx,'evol_reference':reference,'evol_compact':evolved}
    for key,qc in circuits.items():
        with (ROOT/f'data/circuits/{name}_{key}.qpy').open('wb') as f:qpy.dump(qc,f)
    meta={'molecule':name,'dt_atomic_units':dt,'evolution_parameters':times.tolist(),
          'preparation':prep_info,'evolution':evol_info,'circuits':{key:metrics(qc) for key,qc in circuits.items()},
          'method':'Dense real RY/CX state fitting; grouped Pauli-product time optimization. Independent teaching implementation, not Kanno tensor-network sweeps.',
          'hardware_validation':'NOT RUN; gate counts here assume all-to-all connectivity and rz/sx/x/cx basis; actual ISA costs depend on selected backend.',
          'reference_evolution':'6-step first-order product formula; matrix exponential is separate ideal reference',
          'version_order':'Qiskit native little endian; ancilla is highest-index qubit'}
    (ROOT/f'data/circuits/{name}.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(name,json.dumps({'prep':prep_info,'evolution':evol_info,'counts':meta['circuits']}),flush=True)

if __name__=='__main__':
    for name in sys.argv[1:] or ['ethylene','butadiene']:build(name)
