from pathlib import Path
from types import SimpleNamespace
import sys,json
import numpy as np
import pytest
from scipy.linalg import expm
from qiskit.quantum_info import Statevector,Operator
from qiskit import QuantumCircuit
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import workshop as w
from circuit_building import append_pauli,single_ry,cx_permutation,preparation_circuit,prep_operations


@pytest.mark.parametrize('name',['ethylene','butadiene'])
def test_molecular_references(name):
    meta,a=w.load_model(name)
    h=a['hamiltonian']
    assert np.allclose(h,h.T.conj())
    energies=[]
    n=meta['active_spatial_orbitals']
    indices=np.arange(len(h))
    alpha=np.array([(int(i)&((1<<n)-1)).bit_count() for i in indices])
    beta=np.array([(int(i)>>n).bit_count() for i in indices])
    for label in ['singlet','triplet']:
        state=a[label]
        assert np.isclose(np.linalg.norm(state),1)
        energy=np.vdot(state,h@state).real
        energies.append(energy)
        assert np.linalg.norm(h@state-energy*state)<1e-8
        assert np.sum(abs(state[(alpha!=n//2)|(beta!=n//2)])**2)<1e-10
    assert abs(energies[1]-energies[0]-meta['gap_hartree'])<1e-9
    assert np.allclose(meta['s_squared'],[0,2])
    target=np.concatenate([a['singlet'],a['triplet']])/np.sqrt(2)
    fidelity=abs(np.vdot(target,Statevector.from_instruction(w.load_circuit(name,'prep_compact')).data))**2
    assert fidelity>.99
    summary=w.circuit_summary(name)
    assert summary['circuits']['prep_compact']['cx_count']<summary['circuits']['prep_reference']['cx_count']
    assert summary['circuits']['evol_compact']['cx_count']<summary['circuits']['evol_reference']['cx_count']


@pytest.mark.parametrize('label',['XY','YZ','ZI','YY','IX'])
def test_pauli_gate_implementation(label):
    from qiskit.quantum_info import Pauli
    qc=QuantumCircuit(len(label));append_pauli(qc,label,.27)
    target=expm(-1j*.27*Pauli(label).to_matrix())
    np.testing.assert_allclose(Operator(qc).data,target,atol=1e-12)


def test_native_state_order_matches_optimizer():
    n=4;layers=3;p=np.random.default_rng(14).normal(size=n*layers)
    state=np.zeros(1<<n);state[0]=1
    for kind,a,b in prep_operations(n,layers):
        state=single_ry(state,p[b],a) if kind=='ry' else state[cx_permutation(n,a,b)]
    np.testing.assert_allclose(state,Statevector.from_instruction(preparation_circuit(p,n,layers)).data,atol=1e-12)


@pytest.mark.parametrize('variant',['reference','compact'])
def test_final_circuit_preserves_distributed_blocks_and_likelihood(variant):
    for eps in [0.,.175,.5]:
        circuit=w.qpde_circuit('ethylene',eps,2,variant,measure=False)
        actual=abs(Statevector.from_instruction(circuit).data[0])**2
        expected=w.zero_probability('ethylene',eps,2,variant)
        assert abs(actual-expected)<1e-10
    meta,a=w.load_model('ethylene')
    # Independent dense exact preparation/evolution interference formula.
    t=1.3;eps=.2
    u=expm(-1j*a['hamiltonian']*t)
    amp=.5*(np.vdot(a['singlet'],u@a['singlet'])+np.exp(1j*eps*t)*np.vdot(a['triplet'],u@a['triplet']))
    assert abs(abs(amp)**2-(1+np.cos((meta['gap_hartree']-eps)*t))/2)<1e-10


def test_binomial_data_and_sequential_update():
    data=w.read_json(ROOT/'data/example_results/ethylene_compact_depolarized.json')
    history=w.infer_gap(data)
    assert len(history)==3
    for row in history:assert np.isclose(row['posterior'].sum(),1)
    meta,_=w.load_model('ethylene')
    lo,hi=history[-1]['ci95_hartree']
    assert lo<meta['gap_hartree']<hi
    assert history[-1]['ci95_hartree'][1]-history[-1]['ci95_hartree'][0]<.04
    assert w.inference_diagnostics(data,history)['detectable_scan_variation']
    bad=json.loads(json.dumps(data));bad['source']='mitigated_quasiprobabilities'
    with pytest.raises(ValueError):w.infer_gap(bad)
    bad=json.loads(json.dumps(data));bad['rows'][0]['successes']=-1
    with pytest.raises(ValueError):w.validate_counts(bad)
    flat=json.loads(json.dumps(data))
    for row in flat['rows']:row['successes']=0
    assert not w.inference_diagnostics(flat,w.infer_gap(flat))['detectable_scan_variation']


def test_job_submission_recovery_and_duplicate_guard(tmp_path,monkeypatch):
    from qiskit.providers.fake_provider import GenericBackendV2
    import qiskit_ibm_runtime
    backend=GenericBackendV2(num_qubits=5,seed=11)
    jobs={}
    class Counts:
        def get_counts(self):return {'00000':70,'00001':30}
    class Job:
        def __init__(self,id,pubs):self.id=id;self.pubs=pubs
        def job_id(self):return self.id
        def status(self):return 'DONE'
        def result(self):return [SimpleNamespace(data=SimpleNamespace(meas=Counts())) for _ in self.pubs]
    class Sampler:
        def __init__(self,mode,options):pass
        def run(self,pubs,shots):
            id=f'fake-{len(jobs)}';job=Job(id,pubs);jobs[id]=job;return job
    class Service:
        def backend(self,name):return backend
        def job(self,id):return jobs[id]
    monkeypatch.setattr(qiskit_ibm_runtime,'SamplerV2',Sampler)
    plan=w.make_plan(shots=100,step_counts=[1]);plan['rows']=plan['rows'][:2]
    path=tmp_path/'record.json'
    record=w.submit_ibm(plan,Service(),'fake_backend',path)
    assert path.exists() and len(jobs)==1
    assert w.read_json(path)['jobs'][0]['job_id']=='fake-0'
    recovered=w.fetch_ibm(w.read_json(path),Service())
    assert len(recovered['rows'])==2 and recovered['rows'][0]['successes']==70
    with pytest.raises(FileExistsError):w.submit_ibm(plan,Service(),'fake_backend',path)
    assert len(jobs)==1


def test_partial_job_record_can_be_analyzed(monkeypatch):
    plan=w.make_plan(shots=100,step_counts=[1,2])
    record=json.loads(json.dumps(plan))
    record['source']='ibm_submission'
    record['jobs']=[{'job_id':'partial','row_indices':[0,1], 'iteration':1}]
    class Counts:
        def get_counts(self):return {'00000':50,'00001':50}
    class Job:
        def status(self):return 'DONE'
        def result(self):return [SimpleNamespace(data=SimpleNamespace(meas=Counts())) for _ in range(2)]
    class Service:
        def job(self,id):return Job()
    recovered=w.fetch_ibm(record,Service())
    assert recovered['partial_plan']
    assert len(recovered['rows'])==2
    assert recovered['requested_row_count']==18
    assert len(w.infer_gap(recovered))==1
