"""Validate the public snapshot without reading any private research files."""
from pathlib import Path
import hashlib, json
import numpy as np
from scipy.linalg import expm
from qiskit import QuantumCircuit,qpy
from qiskit.quantum_info import Operator,Statevector

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/h4_pec_v1'

def validate():
    manifest=json.loads((DATA/'manifest.json').read_text())
    for filename,digest in manifest['files'].items():
        assert hashlib.sha256((DATA/filename).read_bytes()).hexdigest()==digest,filename
    errors=dict(state=0.,operator=0.,signal=0.,eigenstate=0.)
    circuits=0;signals=0
    assert len(manifest['geometries'])==20
    for geometry in manifest['geometries']:
        directory=DATA/geometry['id']
        with np.load(directory/'model.npz',allow_pickle=False) as saved:
            model={k:saved[k].copy() for k in saved.files}
        cat=json.loads((directory/'catalog.json').read_text())
        def load(key):
            nonlocal circuits
            with (directory/geometry['assets'][key]).open('rb') as stream:q=qpy.load(stream)[0]
            circuits+=1
            return q
        h=model['hamiltonian'];exact=expm(-1j*h*.1)
        for state in ['S0','T']:
            error=float(np.max(abs(h@model[state]-float(model['energy_'+state])*model[state])))
            assert error<1e-10;errors['eigenstate']=max(errors['eigenstate'],error)
        prepared={}
        for pm,info in cat['preparation'].items():
            p0,pt,a=[load(info['keys'][k]) for k in ['S0','T','A']]
            states=[Statevector.from_instruction(p).data for p in [p0,pt]]
            for name,psi in zip(['S0','T'],states):
                np.testing.assert_allclose(abs(np.vdot(model[name],psi))**2,info['fidelities'][name],atol=1e-10)
            error=float(np.max(abs(Statevector.from_instruction(a).data-np.r_[states[0],states[1]]/np.sqrt(2))))
            assert error<1e-10;errors['state']=max(errors['state'],error)
            prepared[pm]=(states,a)
        for em,info in cat['evolution'].items():
            evo=load(info['key']);u=Operator(evo).data
            fid=float(abs(np.trace(exact.conj().T@u)/256)**2)
            np.testing.assert_allclose(fid,info['operator_fidelity'],atol=1e-10)
            if em=='exact':
                error=float(np.max(abs(u-exact)));assert error<1e-12;errors['operator']=max(errors['operator'],error)
            np.testing.assert_allclose(u.conj().T@u,np.eye(256),atol=1e-10)
            for pm,(states,a) in prepared.items():
                for L in [1,3]:
                    trial=float(model['energy_T']-model['energy_S0'])+.17
                    power=np.linalg.matrix_power(u,L)
                    x,y=[np.vdot(v,power@v) for v in states]
                    expected=float(abs(x+np.exp(1j*trial*L*.1)*y)**2/4)
                    qc=QuantumCircuit(9);qc.compose(a,inplace=True)
                    for _ in range(L):qc.compose(evo,range(8),inplace=True)
                    qc.p(trial*L*.1,8);qc.compose(a.inverse(),inplace=True)
                    actual=float(Statevector.from_instruction(qc).probabilities()[0])
                    error=abs(actual-expected);assert error<1e-9
                    errors['signal']=max(errors['signal'],error);signals+=1
                    if L==1:
                        fringe=cat['pairs'][f'{pm}/{em}']['fringe']
                        np.testing.assert_allclose(-np.angle(x.conjugate()*y)/.1,fringe['apparent_gap'],atol=1e-10)
        print('Validated',geometry['id'],flush=True)
    report=dict(status='passed',snapshot_sha256=hashlib.sha256((DATA/'manifest.json').read_bytes()).hexdigest(),
                geometries=20,qpy_circuits=circuits,signal_checks=signals,max_errors=errors,qpu_used=False)
    (ROOT/'data/h4_pec_v1_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)
    return report

if __name__=='__main__':validate()
