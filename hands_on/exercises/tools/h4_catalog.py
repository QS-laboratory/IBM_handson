"""Bundle compatible saved H4 circuits; no optimization or external submission."""
import base64
import hashlib
import io
import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from qiskit import QuantumCircuit, qpy, transpile
from qiskit.circuit.library import UnitaryGate, CXGate
from qiskit.synthesis import TwoQubitBasisDecomposer, TwoQubitWeylDecomposition
from qiskit.quantum_info import Operator, Statevector

BASIS=['rz','sx','x','cx']


def build_catalog(study, assets, provenance):
    prep=study/'h4_mps_preparation_study/results/pec_20/R_2p00000000'
    evo=study/'h4_mpo_orbital_study/corrected_mpo/pec_evolution_20/R_2p00000000'
    source_hashes={}
    def digest(path):
        value=hashlib.sha256(path.read_bytes()).hexdigest()
        source_hashes[str(path.relative_to(study.parent))]=value
        return value
    def loadq(path):
        digest(path)
        with path.open('rb') as stream:return qpy.load(stream)[0]
    def encode(name,qc):
        stream=io.BytesIO();qpy.dump(qc,stream)
        assets[name]=base64.b64encode(stream.getvalue()).decode()
        return name
    def compileq(qc):
        return transpile(qc,basis_gates=BASIS,optimization_level=1,seed_transpiler=7,approximation_degree=1)
    def exact_local(g):
        # Preserve small Weyl components that automatic specialization may drop.
        w=TwoQubitWeylDecomposition(g,fidelity=None)
        parts=TwoQubitBasisDecomposer(CXGate()).decomp3_supercontrolled(w._inner_decomposition)
        q=QuantumCircuit(2)
        for k in range(4):
            q.unitary(parts[2*k],[0]);q.unitary(parts[2*k+1],[1])
            if k<3:q.cx(0,1)
        q.global_phase+=np.angle(np.vdot(Operator(q).data,g))
        np.testing.assert_allclose(Operator(q).data,g,atol=1e-10)
        return q
    def native_local(g, controlled=False):
        raw=QuantumCircuit(3 if controlled else 2)
        if controlled:raw.append(UnitaryGate(g).control(1),[0,1,2])
        else:raw.unitary(g,[0,1])
        expected=Operator(raw).data
        native=compileq(raw)
        if np.max(abs(Operator(native).data-expected))>1e-10:
            raw=QuantumCircuit(3 if controlled else 2)
            if controlled:raw.append(exact_local(g).to_gate().control(1),[0,1,2])
            else:raw.compose(exact_local(g),inplace=True)
            native=compileq(raw)
        np.testing.assert_allclose(Operator(native).data,expected,atol=1e-10,rtol=0)
        return native
    def resources(qc):
        return dict(cx=int(qc.count_ops().get('cx',0)),depth=qc.depth(),qubits=qc.num_qubits)
    def metadata(path):
        digest(path);record=json.loads(path.read_text())
        return {k:v for k,v in record.items() if not isinstance(v,(list,dict))}
    model_hash=digest(prep/'model.npz')
    with np.load(prep/'model.npz') as data:model={k:data[k].copy() for k in data.files}
    with np.load(evo/'chemistry.npz') as data:
        np.testing.assert_allclose(data['orbital_coefficients'],model['orbital_coefficients'],atol=1e-12,rtol=0)
    digest(evo/'chemistry.npz')
    exact=expm(-1j*model['hamiltonian']*.1)
    states_exact=[model['S0'],model['T']]
    catalog={'preparation':{},'evolution':{},'pairs':{},'model_sha256':model_hash}
    preps={};evolutions={}
    for m in [4,6,8]:
        report0=metadata(prep/f'M{m}.json');reportt=metadata(prep/f'T_M{m}.json')
        p0,pt=loadq(prep/f'M{m}_basic.qpy'),loadq(prep/f'T_M{m}_basic.qpy')
        psi=[Statevector.from_instruction(q).data for q in [p0,pt]]
        fid=[float(abs(np.vdot(x,y))**2) for x,y in zip(states_exact,psi)]
        np.testing.assert_allclose(fid,[report0['fidelity_S0'],reportt['fidelity_T']],atol=1e-10)
        if m==6:
            paired=qpy.load(io.BytesIO(base64.b64decode(assets['A'])))[0]
        else:
            with np.load(prep/f'M{m}_state.npz') as data:g0=data['gates'].copy()
            with np.load(prep/f'T_M{m}_state.npz') as data:gt=data['gates'].copy()
            digest(prep/f'M{m}_state.npz');digest(prep/f'T_M{m}_state.npz')
            starts=[s for layer in range(m) for s in range(layer%2,7,2)]
            assert len(g0)==len(gt)==len(starts)
            paired=QuantumCircuit(9,name=f'A_prep_M{m}');paired.h(8)
            assert report0['determinant']==reportt['determinant']==165
            for q in [0,2,5,7]:paired.x(q)
            for u0,ut,s in zip(g0,gt,starts):
                paired.compose(native_local(u0),[s,s+1],inplace=True)
                paired.compose(native_local(ut@u0.conj().T,controlled=True),[8,s,s+1],inplace=True)
            paired=compileq(paired)
        expected=np.block([[Operator(p0).data,Operator(p0).data],
                           [Operator(pt).data,-Operator(pt).data]])/np.sqrt(2)
        np.testing.assert_allclose(Operator(paired).data,expected,atol=2e-10,rtol=0)
        keys=dict(S0=encode(f'prep_{m}_S0',p0),T=encode(f'prep_{m}_T',pt),A=encode(f'prep_{m}_A',paired))
        source=dict(model=model_hash,S0=digest(prep/f'M{m}_basic.qpy'),T=digest(prep/f'T_M{m}_basic.qpy'))
        catalog['preparation'][m]=dict(keys=keys,reports={'S0':report0,'T':reportt},
            fidelities={'S0':fid[0],'T':fid[1],'S0_MPS':report0['fidelity_after_MPS'],'T_MPS':reportt['fidelity_after_MPS']},
            resources={'S0':resources(p0),'T':resources(pt),'A':resources(paired),'inverse':resources(paired.inverse())},sources=source)
        preps[m]=psi
    # Exact signed Fock permutation, interleaved -> blocked. It is a unitary change of mode order.
    desired=[0,2,4,6,1,3,5,7];current=list(range(8));bridge=QuantumCircuit(8)
    for target,mode in enumerate(desired):
        position=current.index(mode)
        while position>target:
            bridge.swap(position-1,position);bridge.cz(position-1,position)
            current[position-1],current[position]=current[position],current[position-1];position-=1
    for m in [8,16,32]:
        report=metadata(evo/f'M{m}/results.json')
        if m==16:
            qc=qpy.load(io.BytesIO(base64.b64decode(assets['evolution'])))[0]
        else:
            raw=loadq(evo/f'M{m}/basic.qpy')
            qc=QuantumCircuit(8,name=f'evolution_M{m}_blocked')
            qc.compose(bridge.inverse(),inplace=True);qc.compose(raw,inplace=True);qc.compose(bridge,inplace=True)
            qc=compileq(qc)
        reference=loadq(evo/f'M{m}/blocked_blocks.qpy')
        u=Operator(qc).data
        np.testing.assert_allclose(u,Operator(reference).data,atol=1e-10,rtol=0)
        np.testing.assert_allclose(u.conj().T@u,np.eye(256),atol=1e-10)
        fidelity=float(abs(np.trace(exact.conj().T@u)/256)**2)
        catalog['evolution'][m]=dict(key=encode(f'evolution_{m}',qc),report=report,
            resources=resources(qc),operator_fidelity=fidelity,
            sources=dict(basic=digest(evo/f'M{m}/basic.qpy'),blocked=digest(evo/f'M{m}/blocked_blocks.qpy')))
        evolutions[m]=u
    for pm,psi in preps.items():
        for em,u in evolutions.items():
            a0,at=[np.vdot(v,u@v) for v in psi]
            spec=dict(model=model_hash,preparation=catalog['preparation'][pm]['sources'],
                      evolution=catalog['evolution'][em]['sources'],preparation_M=pm,evolution_M=em,
                      construction='paired-block-v2; checked-local-synthesis; signed-Fock-v1; rz/sx/x/cx; level1; seed7')
            fingerprint=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()
            catalog['pairs'][f'{pm}/{em}']=dict(digest=fingerprint,
                fringe=dict(apparent_gap=-float(np.angle(a0.conjugate()*at))/.1,
                    Pmax=float((abs(a0)+abs(at))**2/4)))
    # The original 6/16 default is a fixed numerical regression reference.
    np.testing.assert_allclose(catalog['pairs']['6/16']['fringe']['apparent_gap'],
        provenance['results']['fringe']['apparent_gap'],atol=1e-10)
    provenance['catalog_sources']=source_hashes
    return catalog
