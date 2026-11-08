"""Check the public teaching examples, formula/circuit agreement and mode controls."""
from pathlib import Path
from types import SimpleNamespace
import re
import nbformat
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]


def cells(name):
    return nbformat.read(ROOT/name,4).cells

@pytest.fixture(scope='module')
def matrix():
    import matplotlib.pyplot as plt
    import IPython.display
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(plt,'show',lambda:None)
        patch.setattr(IPython.display,'display',lambda *args,**kwargs:None)
        ns={}
        for cell in cells('01_matrix_phase_estimation.ipynb'):
            if cell.cell_type=='code':exec(cell.source,ns)
        plt.close('all')
        return ns


def test_matrix_qpe_and_bayesian_examples(matrix):
    np.testing.assert_array_equal(matrix['PHASES'],[.25,.75])
    for m in [2,3,4]:
        for psi,phase in zip([matrix['zero'],matrix['one']],matrix['PHASES']):
            p=matrix['qpe_distribution'](psi,m)
            np.testing.assert_allclose(p[int(phase*2**m)],1,atol=1e-12)
    for name,truth in [('post_exact',.25),('post_gap',.5),('post_c',.5)]:
        posterior=matrix[name]
        np.testing.assert_allclose(posterior.sum(axis=1),1,atol=1e-12)
        assert abs(matrix['GRID'][np.argmax(posterior[-1])]-truth)<.001
    np.testing.assert_allclose(matrix['post_gap'],matrix['post_c'],atol=1e-12)
    for L in [1,3,8]:
        for beta in [0,.3,-1.2]:
            b=matrix['measurement_probabilities'](matrix['return_circuit'](L=L,beta=beta,kind='B'))
            c=matrix['measurement_probabilities'](matrix['return_circuit'](L=L,beta=beta,kind='C'))
            expected=(1+np.cos(2*np.pi*L*.5+beta))/2
            np.testing.assert_allclose([b[0],b[1],c[1]],expected,atol=1e-12)
            np.testing.assert_allclose(c[0],1,atol=1e-12)


def test_all_explicit_questions_have_answers():
    for path in ROOT.glob('0*.ipynb'):
        for index,cell in enumerate(nbformat.read(path,4).cells):
            if cell.cell_type!='markdown':continue
            prompt=re.search(r'(?m)^#{2,6} .*?(?:課題|演習[0-9０-９]|確認問題)|\*\*(?:実行前の)?予想',cell.source)
            if prompt:
                assert ('<summary>解答例</summary>' in cell.source
                        or '### 解答例' in cell.source),(path.name,index)
    text='\n'.join(c.source for c in cells('01_matrix_phase_estimation.ipynb') if c.cell_type=='markdown')
    assert '|001' not in text and 'P(001)' not in text
    assert '## 3．iterative QPE' not in text


@pytest.fixture(scope='module')
def h4():
    ns={}
    exec(cells('04_h4_circuit_assessment.ipynb')[2].source,ns)
    return ns


def test_h4_same_assets_and_formula_for_both_notebooks(h4):
    other={};exec(cells('05_h4_iterative_bayes.ipynb')[2].source,other)
    assert h4['ASSET_DIGEST']==other['ASSET_DIGEST']
    for case in h4['CASES']:
        for L in [1,9,64]:
            states,u=h4['CASES'][case]
            power=np.linalg.matrix_power(u,L)
            a,b=[np.vdot(psi,power@psi) for psi in states]
            eps=np.array([-.3,h4['GAP'],.4])
            expected=abs(a+np.exp(1j*eps*L*h4['DT'])*b)**2/4
            np.testing.assert_allclose(h4['probability'](case,L,eps),expected,atol=2e-11)
            np.testing.assert_allclose(other['probability'](case,L,eps),expected,atol=2e-11)


def notebook_run(ns,name,mode,case,tmp_path,monkeypatch,shots=128,hardware=False):
    notebook_cells=cells(name)
    config=next(c.source for c in notebook_cells if c.cell_type=='code' and c.source.startswith("MODE = 'exact'"))
    config=config.replace("MODE = 'exact'",f"MODE = '{mode}'")
    config=re.sub(r"CASE = '(?:exact|both)'",f"CASE = '{case}'",config)
    config=config.replace('SHOTS = 1024',f'SHOTS = {shots}')
    if mode=='ibm':
        config=config.replace("BACKEND_NAME = ''","BACKEND_NAME = 'local_test'")
        config=config.replace('RUN_HARDWARE = False',f'RUN_HARDWARE = {hardware}')
        config=config.replace('MAX_HARDWARE_L = 8','MAX_HARDWARE_L = 4')
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(ns['plt'],'show',lambda:None)
    exec(config,ns)
    for cell in notebook_cells:
        if cell.cell_type!='code':continue
        if ('measured=meter.measure(CASE' in cell.source or 'history=iterative_bayes' in cell.source
                or cell.source.startswith('if history:') or cell.source.startswith('ROUNDS =')
                or (name.startswith('05') and cell.source.startswith((
                    'plan1=choose_measurement_settings','FIRST_READY=',
                    'if FIRST_READY:', 'if FIRST_FIT_OK:')))):
            exec(cell.source,ns)
    return ns['meter']

@pytest.mark.parametrize('name',['04_h4_circuit_assessment.ipynb','05_h4_iterative_bayes.ipynb'])
@pytest.mark.parametrize('mode',['exact','shots'])
@pytest.mark.parametrize('case',['exact','prep_only','evolution_only','both'])
def test_notebook_simulation_mode_switches(h4,name,mode,case,tmp_path,monkeypatch):
    def no_connection(*args,**kwargs):raise AssertionError('Simulation attempted IBM connection')
    monkeypatch.setitem(h4,'connect_ibm',no_connection)
    meter=notebook_run(h4,name,mode,case,tmp_path,monkeypatch)
    assert len(meter.records)==(1 if name.startswith('04') else 5)
    for record in meter.records:
        assert record['case']==case
        if mode=='exact':assert record['counts_all_zero'] is None and record['shots'] is None
        else:
            assert record['shots']==128
            np.testing.assert_allclose(np.asarray(record['values'])*128,record['counts_all_zero'])
    # Recreating the settings cell resets the RNG as documented.
    first=[record['values'] for record in meter.records]
    again=notebook_run(h4,name,mode,case,tmp_path,monkeypatch)
    assert first==[record['values'] for record in again.records]

@pytest.mark.parametrize('name',['04_h4_circuit_assessment.ipynb','05_h4_iterative_bayes.ipynb'])
def test_notebook_ibm_branch_with_local_sampler(h4,name,tmp_path,monkeypatch):
    """Exercise the actual public notebook control flow without IBM credentials/QPU."""
    import qiskit_ibm_runtime
    from qiskit.primitives import StatevectorSampler
    jobs={};submissions=[];previews=[]
    class Sampler:
        def __init__(self,mode):pass
        def run(self,circuits,shots):
            submissions.append((len(circuits),shots))
            job=StatevectorSampler(seed=42).run(circuits,shots=shots)
            jobs[job.job_id()]=job
            return job
    backend=SimpleNamespace(name='local_test')
    service=SimpleNamespace(job=lambda job_id:jobs[job_id])
    monkeypatch.setattr(qiskit_ibm_runtime,'SamplerV2',Sampler)
    monkeypatch.setitem(h4,'connect_ibm',lambda *args:(service,backend))
    monkeypatch.setitem(h4,'hardware_preview',lambda backend,L:previews.append(L))
    # Counts come from the actual nine-qubit circuit; real ISA compilation is checked separately.
    monkeypatch.setitem(h4,'transpile',lambda circuits,**kwargs:circuits)
    preview=notebook_run(h4,name,'ibm','both',tmp_path,monkeypatch,hardware=False)
    assert not preview.records and not submissions
    measured=notebook_run(h4,name,'ibm','both',tmp_path,monkeypatch,hardware=True)
    assert len(measured.records)==1  # E5: first L=4; next L is above the limit and stops gracefully.
    assert measured.records[0]['source']=='ibm_raw_counts'
    assert submissions==[(9 if name.startswith('04') else 21,128)]
    cached=notebook_run(h4,name,'ibm','both',tmp_path,monkeypatch,hardware=True)
    assert len(cached.records)==1 and len(submissions)==1
    if name.startswith('05'):
        assert len(h4['history'])==1 and h4['history'][0]['L']==4
        assert (tmp_path/'outputs/exercise5_updates.json').exists()


def test_invalid_hardware_case_rejected_before_connection(h4,tmp_path,monkeypatch):
    monkeypatch.setitem(h4,'connect_ibm',lambda *args:pytest.fail('Invalid case reached authentication'))
    with pytest.raises(ValueError,match="CASE='both'"):
        notebook_run(h4,'05_h4_iterative_bayes.ipynb','ibm','exact',tmp_path,monkeypatch)


def test_hardware_preview_compiles_actual_circuit(h4):
    from qiskit.providers.fake_provider import GenericBackendV2
    backend=GenericBackendV2(9,basis_gates=['rz','sx','x','cx'],
        coupling_map=[[i,j] for i in range(9) for j in range(9) if i!=j],seed=7)
    report=h4['hardware_preview'](backend,L=1)
    assert report['two_qubit_gates']>0 and report['depth']>0


def test_h2_encoding_and_trotter_accuracy(tmp_path,monkeypatch):
    import IPython.display
    import matplotlib.pyplot as plt
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(IPython.display,'display',lambda *args,**kwargs:None)
    monkeypatch.setattr(plt,'show',lambda:None)
    h2={};jw={}
    n=cells('02_h2_one_qubit.ipynb')
    for c in n:
        if c.cell_type=='code' and ('import sys' in c.source or c.source.startswith('H = np.array')):exec(c.source,h2)
    n=cells('03_h2_jw_bpde.ipynb')
    for c in n:
        if c.cell_type=='code' and (c.source.startswith('"""Small dense') or c.source.startswith("H=DATA['hamiltonian']")):exec(c.source,jw)
    np.testing.assert_allclose(h2['E'],[jw['ENERGIES'][k] for k in ['S0','S2']],atol=1e-12)
    np.testing.assert_allclose(h2['ENUC'],jw['ENUC'],atol=1e-12)
    # The general C-readout expression must agree with the implemented 5-qubit circuit,
    # also when evolution is approximated and inputs are no longer eigenstates.
    from qiskit import QuantumCircuit
    from qiskit.circuit.library import PauliEvolutionGate
    from qiskit.synthesis import SuzukiTrotter
    errors=[]
    for slices in [1,4,16]:
        qc=QuantumCircuit(4);qc.append(PauliEvolutionGate(jw['PAULI'],time=jw['DT'],
            synthesis=SuzukiTrotter(order=2,reps=slices)),range(4))
        # Operator(PauliEvolutionGate) evaluates the exact exponential; decompose
        # first so that we assess the selected Trotter synthesis, as the notebook does.
        ut=jw['Operator'](qc.decompose()).data
        errors.append(np.linalg.norm(ut-jw['U']))
        power=np.linalg.matrix_power(ut,32)
        for target in ['S1','S2','T']:
            a,b=[np.vdot(jw['STATES'][j],power@jw['STATES'][j]) for j in ['S0',target]]
            for beta in [0,.3]:
                actual=jw['bpde_probability'](target,jw['theta'],ut,32,beta,True)
                np.testing.assert_allclose(actual,abs(a+np.exp(1j*beta)*b)**2/4,atol=1e-11)
    assert errors[0]>errors[1]>errors[2]
    plt.close('all')


@pytest.fixture(scope='module')
def h2_complete():
    import matplotlib.pyplot as plt
    import IPython.display
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(plt,'show',lambda:None)
        patch.setattr(IPython.display,'display',lambda *args,**kwargs:None)
        ns={}
        for cell in cells('02_h2_one_qubit.ipynb'):
            if cell.cell_type=='code':exec(cell.source,ns)
        plt.close('all')
        return ns


def test_h2_spin_and_configuration_block_against_jw(h2_complete):
    """Independently derive spin from fermionic operators, not the teaching 4x4 S²."""
    import importlib.util
    path=ROOT.parent/'self_study/molecular_jw/tools/jw_core.py'
    spec=importlib.util.spec_from_file_location('h2_jw_spin_check',path)
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    number,sz,s2=core.spin_operators()
    indices=[5,10,6,9]
    ns=h2_complete
    np.testing.assert_allclose(s2[np.ix_(indices,indices)],ns['S2_MS0'],atol=1e-12)
    with np.load(ROOT.parent/'self_study/molecular_jw/data/h2.npz') as data:
        full=data['hamiltonian']
    np.testing.assert_allclose(full[np.ix_(indices,indices)],ns['H_MS0'],atol=1e-12)
    for name,psi in ns['STATES_MS0'].items():
        fullpsi=np.zeros(16);fullpsi[indices]=psi
        energy=ns['ENERGIES_MS0'][name]
        np.testing.assert_allclose(full@fullpsi,energy*fullpsi,atol=1e-12)
        np.testing.assert_allclose(number@fullpsi,2*fullpsi,atol=1e-12)
        np.testing.assert_allclose(sz@fullpsi,0,atol=1e-12)
        np.testing.assert_allclose(s2@fullpsi,(2 if name=='T' else 0)*fullpsi,atol=1e-12)
    # Same occupations, different relative phase and spin.
    np.testing.assert_allclose(abs(ns['STATES_MS0']['S1'])**2,abs(ns['STATES_MS0']['T'])**2,atol=1e-12)
    assert abs(np.vdot(ns['STATES_MS0']['S1'],ns['STATES_MS0']['T']))<1e-12


def test_h2_exact_encoding_and_molecular_bayesian_results(h2_complete):
    from scipy.linalg import expm
    ns=h2_complete;v=ns['V']
    for L in [1,3,32]:
        full=expm(-1j*ns['H_MS0']*L*ns['DT'])
        logical=ns['Operator'](ns['evolution_gate'](L)).data
        np.testing.assert_allclose(full@v,v@logical,atol=1e-12)
        np.testing.assert_allclose((np.eye(4)-v@v.T)@full@v,0,atol=1e-12)
    np.testing.assert_allclose(ns['Statevector'].from_instruction(ns['prep_g']).data,ns['g'],atol=1e-12)
    np.testing.assert_allclose(ns['Statevector'].from_instruction(ns['prep_e']).data,ns['e'],atol=1e-12)
    for key,truth in [('post_exact',ns['PHASES'][0]),('post_gap',ns['GAP']),('post_c',ns['GAP'])]:
        p=ns[key]
        np.testing.assert_allclose(p.sum(axis=1),1,atol=1e-12)
        estimate=ns['GRID'][np.argmax(p[-1])]
        assert abs(estimate-truth)<.001
    for key in ['post_gap','post_c']:
        estimate=ns['GRID'][np.argmax(ns[key][-1])]
        assert abs(ns['phase_to_gap'](estimate)-ns['ENERGY_GAP'])<.005


@pytest.fixture(scope='module')
def jw_complete():
    import matplotlib.pyplot as plt
    import IPython.display
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(plt,'show',lambda:None)
        patch.setattr(IPython.display,'display',lambda *args,**kwargs:None)
        ns={}
        for cell in cells('03_h2_jw_bpde.ipynb'):
            if cell.cell_type=='code':exec(cell.source,ns)
        plt.close('all')
        return ns


def test_manual_pauli_rotations_include_y_and_global_phase(jw_complete):
    from scipy.linalg import expm
    ns=jw_complete
    for label in ['IIII','XXYY','YIYX','ZYIY','IIIZ','IZZI','IYII']:
        for coefficient,time in [(.37,.23),(-.13,.1)]:
            qc=ns['pauli_rotation'](label,coefficient,time)
            expected=expm(-1j*coefficient*time*ns['SparsePauliOp'](label).to_matrix())
            np.testing.assert_allclose(ns['Operator'](qc).data,expected,atol=1e-12)


def test_jw_trotter_product_and_physical_bpde(jw_complete):
    ns=jw_complete
    for order in [1,2]:
        errors=[]
        for slices in [1,4,16]:
            actual=ns['Operator'](ns['build_trotter'](slices,order)).data
            np.testing.assert_allclose(actual,ns['dense_product'](slices,order),atol=1e-11)
            errors.append(np.linalg.norm(actual-ns['U']))
        assert errors[0]>errors[1]>errors[2]
    for target in ['S1','S2','T']:
        np.testing.assert_allclose(ns['STATES'][target],ns['EXPECTED_STATES'][target],atol=1e-12)
        for steps in [1,3,17]:
            power=np.linalg.matrix_power(ns['UT'],steps)
            a0,aj=[np.vdot(ns['STATES'][k],power@ns['STATES'][k]) for k in ['S0',target]]
            for gap in [ns['GAPS'][target],ns['GAPS'][target]+.2]:
                qc=ns['physical_bpde'](target,steps,gap)
                actual=ns['Statevector'].from_instruction(qc).probabilities()[0]
                expected=abs(a0+np.exp(1j*gap*steps*ns['DT'])*aj)**2/4
                np.testing.assert_allclose(actual,expected,atol=1e-11)


def test_h2_bridge_answers_are_visible_and_numbered():
    text=cells('02_h2_one_qubit.ipynb')[-1].source
    assert '### 解答例（課題9：演習3への橋渡し）' in text
    question,answer=text.split('### 解答例',1)
    for index in [1,2,3]:
        assert re.search(rf'(?m)^{index}\. ',question)
        assert re.search(rf'(?m)^{index}\. ',answer)


@pytest.fixture(scope='module')
def h4_assessed(tmp_path_factory):
    import matplotlib.pyplot as plt
    with pytest.MonkeyPatch.context() as patch:
        patch.chdir(tmp_path_factory.mktemp('h4-assessed'))
        patch.setattr(plt,'show',lambda:None)
        ns={}
        for cell in cells('04_h4_circuit_assessment.ipynb'):
            if cell.cell_type=='code':exec(cell.source,ns)
        plt.close('all')
        return ns


def test_h4_resources_and_tensor_fit_are_distinct(h4_assessed):
    ns=h4_assessed;costs=ns['COSTS']
    np.testing.assert_allclose(ns['pauli'].to_matrix(),ns['H'],atol=1e-12,rtol=0)
    assert any(1e-12<abs(c)<1e-5 for c in ns['pauli'].coeffs)
    assert costs['Approx S0 preparation']['cx']==60
    assert costs['Approx T preparation']['cx']==60
    assert costs['Approx coherent A']['cx']==258
    assert costs['Approx evolution with conversion']['cx']==216
    assert costs['M6 prep + M16 evolution L=1']['cx']==732
    assert costs['Exact coherent A']['cx']>costs['Approx coherent A']['cx']
    assert costs['Trotter step']['cx']>costs['Approx evolution with conversion']['cx']
    assert costs['Exact prep + Trotter L=1']['cx']>costs['M6 prep + M16 evolution L=1']['cx']
    assert ns['PROVENANCE']['results']['fidelities']['S0_MPS']>1-1e-12
    assert .93<ns['STATE_FIDELITIES']['S0']<.95
    assert ns['STATE_FIDELITIES']['T']>.999
    assert ns['OPERATOR_FIDELITY']>.99999
    assert ns['FIT_INFO']['full_convergence_claimed'] is False
    native=ns['COMPILED_CIRCUITS']['M6 prep + M16 evolution L=1']
    actual=ns['Statevector'].from_instruction(native).probabilities()[0]
    np.testing.assert_allclose(actual,ns['probability']('both',1,ns['GAP']),atol=1e-10)


def test_h4_l1_curves_match_physical_readout_in_all_cases(h4_assessed):
    ns=h4_assessed
    from qiskit import QuantumCircuit
    for case,(states,u) in ns['CASES'].items():
        a=ns['A_APPROX'] if case in ['prep_only','both'] else ns['A_EXACT']
        for gap in [ns['GAP']-.05,ns['GAP'],ns['GAP']+.05]:
            qc=QuantumCircuit(9);qc.compose(a,inplace=True)
            qc.unitary(u,range(8));qc.p(gap*ns['DT'],8);qc.compose(a.inverse(),inplace=True)
            actual=ns['Statevector'].from_instruction(qc).probabilities()[0]
            np.testing.assert_allclose(actual,ns['probability'](case,1,gap),atol=1e-10)
    metrics=ns['L1_METRICS']
    np.testing.assert_allclose(metrics['exact']['peak'],ns['GAP'],atol=1e-10)
    np.testing.assert_allclose(metrics['exact']['pmax'],1,atol=1e-10)
    assert metrics['both']['pmax']>.9999
    assert -.007<metrics['both']['bias']<-.006


def test_h4_bridge_contains_explicit_visible_answers():
    bridge=cells('04_h4_circuit_assessment.ipynb')[-1].source
    assert '### 解答例（課題7：演習5への接続）' in bridge
    assert 'Gaussian' in bridge and '二項尤度' in bridge
    assert '同じH、同じMの選択回路' in bridge


@pytest.fixture(scope='module')
def h4_bayes_complete(tmp_path_factory):
    import matplotlib.pyplot as plt
    work=tmp_path_factory.mktemp('h4-bayes-complete')
    with pytest.MonkeyPatch.context() as patch:
        patch.chdir(work);patch.setattr(plt,'show',lambda:None)
        ns={}
        for cell in cells('05_h4_iterative_bayes.ipynb'):
            if cell.cell_type=='code':exec(cell.source,ns)
        plt.close('all')
        ns['_workdir']=work
        return ns


def test_bayes_walkthrough_matches_first_loop_round(h4_bayes_complete):
    ns=h4_bayes_complete;first=ns['history'][0]
    assert ns['FIRST_FIT_OK'] and len(ns['history'])==5
    np.testing.assert_allclose(first['observed'],ns['observed_1'],atol=0)
    np.testing.assert_allclose([first['posterior_mu'],first['posterior_variance']],
                              [ns['mu1'],ns['v1']],atol=1e-12)
    assert ns['history'][1]['L']==ns['plan2']['L']
    assert ns['history'][1]['prior_mu']==ns['mu1']
    assert ns['history'][1]['prior_variance']==ns['v1']
    for row in ns['history']:
        residual=np.asarray(row['observed'])-np.asarray(row['fitted_probability'])
        np.testing.assert_allclose(row['residual_rmse'],np.sqrt(np.mean(residual**2)),atol=1e-15)
        assert row['source']=='exact_probabilities' and row['shots'] is None


def test_bias_controls_preserve_primary_data_and_use_no_shots(h4_bayes_complete):
    import json
    ns=h4_bayes_complete
    primary=json.loads((ns['_workdir']/'outputs/exercise5_updates.json').read_text())
    assert primary==ns['history']
    assert set(ns['comparison'])=={'exact','prep_only','evolution_only','both'}
    for case,rows in ns['comparison'].items():
        assert len(rows)>=6
        assert all(row['mode']=='exact' and row['case']==case for row in rows)
        assert all(record['source']=='exact_probabilities' and record['shots'] is None
                   for record in ns['comparison_meters'][case].records)
    # With exact preparation, the saved evolution circuit gives a repeatable phase bias.
    evo=ns['comparison']['evolution_only'][-1]
    assert .001<evo['posterior_mu']-ns['GAP']<.003


def test_round_figure_exposes_plan_fit_residual_and_update(h4_bayes_complete,monkeypatch):
    ns=h4_bayes_complete;captured=[]
    monkeypatch.setattr(ns['plt'],'show',lambda:captured.append(ns['plt'].gcf()))
    ns['plot_bayes_round'](ns['history'][0],'exact','exact')
    assert len(captured)==1 and len(captured[0].axes)==4
    for prefix,ax in zip(['A.','B.','C.','D.'],captured[0].axes):assert ax.get_title().startswith(prefix)


def test_h2_phase_energy_sign_and_readout_match_jw_at_same_time(h2_complete,jw_complete):
    """02 uses DT=1; 03 uses DT=.1. Match physical time, not the integer L."""
    one,jw=h2_complete,jw_complete
    np.testing.assert_allclose(one['ENERGY_GAP'],jw['GAPS']['S2'],atol=1e-12)
    for l_one in [1,3]:
        l_jw=10*l_one;t=l_one*one['DT']
        np.testing.assert_allclose(t,l_jw*jw['DT'],atol=1e-12)
        for energy_candidate in [0.,one['ENERGY_GAP'],one['ENERGY_GAP']+.3]:
            phase_candidate=(-energy_candidate*one['DT']/(2*np.pi))%1
            phase_beta=-2*np.pi*l_one*phase_candidate
            energy_beta=energy_candidate*t
            np.testing.assert_allclose(np.exp(1j*phase_beta),np.exp(1j*energy_beta),atol=1e-12)
            expected=(1+np.cos((energy_candidate-one['ENERGY_GAP'])*t))/2
            small=one['readout'](one['return_circuit'](L=l_one,beta=energy_beta,kind='C'))[1]
            full=jw['bpde_probability']('S2',jw['theta'],jw['U'],l_jw,energy_beta,True)
            np.testing.assert_allclose([small,full],expected,atol=1e-11)


def test_c_ancilla_is_constant_only_for_exact_preparation(matrix,h2_complete):
    for ns,readout in [(matrix,'measurement_probabilities'),(h2_complete,'readout')]:
        for eta in [0.,.25]:
            actual=ns[readout](ns['return_circuit'](eta=eta,L=1,beta=.31,kind='C'))[0]
            expected=1-eta*(1-eta)*(1-np.cos(2*np.pi*ns['GAP']))
            np.testing.assert_allclose(actual,expected,atol=1e-12)
            if eta==0:np.testing.assert_allclose(actual,1,atol=1e-12)
            else:assert actual<.99
