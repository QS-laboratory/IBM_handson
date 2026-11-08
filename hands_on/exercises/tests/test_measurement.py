"""Numerical and hardware-adapter checks; never contact IBM."""
import json
from pathlib import Path
from types import SimpleNamespace
import nbformat
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]

@pytest.fixture(scope='module')
def runtime():
    notebook=nbformat.read(ROOT/'04_h4_circuit_assessment.ipynb',4)
    ns={}
    exec(notebook.cells[2].source,ns)
    return ns

@pytest.mark.parametrize('L',[1,3,16])
def test_repeated_circuit_matches_probability(runtime,L):
    for gap in [runtime['GAP'],runtime['GAP']+.37]:
        qc=runtime['experiment_circuit'](L,gap)
        p=runtime['Statevector'].from_instruction(qc).probabilities()[0]
        np.testing.assert_allclose(p,runtime['probability']('both',L,gap),atol=1e-10,rtol=0)


def test_exact_and_shots_are_distinct(runtime):
    cls=runtime['Measurement'];eps=np.linspace(-2,2,21)
    exact=cls('exact').measure('both',4,eps,'scan')
    np.testing.assert_allclose(exact,runtime['probability']('both',4,eps),atol=0)
    np.testing.assert_array_equal(exact,cls('exact',shots=None).measure('both',4,eps,'scan'))
    np.testing.assert_array_equal(exact,cls('exact',shots=0).measure('both',4,eps,'scan'))
    sampled=cls('shots',shots=31,seed=19).measure('both',4,eps,'scan')
    np.testing.assert_allclose(sampled*31,np.rint(sampled*31),atol=1e-12)
    np.testing.assert_array_equal(sampled,cls('shots',shots=31,seed=19).measure('both',4,eps,'scan'))
    assert not np.array_equal(sampled,exact)
    with pytest.raises(ValueError):cls('shots',shots=0)
    for invalid in [0,-1,1.5,True]:
        with pytest.raises(ValueError,match='positive integer'):
            cls('exact').measure('both',invalid,eps,'scan')
        with pytest.raises(ValueError,match='positive integer'):
            runtime['experiment_circuit'](invalid,0)


@pytest.mark.parametrize('case',['exact','both'])
def test_iterative_handoff_for_both_simulation_modes(runtime,tmp_path,monkeypatch,case):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(runtime['plt'],'show',lambda:None)
    for mode in ['exact','shots']:
        meter=runtime['Measurement'](mode,shots=4096,seed=42)
        rows=runtime['iterative_bayes'](meter,case=case,rounds=5)
        assert len(rows)==5
        for previous,current in zip(rows,rows[1:]):
            assert current['prior_mu']==previous['posterior_mu']
            assert current['prior_variance']==previous['posterior_variance']
            assert current['L']>previous['L']
        assert rows[-1]['posterior_variance']<rows[0]['prior_variance']
        assert len(meter.records)==5
        if mode=='exact':assert all(r['counts_all_zero'] is None for r in meter.records)
        else:assert all(r['shots']==4096 for r in meter.records)
    assert json.loads((tmp_path/'outputs/exercise5_updates.json').read_text())


def test_hardware_submission_recovery_and_guards(runtime,tmp_path,monkeypatch):
    import qiskit_ibm_runtime
    submitted=[];recovered=[]
    counts=[{'000000000':7,'111111111':3},{'000000000':4,'111111111':6}]
    fail={'result':True}
    class Job:
        def job_id(self):return 'test-job-42'
        def result(self):
            if fail['result']:raise TimeoutError('interrupted while waiting')
            return [SimpleNamespace(data=SimpleNamespace(meas=SimpleNamespace(get_counts=lambda c=c:c))) for c in counts]
    class Sampler:
        def __init__(self,mode):pass
        def run(self,circuits,shots):
            submitted.append((circuits,shots));return Job()
    service=SimpleNamespace(job=lambda job_id:(recovered.append(job_id) or Job()))
    backend=SimpleNamespace(name='fake_backend')
    monkeypatch.setattr(qiskit_ibm_runtime,'SamplerV2',Sampler)
    monkeypatch.setitem(runtime,'transpile',lambda circuits,**kwargs:circuits)
    meter=runtime['Measurement']('ibm',shots=10,service=service,backend=backend,
        directory=tmp_path,run_hardware=False,max_hardware_L=8)
    with pytest.raises(ValueError,match='RUN_HARDWARE'):meter.measure('both',1,[0,.1],'round_1')
    assert not submitted
    meter.run_hardware=True
    with pytest.raises(ValueError,match='MAX_HARDWARE_L'):meter.measure('both',9,[0,.1],'round_1')
    with pytest.raises(ValueError,match='CASE=both'):meter.measure('exact',1,[0,.1],'round_1')
    with pytest.raises(TimeoutError):meter.measure('both',1,[0,.1],'round_1')
    assert len(submitted)==1
    saved=json.loads(next(tmp_path.glob('*.json')).read_text());assert saved['job_id']=='test-job-42'
    fail['result']=False
    np.testing.assert_allclose(meter.measure('both',1,[0,.1],'round_1'),[.7,.4])
    assert recovered==['test-job-42'] and len(submitted)==1
    np.testing.assert_allclose(meter.measure('both',1,[0,.1],'round_1'),[.7,.4])
    assert len(submitted)==1 and recovered==['test-job-42']
    assert meter.records[-1]['source']=='ibm_raw_counts'


def test_unresolved_submission_marker_is_not_resubmitted(runtime,tmp_path,monkeypatch):
    import qiskit_ibm_runtime
    class Sampler:
        def __init__(self,mode):pass
        def run(self,*args,**kwargs):raise ConnectionError('uncertain submission')
    monkeypatch.setattr(qiskit_ibm_runtime,'SamplerV2',Sampler)
    monkeypatch.setitem(runtime,'transpile',lambda circuits,**kwargs:circuits)
    meter=runtime['Measurement']('ibm',service=object(),backend=SimpleNamespace(name='fake'),
        run_hardware=True,directory=tmp_path)
    with pytest.raises(ConnectionError):meter.measure('both',1,[0],'test')
    with pytest.raises(RuntimeError,match='Unresolved submission marker'):meter.measure('both',1,[0],'test')


def test_measurement_plan_uses_prior_variance_and_quantizes_time(runtime):
    plan=runtime['choose_measurement_settings'](0.,4.,21)
    assert plan['requested_time']==.45 and plan['L']==4 and plan['time']==.4
    np.testing.assert_allclose(plan['trial_gaps'],np.linspace(-4,4,21),atol=0)
    shifted=runtime['choose_measurement_settings'](.7,.25,5)
    np.testing.assert_allclose(shifted['trial_gaps'],np.linspace(.45,.95,5),atol=1e-15)
    assert shifted['L']==int((1.8/.25)/runtime['DT'])
    for v in [0,-1,np.inf]:
        with pytest.raises(ValueError):runtime['choose_measurement_settings'](0,v)


@pytest.mark.parametrize('mu,v,m,w',[(0,4,.02,11),(.1,.04,.12,.01),(-.3,.5,.4,.7)])
def test_gaussian_product_matches_numerically_normalized_density(runtime,mu,v,m,w):
    post_mu,post_v=runtime['update_gaussian_prior'](mu,v,m,w)
    x=np.linspace(post_mu-9*np.sqrt(post_v),post_mu+9*np.sqrt(post_v),30001)
    product=runtime['gaussian'](x,.31,mu,np.sqrt(v))*runtime['gaussian'](x,.87,m,np.sqrt(w))
    normalized=product/np.trapezoid(product,x)
    mean=np.trapezoid(normalized*x,x)
    variance=np.trapezoid(normalized*(x-mean)**2,x)
    np.testing.assert_allclose([post_mu,post_v],[mean,variance],atol=1e-12)
    assert 0<post_v<min(v,w)


def test_gaussian_fit_uses_signal_width_not_parameter_covariance(runtime):
    x=np.linspace(-.7,1.1,41)
    y=runtime['gaussian'](x,.9,.23,.3)
    result=runtime['fit_gaussian_signal'](x,y,0.,1.)
    np.testing.assert_allclose([result['amplitude'],result['fit_mu'],result['fit_variance']],
                              [.9,.23,.09],atol=1e-10)
    assert result['residual_rmse']<1e-12
    assert result['fit_variance']>.08  # A noiseless fit still has a nonzero signal width.


def test_fit_failure_keeps_observed_data_and_last_valid_posterior(runtime,tmp_path,monkeypatch):
    monkeypatch.chdir(tmp_path)
    meter=runtime['Measurement']('exact')
    original=runtime['fit_gaussian_signal'];calls=[]
    def fail_second(*args):
        calls.append(1)
        if len(calls)==2:raise RuntimeError('injected fit failure')
        return original(*args)
    monkeypatch.setitem(runtime,'fit_gaussian_signal',fail_second)
    rows=runtime['iterative_bayes'](meter,rounds=3,show_plots=False)
    assert len(rows)==1 and len(meter.records)==2
    assert len(json.loads((tmp_path/'outputs/exercise5_updates.json').read_text()))==1
    assert len(json.loads((tmp_path/'outputs/exercise5_measurements.json').read_text()))==2


def test_no_output_control_does_not_overwrite_main_history(runtime,tmp_path,monkeypatch):
    monkeypatch.chdir(tmp_path)
    main=runtime['iterative_bayes'](runtime['Measurement']('exact'),rounds=2,show_plots=False)
    path=tmp_path/'outputs/exercise5_updates.json';before=path.read_bytes()
    control=runtime['iterative_bayes'](runtime['Measurement']('exact'),case='both',rounds=3,
        show_plots=False,output_prefix=None,verbose=False)
    assert len(main)==2 and len(control)==3 and path.read_bytes()==before
