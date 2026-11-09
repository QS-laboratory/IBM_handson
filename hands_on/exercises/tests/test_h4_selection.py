"""Saved-M selection must change the circuits, signals and provenance together."""
from pathlib import Path
import nbformat
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]

@pytest.fixture(scope='module')
def runtime():
    ns={}
    exec(next(c.source for c in nbformat.read(ROOT/'04_h4_circuit_assessment.ipynb',4).cells if c.cell_type == 'code' and c.source.startswith('ASSETS=')),ns)
    return ns

@pytest.mark.parametrize('pm',[4,6,8])
@pytest.mark.parametrize('em',[8,16,32])
def test_all_saved_combinations(runtime,pm,em):
    ns=runtime
    ns['select_circuits'](pm,em)
    assert ns['SELECTED_PREP_M']==pm and ns['SELECTED_EVOLUTION_M']==em
    p0,pt=[ns['Statevector'].from_instruction(ns['load_circuit_asset'](k)).data for k in ['S0','T']]
    np.testing.assert_allclose(ns['Statevector'].from_instruction(ns['A_APPROX']).data,
                               np.r_[p0,pt]/np.sqrt(2),atol=2e-10)
    for L in [1,3]:
        for eps in [ns['GAP'],.37]:
            measured=ns['Statevector'].from_instruction(ns['experiment_circuit'](L,eps)).probabilities()[0]
            np.testing.assert_allclose(measured,ns['probability']('both',L,eps),atol=2e-10)
    for k,psi,ref in zip(['S0','T'],[p0,pt],ns['EXACT_STATES']):
        np.testing.assert_allclose(abs(np.vdot(psi,ref))**2,ns['SELECTED_RESULTS']['fidelities'][k],atol=1e-10)
    eps=np.linspace(-1,1,21)
    for mode in ['exact','shots']:
        meter=ns['Measurement'](mode,shots=2048,seed=42)
        values=meter.measure('both',1,eps,'scan')
        if mode=='exact':np.testing.assert_allclose(values,ns['probability']('both',1,eps),atol=0)
        else:np.testing.assert_allclose(values*2048,np.rint(values*2048),atol=1e-12)
        assert meter.records[-1]['prep_m']==pm and meter.records[-1]['evolution_m']==em
        history=ns['iterative_bayes'](meter,'both',rounds=2,show_plots=False,output_prefix=None,verbose=False)
        assert len(history)==2
        assert history[1]['prior_mu']==history[0]['posterior_mu']
        assert all(row['prep_m']==pm and row['evolution_m']==em for row in history)
    ns['select_circuits'](6,16)


def test_invalid_selection_is_atomic_and_old_meter_is_rejected(runtime):
    ns=runtime;ns['select_circuits'](6,16)
    meter=ns['Measurement']('exact')
    for args in [(5,16),(6,12),(True,16),(4,16.),(None,16)]:
        with pytest.raises(ValueError):ns['select_circuits'](*args)
        assert ns['SELECTED_PREP_M']==6 and ns['SELECTED_EVOLUTION_M']==16
    ns['select_circuits'](4,8)
    with pytest.raises(ValueError,match='recreate Measurement'):meter.measure('both',1,[0.],'stale')
    ns['select_circuits'](6,16)


def test_selection_controls_and_embedded_catalog_match_between_notebooks():
    notebooks=[nbformat.read(ROOT/name,4) for name in ['04_h4_circuit_assessment.ipynb','05_h4_iterative_bayes.ipynb']]
    namespaces=[]
    for book in notebooks:
        ns={};exec(next(c.source for c in book.cells if c.cell_type == 'code' and c.source.startswith('ASSETS=')),ns);namespaces.append(ns)
    assert namespaces[0]['CIRCUIT_CATALOG']==namespaces[1]['CIRCUIT_CATALOG']
    assert namespaces[0]['ASSET_DIGEST']==namespaces[1]['ASSET_DIGEST']
    for book in notebooks:
        controls=[c for c in book.cells if c.cell_type=='code' and c.source.startswith('PREP_M =')]
        assert len(controls)==1
        assert 'EVOLUTION_M = 16' in controls[0].source
        assert 'select_circuits(PREP_M, EVOLUTION_M)' in controls[0].source


def test_hardware_cache_separates_m_and_reuses_same_selection(runtime,tmp_path,monkeypatch):
    from types import SimpleNamespace
    import json
    import qiskit_ibm_runtime
    ns=runtime;submitted=[]
    class Job:
        def job_id(self):return f'local-{len(submitted)}'
        def result(self):
            return [SimpleNamespace(data=SimpleNamespace(meas=SimpleNamespace(
                get_counts=lambda:{'000000000':7,'111111111':3})))]
    class Sampler:
        def __init__(self,mode):pass
        def run(self,circuits,shots):submitted.append(circuits);return Job()
    monkeypatch.setattr(qiskit_ibm_runtime,'SamplerV2',Sampler)
    monkeypatch.setitem(ns,'transpile',lambda circuits,**kwargs:circuits)
    for pm,em in [(4,8),(8,32),(4,8)]:
        ns['select_circuits'](pm,em)
        meter=ns['Measurement']('ibm',shots=10,service=object(),backend=SimpleNamespace(name='local'),
                                 run_hardware=True,directory=tmp_path)
        np.testing.assert_allclose(meter.measure('both',1,[0.],'scan'),[.7])
    assert len(submitted)==2
    specs=[json.loads(p.read_text())['specification'] for p in tmp_path.glob('*.json')]
    assert {(s['prep_m'],s['evolution_m']) for s in specs}=={(4,8),(8,32)}
    assert len({s['circuit_digest'] for s in specs})==2
    assert submitted[0][0].count_ops()!=submitted[1][0].count_ops()
    ns['select_circuits'](6,16)


@pytest.mark.parametrize('name',['04_h4_circuit_assessment.ipynb','05_h4_iterative_bayes.ipynb'])
@pytest.mark.parametrize('pm,em',[(4,8),(8,32)])
def test_complete_notebook_with_nondefault_m(name,pm,em,tmp_path,monkeypatch):
    import matplotlib.pyplot as plt
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(plt,'show',lambda:None)
    ns={}
    for cell in nbformat.read(ROOT/name,4).cells:
        if cell.cell_type!='code':continue
        source=cell.source
        if source.startswith('PREP_M ='):
            source=source.replace('PREP_M = 6',f'PREP_M = {pm}').replace('EVOLUTION_M = 16',f'EVOLUTION_M = {em}')
        exec(source,ns)
    assert ns['SELECTED_PREP_M']==pm and ns['SELECTED_EVOLUTION_M']==em
    if name.startswith('04'):
        assert f'M{pm} prep + M{em} evolution L=1' in ns['COSTS']
        np.testing.assert_allclose(ns['L1_METRICS']['both']['peak'],ns['SELECTED_RESULTS']['fringe']['apparent_gap'],atol=1e-10)
        assert len(ns['M_ASSESSMENT'])==9
        for row in ns['M_ASSESSMENT']:
            prep=ns['CIRCUIT_CATALOG']['preparation'][row['prep_m']]
            evo=ns['CIRCUIT_CATALOG']['evolution'][row['evolution_m']]
            assert row['cx_L1']==prep['resources']['A']['cx']+prep['resources']['inverse']['cx']+evo['resources']['cx']
            assert row['state_S0']==prep['fidelities']['S0']
            assert row['operator_fidelity']==evo['operator_fidelity']
            assert row['curve'].shape==ns['M_SCAN'].shape
        chosen=next(r for r in ns['M_ASSESSMENT'] if (r['prep_m'],r['evolution_m'])==(pm,em))
        np.testing.assert_allclose(chosen['curve'],ns['probability']('both',1,ns['M_SCAN']),atol=0)

    else:
        assert ns['history']
        assert all(row['prep_m']==pm and row['evolution_m']==em for row in ns['history'])
        assert all(row['prep_m']==pm and row['evolution_m']==em
                   for rows in ns['comparison'].values() for row in rows)
        assert set(ns['M_BAYES_COMPARISON'])=={(p,e) for p in [4,6,8] for e in [8,16,32]}
        for pair,rows in ns['M_BAYES_COMPARISON'].items():
            assert rows
            assert all((r['prep_m'],r['evolution_m'])==pair and r['mode']=='exact' and r['case']=='both'
                       and r['shots'] is None for r in rows)
            assert rows[0]['prior_mu']==ns['INITIAL_MU']
            assert rows[0]['prior_variance']==ns['INITIAL_VARIANCE']
            for a,b in zip(rows,rows[1:]):
                assert a['posterior_mu']==b['prior_mu'] and a['posterior_variance']==b['prior_variance']
        assert ns['M_BAYES_REFERENCE']==ns['comparison']['exact']
        assert ns['meter'].selection['circuit_digest']==ns['CIRCUIT_DIGEST']
        import json
        assert json.loads((tmp_path/'outputs/exercise5_updates.json').read_text())==ns['history']

