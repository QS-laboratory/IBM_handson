"""PEC peaks must come from phase scans with explicit branch information."""
from pathlib import Path
import json,runpy
import nbformat
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[1]

@pytest.fixture(scope='module')
def runtime():
    ns={}
    n=nbformat.read(ROOT/'04_h4_circuit_assessment.ipynb',4)
    exec(next(c.source for c in n.cells if c.cell_type=='code' and c.source.startswith('H4_DATA_COMMIT =')),ns)
    exec((ROOT/'tools/h4_pec_scan.py').read_text(),ns)
    return ns


def test_static_inventory_distinguishes_saved_chi_and_gate_scopes():
    module=runpy.run_path(str(ROOT/'tools/h4_inventory.py'))
    text=module['circuit_inventory']()
    assert 'χ=5、χ=4/M=6' in text and '該当する保存回路なし' in text
    assert 'MPO化' in text and '256×256密行列unitary' in text
    cat=json.loads((ROOT/'data/h4_pec_v1/R_2p00000000/catalog.json').read_text())
    for m in ['4','6','8']:
        for state in ['S0','T']:
            assert f"{cat['preparation'][m]['fidelities'][state]:.10f}" in text
    for m in ['8','16','32']:
        assert f"{cat['evolution'][m]['operator_fidelity']:.12f}" in text
    assert '未分解・未評価' in text and 'IBM backend' in text


def test_scan_formulas_match_bpe_and_bpde_circuits(runtime):
    ns=runtime;ns['select_circuits'](6,16)
    L=3;u=np.linalg.matrix_power(ns['U_APPROX'],L)
    a,b=[np.vdot(psi,u@psi) for psi in ns['APPROX_STATES']]
    controlled=np.block([[np.eye(256),np.zeros((256,256))],[np.zeros((256,256)),u]])
    for beta in [0.,.4,3.4]:
        qc=ns['QuantumCircuit'](9)
        qc.compose(ns['load_circuit_asset']('S0'),range(8),inplace=True)
        qc.h(8);qc.unitary(controlled,range(9));qc.p(beta,8);qc.h(8)
        p=ns['Statevector'].from_instruction(qc).probabilities()[:256].sum()
        np.testing.assert_allclose(p,(1+np.real(np.exp(1j*beta)*a))/2,atol=1e-10)
        actual=ns['Statevector'].from_instruction(ns['experiment_circuit'](L,beta/(L*ns['DT']))).probabilities()[0]
        np.testing.assert_allclose(actual,abs(a+np.exp(1j*beta)*b)**2/4,atol=1e-10)


@pytest.mark.parametrize('pm,em,L',[('exact','exact',64),(4,'exact',64),(8,32,1024)])
def test_twenty_point_scan_and_branch_certificates(runtime,pm,em,L,tmp_path,monkeypatch):
    ns=runtime;monkeypatch.setenv('H4_DATA_DIR',str(ROOT/'data/h4_pec_v1'))
    ns['select_geometry']('R_2p00000000');ns['select_circuits'](6,16)
    meter=ns['Measurement']('exact');digest=ns['CIRCUIT_DIGEST']
    path=tmp_path/'scan.json'
    rows=ns['run_pec_scan'](pm,em,L,path)
    assert len(rows)==20 and json.loads(path.read_text())==rows
    assert ns['SELECTED_GEOMETRY']=='R_2p00000000' and ns['CIRCUIT_DIGEST']==digest
    meter.measure('both',1,[0.],'after_scan')  # The original meter remains usable.
    for row in rows:
        assert row['prep_m']==pm and row['evolution_m']==em and row['L']==L
        assert row['phase_certificate']['certified'] and row['phase_certificate']['hD_over_minimum']<=.2
        assert row['max_fit_residual']<1e-12
        np.testing.assert_allclose(row['T_peak_total']-row['S0_peak_total'],row['gap_peak'],atol=1e-12)
        beta=np.asarray(row['beta']);principal,res=ns['fit_phase_scan'](beta,row['BPE_probability'])
        np.testing.assert_allclose(principal,row['principal_peak'][0],atol=1e-12)
        if pm==em=='exact':
            np.testing.assert_allclose([row['S0_peak_total'],row['T_peak_total'],row['gap_peak']],
                [row['S0_FCI_total'],row['T_FCI_total'],row['gap_FCI']],atol=1e-10)
    if pm==4:
        assert max(abs(r['S0_expect_total']-r['S0_FCI_total']) for r in rows)>.01
    if em==32:
        assert max(abs(r['S0_peak_total']-r['S0_FCI_total']) for r in rows)>1e-5


def test_invalid_L_rejected_without_changing_selection(runtime):
    for value in [0,-1,1.5,True,1025]:
        with pytest.raises(ValueError,match='integer in 1..1024'):
            runtime['run_pec_scan'](L=value)
