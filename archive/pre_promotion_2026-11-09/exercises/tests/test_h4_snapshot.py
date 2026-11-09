"""Geometry and exact/approx controls use the frozen public data, with no QPU."""
from pathlib import Path
import json, shutil
import nbformat
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/h4_pec_v1'
INDEX=json.loads((DATA/'manifest.json').read_text())

@pytest.fixture(scope='module')
def runtime():
    ns={}
    book=nbformat.read(ROOT/'04_h4_circuit_assessment.ipynb',4)
    exec(next(c.source for c in book.cells if c.cell_type=='code' and c.source.startswith('H4_DATA_COMMIT =')),ns)
    return ns

@pytest.mark.parametrize('geometry',[r['id'] for r in INDEX['geometries']])
def test_each_geometry_runtime_and_exact_reference(runtime,geometry,monkeypatch):
    ns=runtime;monkeypatch.setenv('H4_DATA_DIR',str(DATA))
    ns['select_geometry'](geometry)
    row=next(r for r in INDEX['geometries'] if r['id']==geometry)
    assert ns['SELECTED_GEOMETRY']==geometry and float(ns['MODEL']['R'])==row['R']
    for pm,em in [('exact','exact'),('exact',32),(8,'exact'),(6,16)]:
        ns['select_circuits'](pm,em)
        meter=ns['Measurement']('exact')
        eps=np.array([ns['GAP'],ns['GAP']+.3])
        actual=ns['Statevector'].from_instruction(ns['experiment_circuit'](3,eps[1])).probabilities()[0]
        np.testing.assert_allclose(actual,ns['probability']('both',3,eps[1]),atol=1e-9)
        if pm==em=='exact':
            np.testing.assert_allclose(ns['probability']('both',3,eps),(1+np.cos((eps-ns['GAP'])*.3))/2,atol=1e-10)
        meter.measure('both',3,eps,'test')
        assert meter.records[-1]['geometry']==geometry
    ns['select_geometry']('R_2p00000000');ns['select_circuits'](6,16)


def test_invalid_or_corrupt_geometry_does_not_replace_active_data(runtime,tmp_path,monkeypatch):
    ns=runtime;monkeypatch.setenv('H4_DATA_DIR',str(DATA));ns['select_geometry']('R_2p00000000')
    with pytest.raises(ValueError,match='Unknown GEOMETRY'):ns['select_geometry']('R_unknown')
    geometry=INDEX['geometries'][0]['id']
    shutil.copyfile(DATA/'manifest.json',tmp_path/'manifest.json')
    shutil.copytree(DATA/geometry,tmp_path/geometry)
    (tmp_path/geometry/'model.npz').write_bytes(b'corrupt')
    monkeypatch.setenv('H4_DATA_DIR',str(tmp_path))
    with pytest.raises(ValueError,match='checksum mismatch'):ns['select_geometry'](geometry)
    assert ns['SELECTED_GEOMETRY']=='R_2p00000000'


def test_geometry_change_rejects_stale_meter_and_exact_evolution_on_ibm(runtime,monkeypatch):
    ns=runtime;monkeypatch.setenv('H4_DATA_DIR',str(DATA))
    ns['select_geometry']('R_2p00000000');ns['select_circuits'](6,16)
    old=ns['Measurement']('exact')
    ns['select_geometry']('R_3p00000000')
    with pytest.raises(ValueError,match='recreate Measurement'):old.measure('both',1,[0.],'old')
    ns['select_circuits']('exact','exact')
    meter=ns['Measurement']('ibm',run_hardware=True)
    with pytest.raises(ValueError,match='simulation reference'):meter.measure('both',1,[0.],'exact')
    with pytest.raises(ValueError,match='simulation reference'):ns['hardware_preview'](None)
    ns['select_geometry']('R_2p00000000');ns['select_circuits'](6,16)


@pytest.mark.parametrize('book',['04_h4_circuit_assessment.ipynb','05_h4_iterative_bayes.ipynb'])
@pytest.mark.parametrize('pm,em',[('exact',32),(8,'exact'),('exact','exact')])
def test_full_notebook_other_geometry_and_exact_choices(book,pm,em,tmp_path,monkeypatch):
    import matplotlib.pyplot as plt
    monkeypatch.chdir(tmp_path);monkeypatch.setenv('H4_DATA_DIR',str(DATA))
    monkeypatch.setattr(plt,'show',lambda:None)
    ns={}
    for cell in nbformat.read(ROOT/book,4).cells:
        if cell.cell_type!='code':continue
        source=cell.source
        if source.startswith('PREP_M ='):
            source=source.replace('PREP_M = 6',f'PREP_M = {pm!r}').replace('EVOLUTION_M = 16',f'EVOLUTION_M = {em!r}')
            source=source.replace("GEOMETRY = 'R_2p00000000'","GEOMETRY = 'R_3p00000000'")
        if source.startswith('RUN_PEC_SCAN =') and (pm,em)!=('exact',32):
            source=source.replace('RUN_PEC_SCAN = True','RUN_PEC_SCAN = False')
        exec(source,ns)
    assert (ns['SELECTED_GEOMETRY'],ns['SELECTED_PREP_M'],ns['SELECTED_EVOLUTION_M'])==('R_3p00000000',pm,em)
    if book.startswith('05'):
        assert ns['history']
        assert all(r['geometry']=='R_3p00000000' and r['prep_m']==pm and r['evolution_m']==em for r in ns['history'])
    else:
        if em=='exact':assert ns['COSTS'][f'M{pm} prep + Mexact evolution L=1']['cx'] is None


def test_pauli_reconstruction_for_all_frozen_geometries():
    from qiskit.quantum_info import SparsePauliOp,Operator
    for row in INDEX['geometries']:
        with np.load(DATA/row['id']/'model.npz') as saved:h=saved['hamiltonian']
        pauli=SparsePauliOp.from_operator(Operator(h),atol=1e-14,rtol=0)
        np.testing.assert_allclose(pauli.to_matrix(),h,atol=1e-12,rtol=0)
