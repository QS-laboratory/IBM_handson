"""Check the promoted collection and the real 06 submission guard without IBM access."""
from pathlib import Path
import re
import nbformat
import pytest
ROOT=Path(__file__).resolve().parents[1]


def load(number):
    paths=list(ROOT.glob(f'{number:02d}_*.ipynb'))
    assert len(paths)==1
    return nbformat.read(paths[0],4)


def test_complete_six_notebook_collection_and_navigation():
    assert len(list(ROOT.glob('0[1-6]_*.ipynb')))==6
    assert not list(ROOT.glob('tmp_0*.ipynb'))
    readme=(ROOT/'README.md').read_text()
    for number in range(1,7):
        book=load(number);nbformat.validate(book)
        path=next(ROOT.glob(f'{number:02d}_*.ipynb'))
        assert path.name in readme
        assert f'/blob/main/hands_on/exercises/{path.name}' in readme
        assert not any('tmp_0' in c.source for c in book.cells)


def test_h2_jw_version_keeps_the_same_gaussian_update_as_one_qubit():
    def helper(book):
        return next(c.source for c in book.cells if c.cell_type=='code' and 'def gaussian_updates(' in c.source)
    assert helper(load(2))==helper(load(3))


def test_hardware_setup_skips_authentication_by_default(monkeypatch):
    import builtins,subprocess
    monkeypatch.setattr(builtins,'input',lambda *a:pytest.fail('Unexpected prompt'))
    monkeypatch.setattr(subprocess,'check_call',lambda *a,**k:pytest.fail('Unexpected installation'))
    setup=next(c.source for c in load(6).cells if c.cell_type=='code' and c.source.startswith('USE_IBM ='))
    ns={};exec(setup,ns)
    assert ns['IBM_SERVICE'] is None


def test_hardware_submit_guard_rejects_missing_preview():
    source=next(c.source for c in load(6).cells if c.cell_type=='code' and c.source.startswith('RUN_HARDWARE ='))
    def no_measurement(*args,**kwargs):pytest.fail('Unpreviewed hardware access')
    ns={'backend':None,'preview':None,'Measurement':no_measurement}
    exec(source,ns)  # Default OFF does not call Measurement.
    with pytest.raises(RuntimeError,match='プレビュー'):
        exec(source.replace('RUN_HARDWARE = False','RUN_HARDWARE = True',1),ns)
