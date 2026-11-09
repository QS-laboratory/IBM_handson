"""Mock remote files to check lazy download, cache reuse, offline use and integrity."""
from pathlib import Path
import io,urllib.request
import nbformat
import pytest
ROOT=Path(__file__).resolve().parents[1]
SOURCE=(ROOT/'tools/h4_bootstrap.py').read_text()


def test_lazy_download_reuse_and_offline_restart(tmp_path,monkeypatch):
    cache=tmp_path/'data';monkeypatch.setenv('H4_DATA_DIR',str(cache));calls=[]
    def remote(url,timeout):
        relative=url.split('/hands_on/exercises/',1)[1]
        calls.append(relative)
        return io.BytesIO((ROOT/relative).read_bytes())
    monkeypatch.setattr(urllib.request,'urlopen',remote)
    ns={};exec(SOURCE,ns)
    initial=len(calls)
    assert initial>0 and (cache/'manifest.json').is_file()
    assert not (cache/'R_2p00000000/evolution_exact.qpy').exists()
    assert not (cache/'R_3p00000000').exists()
    assert len(SOURCE)<10000 and "ASSETS={'" not in SOURCE
    ns['select_geometry']('R_3p00000000')
    assert len(calls)>initial and not (cache/'R_3p00000000/evolution_exact.qpy').exists()
    before=len(calls);ns['select_geometry']('R_2p00000000');assert len(calls)==before
    def no_network(*args,**kwargs):pytest.fail('Offline restart attempted network access')
    monkeypatch.setattr(urllib.request,'urlopen',no_network)
    offline={};exec(SOURCE.replace('H4_OFFLINE = False','H4_OFFLINE = True'),offline)
    offline['select_geometry']('R_3p00000000')
    with pytest.raises(FileNotFoundError,match='Offline: missing'):
        offline['select_circuits'](4,16)
    assert offline['SELECTED_PREP_M']==6
    (cache/'R_2p00000000/model.npz').write_bytes(b'corrupt')
    with pytest.raises(ValueError,match='checksum mismatch'):
        offline['_h4_data']('R_2p00000000/model.npz')


def test_notebooks_contain_bootstrap_not_binary_payloads():
    for path in ROOT.glob('0[45]*.ipynb'):
        book=nbformat.read(path,4)
        init=next(c for c in book.cells if c.cell_type=='code' and c.source.startswith('H4_DATA_COMMIT ='))
        assert len(init.source)<10000
        assert "ASSETS={'" not in '\n'.join(c.source for c in book.cells)
        assert init.metadata['jupyter']['source_hidden']
