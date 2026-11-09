"""Mock remote files to check lazy download, cache reuse, offline use and integrity."""
from pathlib import Path
import io,json,urllib.request
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
    for path in ROOT.glob('0[56]*.ipynb'):
        book=nbformat.read(path,4)
        init=next(c for c in book.cells if c.cell_type=='code' and c.source.startswith('H4_DATA_COMMIT ='))
        assert len(init.source)<10000
        assert "ASSETS={'" not in '\n'.join(c.source for c in book.cells)
        assert init.metadata['jupyter']['source_hidden']


def pec_namespace():
    ns={'H4_INDEX':json.loads((ROOT/'data/h4_pec_v1/manifest.json').read_text())}
    exec((ROOT/'tools/h4_pec.py').read_text(),ns)
    return ns


def mock_pec_remote(monkeypatch):
    calls=[]
    def remote(url,timeout):
        relative=url.split('/data/h4_pec_v1/',1)[1]
        calls.append(relative)
        return io.BytesIO((ROOT/'data/h4_pec_v1'/relative).read_bytes())
    monkeypatch.setattr(urllib.request,'urlopen',remote)
    return calls


def test_pec_explicit_cache_download_and_restart(tmp_path,monkeypatch):
    monkeypatch.delenv('H4_DATA_DIR',raising=False)
    calls=mock_pec_remote(monkeypatch)
    cache=tmp_path/'local cache'
    ns=pec_namespace();ns['H4_CACHE_DIR']=str(cache)
    ns['saved_states']('R_2p00000000',6)
    ns['evolution_spectrum']('R_2p00000000',16)
    assert len(calls)==4  # model, two prepared states, one evolution circuit
    assert (cache/'R_2p00000000/evolution_16.qpy').is_file()
    def offline(*args,**kwargs):pytest.fail('Cache restart attempted a download')
    monkeypatch.setattr(urllib.request,'urlopen',offline)
    restart=pec_namespace();restart['H4_CACHE_DIR']=str(cache)
    restart['saved_states']('R_2p00000000',6)
    restart['evolution_spectrum']('R_2p00000000',16)
    path=cache/'R_2p00000000/model.npz';path.write_bytes(b'corrupt')
    with pytest.raises(ValueError,match='checksum mismatch'):
        restart['data_bytes']('R_2p00000000/model.npz')


def test_pec_destination_priority_and_defaults(tmp_path,monkeypatch):
    calls=mock_pec_remote(monkeypatch)
    repo=tmp_path/'repo';repo.mkdir();monkeypatch.chdir(repo)
    local=repo/'data/h4_pec_v1';local.mkdir(parents=True)
    (local/'manifest.json').write_text('{}')
    env=tmp_path/'environment';cache=tmp_path/'chosen'
    monkeypatch.setenv('H4_DATA_DIR',str(env))
    ns=pec_namespace();ns['H4_CACHE_DIR']=str(cache)
    ns['data_bytes']('R_2p00000000/model.npz')
    assert (env/'R_2p00000000/model.npz').is_file() and not cache.exists()
    monkeypatch.delenv('H4_DATA_DIR')
    ns['data_bytes']('R_2p00000000/model.npz')
    assert (cache/'R_2p00000000/model.npz').is_file() and len(calls)==2
    ns['H4_CACHE_DIR']=None
    assert ns['data_root']()==local
    outside=tmp_path/'outside';outside.mkdir();monkeypatch.chdir(outside)
    monkeypatch.setattr(Path,'home',classmethod(lambda cls:tmp_path/'home'))
    assert ns['data_root']()==tmp_path/'home/.cache/ftqc_handson/h4_pec_v1'


def test_pec_drive_mount_save_and_restart(tmp_path,monkeypatch):
    """Map Colab's mounted filesystem into a temporary directory, without Google auth."""
    monkeypatch.delenv('H4_DATA_DIR',raising=False)
    calls=mock_pec_remote(monkeypatch)
    drive='/content/drive/MyDrive';mount=tmp_path/'mounted MyDrive'
    class DrivePath(type(Path())):
        def __fspath__(self):
            value=str(self)
            return str(mount)+value[len(drive):] if value==drive or value.startswith(drive+'/') else value
        def resolve(self,*args,**kwargs):
            if str(self)==drive or str(self).startswith(drive+'/'):return self
            return super().resolve(*args,**kwargs)
    ns=pec_namespace();ns['Path']=DrivePath
    ns['H4_CACHE_DIR']=drive+'/FTQC_HandsOn/h4_pec_v1'
    with pytest.raises(FileNotFoundError,match='drive.mount'):
        ns['data_bytes']('R_2p00000000/model.npz')
    assert calls==[]
    mount.mkdir()
    ns['saved_states']('R_2p00000000',6)
    assert (mount/'FTQC_HandsOn/h4_pec_v1/R_2p00000000/prep_6_S0.qpy').is_file()
    def offline(*args,**kwargs):pytest.fail('Drive restart attempted a download')
    monkeypatch.setattr(urllib.request,'urlopen',offline)
    restart=pec_namespace();restart['Path']=DrivePath
    restart['H4_CACHE_DIR']=ns['H4_CACHE_DIR']
    restart['saved_states']('R_2p00000000',6)


def test_pec_invalid_download_is_not_saved(tmp_path,monkeypatch):
    monkeypatch.delenv('H4_DATA_DIR',raising=False)
    monkeypatch.setattr(urllib.request,'urlopen',lambda *a,**k:io.BytesIO(b'invalid'))
    ns=pec_namespace();ns['H4_CACHE_DIR']=str(tmp_path/'cache')
    with pytest.raises(ValueError,match='Downloaded snapshot checksum mismatch'):
        ns['data_bytes']('R_2p00000000/model.npz')
    assert not (tmp_path/'cache/R_2p00000000/model.npz').exists()


def test_pec_notebook_embeds_current_loader():
    book=nbformat.read(ROOT/'04_h4_circuit_assessment.ipynb',4)
    init=next(c for c in book.cells if c.cell_type=='code' and 'H4_INDEX = ' in c.source)
    assert 'H4_CACHE_DIR = None' in init.source
    assert init.source.endswith((ROOT/'tools/h4_pec.py').read_text())
