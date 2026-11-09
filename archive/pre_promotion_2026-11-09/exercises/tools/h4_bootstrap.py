H4_DATA_COMMIT = '4ab6ef3d06f4256b5c818feb8c3fa2af4381dd43'
# 保存先を指定する場合の例（Google Driveは先にmountする）:
# H4_CACHE_DIR = '/content/drive/MyDrive/FTQC_HandsOn/h4_pec_v1'
H4_CACHE_DIR = None  # None: H4_DATA_DIR環境変数 → 取得済みリポジトリ → ユーザーキャッシュ
H4_OFFLINE = False  # True: ダウンロードせず、既存ファイルのみを使用

import os, json, hashlib, base64, tempfile
from pathlib import Path
from urllib.request import urlopen

SNAPSHOT_DIGEST = '6910f1f0aa7d0f2e43499683ba3a50f83db7aee09f3c64b07c91a6815ebdd4b3'
ASSET_DIGEST = SNAPSHOT_DIGEST
_H4_URL = f'https://raw.githubusercontent.com/QS-laboratory/IBM_handson/{H4_DATA_COMMIT}/hands_on/exercises/'
_H4_RUNTIME = {'h4_measurement.py':'d4dea4170f5e925599dcbbbc3b2c0ee7f2da1e25d7844ec024759e415c3a35f4',
               'h4_snapshot.py':'163691a0c005c844e0e449411c59eb5192218402df999c77056cbd7b223be997'}


def snapshot_root():
    explicit=os.environ.get('H4_DATA_DIR') or H4_CACHE_DIR
    root=Path(explicit).expanduser().resolve() if explicit else None
    if root is None:
        for parent in [Path.cwd(),*Path.cwd().parents]:
            for candidate in [parent/'data/h4_pec_v1',parent/'hands_on/exercises/data/h4_pec_v1']:
                if (candidate/'manifest.json').is_file():root=candidate;break
            if root is not None:break
    if root is None:root=Path.home()/'.cache/ftqc_handson/h4_pec_v1'
    if str(root).startswith('/content/drive/MyDrive/') and not Path('/content/drive/MyDrive').is_dir():
        raise FileNotFoundError('先に drive.mount("/content/drive") でGoogle Driveを接続してください')
    manifest=root/'manifest.json'
    if manifest.is_file() and hashlib.sha256(manifest.read_bytes()).hexdigest()!=SNAPSHOT_DIGEST:
        raise ValueError('Snapshot manifest differs from the pinned h4-pec-v1 manifest')
    return root


def _h4_fetch(path, url, digest):
    path=Path(path)
    if path.is_file():
        raw=path.read_bytes()
    else:
        if H4_OFFLINE:raise FileNotFoundError(f'Offline: missing {path}; pre-download the required data')
        try:
            with urlopen(url,timeout=60) as response:raw=response.read()
        except Exception as exc:
            raise RuntimeError(f'H4データを取得できません: {path.name}. 接続またはH4_DATA_DIR/H4_CACHE_DIRを確認してください。') from exc
        if hashlib.sha256(raw).hexdigest()!=digest:raise ValueError(f'Snapshot checksum mismatch: {path.name}')
        path.parent.mkdir(parents=True,exist_ok=True)
        tmp=None
        try:
            with tempfile.NamedTemporaryFile(dir=path.parent,delete=False) as stream:
                tmp=Path(stream.name);stream.write(raw)
            os.replace(tmp,path)
        finally:
            if tmp is not None:tmp.unlink(missing_ok=True)
    if hashlib.sha256(raw).hexdigest()!=digest:raise ValueError(f'Snapshot checksum mismatch: {path.name}')
    return raw


def _h4_data(relative):
    digest=SNAPSHOT_DIGEST if relative=='manifest.json' else SNAPSHOT_INDEX['files'][relative]
    return _h4_fetch(snapshot_root()/relative,_H4_URL+'data/h4_pec_v1/'+relative,digest)


SNAPSHOT_INDEX=json.loads(_h4_data('manifest.json'))


class _H4Assets:
    """Read only the requested asset; no notebook-embedded binary strings."""
    def __init__(self, row):self.row=row
    def __getitem__(self, key):
        raw=_h4_data(self.row['id']+'/'+self.row['assets'][key])
        return base64.b64encode(raw).decode()


def _h4_geometry(geometry):
    row=next((r for r in SNAPSHOT_INDEX['geometries'] if r['id']==geometry),None)
    if row is None:raise ValueError('Unknown GEOMETRY; choose an ID from list_geometries()')
    cat=json.loads(_h4_data(geometry+'/catalog.json'))
    for name in ['preparation','evolution']:
        cat[name]={int(k) if k!='exact' else k:v for k,v in cat[name].items()}
    return _H4Assets(row),cat,json.loads(_h4_data(geometry+'/provenance.json'))


ASSETS,CIRCUIT_CATALOG,PROVENANCE=_h4_geometry('R_2p00000000')
exec(compile(_h4_data('exact_preparation.py'),'h4_pec_v1/exact_preparation.py','exec'),globals())
_h4_cache_root=snapshot_root
for _name,_digest in _H4_RUNTIME.items():
    _local=snapshot_root().parent.parent/'tools'/_name
    _path=_local if _local.is_file() else snapshot_root().parent/(snapshot_root().name+'_runtime')/H4_DATA_COMMIT/_name
    exec(compile(_h4_fetch(_path,_H4_URL+'tools/'+_name,_digest),_name,'exec'),globals())
# The pinned runtime supplies measurement/inference helpers. Replace its eager
# snapshot reader with the verified, per-file cache above.
snapshot_geometry=_h4_geometry
snapshot_root=_h4_cache_root


def download_h4_data(geometries='all',prep_m=(4,6,8,'exact'),evolution_m=(8,16,32,'exact')):
    """Optional prefetch for offline use. This never changes active circuits."""
    ids=[r['id'] for r in SNAPSHOT_INDEX['geometries']] if geometries=='all' else list(geometries)
    count=0
    for geometry in ids:
        row=next((r for r in SNAPSHOT_INDEX['geometries'] if r['id']==geometry),None)
        if row is None:raise ValueError('Unknown GEOMETRY')
        assets,cat,_=_h4_geometry(geometry)
        keys={'model'}
        for m in prep_m:keys.update(cat['preparation'][m]['keys'].values())
        for m in evolution_m:keys.add(cat['evolution'][m]['key'])
        for key in sorted(keys):assets[key];count+=1
    print(f'Cached {len(ids)} geometries / {count} assets in {snapshot_root()}')

print('H4 data cache:',snapshot_root(),'— missing files are downloaded on demand; existing files are verified and reused.')
