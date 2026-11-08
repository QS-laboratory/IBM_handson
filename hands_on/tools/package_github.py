"""Package local materials for the designated GitHub repository; never uploads."""
from pathlib import Path
import zipfile,json
ROOT=Path(__file__).resolve().parents[2]
output=ROOT/'output';output.mkdir(exist_ok=True)
archive=output/'IBM_handson-v2-preview.zip'
files=[ROOT/'README.md',ROOT/'.gitignore']
for folder in [ROOT/'hands_on', ROOT/'versions/v1_2026-11-04']:
    for p in folder.rglob('*'):
        relative=p.relative_to(folder)
        if p.is_file() and not any(part in {'__pycache__','.venv','.ipynb_checkpoints','outputs','output'} for part in relative.parts):files.append(p)
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as bundle:
    for path in sorted(files):bundle.write(path,path.relative_to(ROOT))
with zipfile.ZipFile(archive) as bundle:
    assert bundle.testzip() is None
    names=bundle.namelist()
    assert 'hands_on/data/circuits/ethylene_prep_compact.qpy' in names
    assert not any(name.endswith('.pdf') for name in names)
    assert all('__pycache__' not in name for name in names)
print(archive)
print('Files:',len(files),'Bytes:',archive.stat().st_size)
