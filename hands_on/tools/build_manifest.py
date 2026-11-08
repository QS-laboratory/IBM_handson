from pathlib import Path
import hashlib,json
from importlib.metadata import version
ROOT=Path(__file__).resolve().parents[1]
files={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((ROOT/'data').rglob('*')) if p.is_file()}
manifest={'format_version':1,'molecular_scope':'fixed planar teaching geometries; STO-3G; pi-CASCI; vertical S0-T1 gaps',
          'software':{p:version(p) for p in ['qiskit','qiskit-ibm-runtime','numpy','scipy','matplotlib','pyscf']},
          'hardware_results_included':False,'sha256':files}
(ROOT/'data/MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Hashed',len(files),'data files')
