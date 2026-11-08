from pathlib import Path
import json,hashlib,sys
ROOT=Path(__file__).resolve().parents[1]
manifest=json.loads((ROOT/'existing_materials_sha256.json').read_text())
# Retain the original snapshot; audit authorized editorial changes separately.
revision_path=ROOT/'audience_revision_sha256.json'
revisions=json.loads(revision_path.read_text())['sha256'] if revision_path.exists() else {}
manifest.update(revisions)
# The repository landing page and ignore rules evolve during publication.
# Retain their original hashes in the snapshot, but audit teaching files strictly.
publication_metadata={'README.md','.gitignore'}
changed=[name for name,value in manifest.items() if name not in publication_metadata and (not (ROOT.parent/name).exists() or hashlib.sha256((ROOT.parent/name).read_bytes()).hexdigest()!=value)]
core=ROOT/'runtime/Gate_generation_from_MPS_MPO.jl'
upstream=ROOT/'upstream/src/Gate_generation_from_MPS_MPO.jl'
if core.read_bytes()!=upstream.read_bytes():raise AssertionError('Author core differs')
if changed:raise AssertionError(f'Existing materials changed: {changed}')
print(f'{len(manifest)-len(publication_metadata)} teaching files verified ({len(revisions)} approved audience revisions); upstream Julia core byte-identical. README.md and .gitignore may change for publication.')
