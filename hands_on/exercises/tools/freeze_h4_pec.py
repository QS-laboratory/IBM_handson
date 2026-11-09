"""Explicit instructor-only export. Never run implicitly during notebook generation."""
from pathlib import Path
import base64, hashlib, json, shutil, sys
import numpy as np
from h4_catalog import build_catalog

ROOT=Path(__file__).resolve().parents[1]
STUDY=ROOT.parents[1]/'archive/self_study'
DEST=ROOT/'data/h4_pec_v1'

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    if (DEST/'manifest.json').exists():
        raise SystemExit('Snapshot already frozen; export a new version instead of overwriting it.')
    prep=STUDY/'h4_mps_preparation_study/results/pec_20'
    old=STUDY/'h4_mps_preparation_study/results/both_states_M6_evolution_M16_L1'
    index=dict(version='h4-pec-v1',dt=.1,modes='alpha0..3,beta0..3',
        geometry='linear equally spaced H4; R is adjacent distance in angstrom',
        orbitals='Boys + ordering; same saved coefficients per geometry',geometries=[],files={})
    for directory in sorted(prep.glob('R_*')):
        geometry=directory.name
        target=DEST/geometry
        if (target/'catalog.json').exists():
            catalog=json.loads((target/'catalog.json').read_text())
            prior=json.loads((target/'provenance.json').read_text())
            assert all(sha(STUDY.parent/name)==value for name,value in prior['catalog_sources'].items())
            filenames={p.stem:p.name for p in target.iterdir() if p.suffix in ['.qpy','.npz']}
        else:
            paths={'model':directory/'model.npz','A':old/'A_basic.qpy',
                   'evolution':old/'evolution_M16_blocked_basic.qpy'}
            assets={k:base64.b64encode(p.read_bytes()).decode() for k,p in paths.items()}
            provenance={'results':json.loads((old/'results.json').read_text())}
            catalog=build_catalog(STUDY,assets,provenance,geometry,include_exact=True)
            # Source paths are relative to hands_on; retain hashes, never local absolute paths.
            provenance.pop('results')
            target=DEST/geometry;target.mkdir(parents=True,exist_ok=True)
            used={'model'}
            for info in catalog['preparation'].values():used.update(info['keys'].values())
            for info in catalog['evolution'].values():used.add(info['key'])
            filenames={}
            for key in sorted(used):
                filename=key+('.npz' if key=='model' else '.qpy')
                (target/filename).write_bytes(base64.b64decode(assets[key]));filenames[key]=filename
            (target/'catalog.json').write_text(json.dumps(catalog,indent=2)+'\n')
            (target/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
        with np.load(directory/'model.npz') as data:
            row=dict(id=geometry,R=float(data['R']),energy_S0=float(data['energy_S0']),
                     energy_T=float(data['energy_T']),enuc=float(data['enuc']),assets=filenames)
        index['geometries'].append(row)
        print('Frozen',geometry,flush=True)
    assert len(index['geometries'])==20
    core=(STUDY/'molecular_h4_jw/tools/h4_core.py').read_text()
    (DEST/'exact_preparation.py').write_text(core[:core.index('class Propagator:')])
    # Remove machine-specific path prefixes before pinning the public snapshot.
    for path in DEST.glob('R_*/catalog.json'):
        content=path.read_text().replace(str(STUDY.parent.parent)+'/', '')
        path.write_text(content)
    shutil.copyfile(STUDY/'molecular_h4_jw/UPSTREAM_LICENSE',DEST/'UPSTREAM_LICENSE')
    index['exporter_sha256']=sha(Path(__file__))
    index['catalog_builder_sha256']=sha(ROOT/'tools/h4_catalog.py')
    index['exact_helpers_source_sha256']=sha(STUDY/'molecular_h4_jw/tools/h4_core.py')
    for p in sorted(DEST.rglob('*')):
        if p.is_file():index['files'][str(p.relative_to(DEST))]=sha(p)
    (DEST/'manifest.json').write_text(json.dumps(index,indent=2)+'\n')
    print('Dataset MiB:',sum(p.stat().st_size for p in DEST.rglob('*') if p.is_file())/2**20)

if __name__=='__main__':main()
