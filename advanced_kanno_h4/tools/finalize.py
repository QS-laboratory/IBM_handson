"""Validate, render and package only after the author computations have completed."""
from pathlib import Path
import sys,json,hashlib,subprocess,zipfile
from importlib.metadata import version
ROOT=Path(__file__).resolve().parents[1]

def run(*args):subprocess.run([sys.executable,*args],cwd=ROOT,check=True)

if __name__=='__main__':
    if '--package-only' not in sys.argv:
        run('tools/analyze_all.py')
        run('-m','pytest','tests','-q','-o',f'cache_dir={ROOT}/.pytest_cache')
        run('tools/build_notebook.py')
        run('tools/run_notebook.py')
        run('tools/check_preservation.py')
    summary=json.loads((ROOT/'results/geometry_0/comparison.json').read_text())
    readme=ROOT/'README.md'
    text=readme.read_text().split('\n<!-- GENERATED_RESULTS -->')[0]
    text+='\n<!-- GENERATED_RESULTS -->\n## 実行済みの比較（直鎖ID 0）\n\n'
    text+='| 軌道 | 準備の重なり f | 状態忠実度 f² | 時間発展 δ | ショット版の最終推定 / Ha |\n|---|---:|---:|---:|---:|\n'
    for row in summary:
        name=row['orbitals']
        result=json.loads((ROOT/f'results/geometry_0/{name}_qpde_shots.json').read_text())
        valid=[r for r in result['iterations'] if 'posterior_mean' in r]
        estimate=valid[-1]['posterior_mean'] if valid else float('nan')
        text+=f"| {name} | {row['preparation_overlap']:.6f} | {row['preparation_fidelity_squared']:.6f} | {row['evolution_delta_recomputed']:.6f} | {estimate:.6f} |\n"
    text+='\n厳密ギャップは約0.015905 Haです。推定値のずれを隠さず、回路近似とGaussian推定の限界を議論します。\n'
    readme.write_text(text)
    upstream=ROOT/'runtime/Gate_generation_from_MPS_MPO.jl'
    provenance={'upstream_repository':'https://github.com/sk888ks/TQPDE_open',
                'upstream_commit':'9cd88464c5c50dda783150661813609092d7f2b6',
                'upstream_julia_core_modified':False,'core_sha256':hashlib.sha256(upstream.read_bytes()).hexdigest(),
                'software':{p:version(p) for p in ['numpy','scipy','qiskit','qiskit-nature','pyscf']},
                'julia_version':'1.11.1','ITensors':'0.6.16','PyCall':'1.96.4','Infiltrator':'1.7.0',
                'python_compatibility_packages':{'deap':'1.4.3','qiskit-algorithms':'0.4.0'},
                'dmrg_sweeps':200,'preparation_optimization_sweeps':1000,'evolution_optimization_sweeps':1000,
                'depth_preparation':8,'depth_evolution':8,'dt':.1,'reference_trotter_slices':100,
                'all_paper_h4_exact_gaps_match_printed_precision':True,
                'hardware_executed':False,'fire_opal_executed':False,'colab_tested':False,
                'existing_materials_unchanged':True,
                'notebook_execution':json.loads((ROOT/'output/qa/execution.json').read_text())}
    (ROOT/'PROVENANCE.json').write_text(json.dumps(provenance,indent=2)+'\n')
    allowed=[]
    excluded={'.julia_depot','.python_deps','.upstream_git_metadata','.venv','__pycache__','.git','.pytest_cache','.ipynb_checkpoints','output'}
    for path in ROOT.rglob('*'):
        if not path.is_file():continue
        rel=path.relative_to(ROOT)
        if any(part in excluded for part in rel.parts):continue
        if path.suffix=='.dat':continue
        # The finalizer's stdout log is still open while packaging.
        if path.name=='final_validation.log':continue
        if rel.parts[0]=='upstream' and not (rel.parts[1] in {'src','examples'} or path.name in {'README.md','LICENCE'}):continue
        allowed.append(path)
    manifest={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(allowed) if p.name!='MANIFEST.sha256.json'}
    (ROOT/'MANIFEST.sha256.json').write_text(json.dumps(manifest,indent=2)+'\n')
    if ROOT/'MANIFEST.sha256.json' not in allowed:allowed.append(ROOT/'MANIFEST.sha256.json')
    output=ROOT/'output';output.mkdir(exist_ok=True)
    package=output/'advanced_kanno_h4.zip'
    with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(allowed):archive.write(path,Path(ROOT.name)/path.relative_to(ROOT))
    with zipfile.ZipFile(package) as archive:assert archive.testzip() is None
    print('Ready:',package,'files',len(allowed),'bytes',package.stat().st_size)
