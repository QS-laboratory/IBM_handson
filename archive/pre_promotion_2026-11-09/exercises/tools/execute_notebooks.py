"""Execute every notebook in a fresh Python/IPython process and temporary cwd.

This captures notebook outputs without requiring kernel sockets. It is not a
Colab UI or IBM hardware test. --worker is an internal execution option.
"""
import argparse
import atexit
import base64
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import nbformat

ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser()
parser.add_argument('--worker',type=Path)
args=parser.parse_args()
if args.worker is None:
    for path in sorted(ROOT.glob('0*.ipynb')):
        subprocess.run([sys.executable,__file__,'--worker',str(path)],check=True,
            env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MPLBACKEND':'Agg'})
else:
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output
    from traitlets.config import Config
    path=args.worker.resolve()
    os.environ.setdefault('H4_DATA_DIR',str(ROOT/'data/h4_pec_v1'))
    qa=ROOT/'outputs/qa'/path.stem
    qa.mkdir(parents=True,exist_ok=True)
    os.environ['MPLCONFIGDIR']=str(qa/'mplconfig')
    os.environ.setdefault('MPLBACKEND','Agg')
    notebook=nbformat.read(path,4)
    started=time.monotonic()
    with tempfile.TemporaryDirectory(prefix='ftqc-exercise-') as work:
        os.environ['IPYTHONDIR']=str(Path(work)/'ipython')
        shell=InteractiveShell.instance(config=Config({'HistoryManager':{'enabled':False}}))
        atexit.register(os.chdir, str(ROOT))
        import matplotlib.pyplot as plt
        from matplotlib_inline.backend_inline import configure_inline_support
        plt.switch_backend('module://matplotlib_inline.backend_inline')
        configure_inline_support(shell,'inline')
        os.chdir(work)
        count=0;images=0;stderr=[]
        for index,cell in enumerate(notebook.cells):
            if cell.cell_type!='code':continue
            count+=1
            with capture_output(stdout=True,stderr=True,display=True) as captured:
                result=shell.run_cell(cell.source,store_history=False)
            if result.error_in_exec or result.error_before_exec:
                print(captured.stdout,captured.stderr)
                raise RuntimeError(f'{path.name}: cell {index} failed') from (result.error_in_exec or result.error_before_exec)
            cell.execution_count=count;cell.outputs=[]
            if captured.stdout:cell.outputs.append(nbformat.v4.new_output('stream',name='stdout',text=captured.stdout))
            if captured.stderr:
                cell.outputs.append(nbformat.v4.new_output('stream',name='stderr',text=captured.stderr))
                stderr.append({'cell':index,'text':captured.stderr})
            for out in captured.outputs:
                cell.outputs.append(nbformat.v4.new_output('display_data',data=out.data,metadata=out.metadata))
                if 'image/png' in out.data:
                    (qa/f'figure_{images:02d}_cell_{index:02d}.png').write_bytes(base64.b64decode(out.data['image/png']))
                    images+=1
            print(path.name,f'cell {index}: passed',flush=True)
        # Preserve numerical/measurement artifacts produced in the temporary cwd.
        import shutil
        if Path('outputs').exists():shutil.copytree('outputs',qa/'results',dirs_exist_ok=True)
        os.chdir(ROOT)
    report=dict(notebook=path.name,status='passed',code_cells=count,images=images,stderr=stderr,
        elapsed_seconds=round(time.monotonic()-started,2),qpu_used=False,colab_tested=False,
        method='Fresh Python/IPython process, temporary cwd; no Jupyter kernel/UI test')
    notebook.metadata['verification']={k:v for k,v in report.items() if k!='stderr'}
    nbformat.validate(notebook);nbformat.write(notebook,path)
    (qa/'execution_report.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report),flush=True)
