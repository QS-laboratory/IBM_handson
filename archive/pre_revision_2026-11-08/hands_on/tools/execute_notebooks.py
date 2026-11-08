"""Execute saved notebooks sequentially in-process without Jupyter socket access.
Hardware/connect flags stay False. Rich figures are captured into notebook outputs.
"""
import os
os.environ.setdefault('MPLBACKEND','Agg')
os.environ.setdefault('MPLCONFIGDIR','/tmp/ftqc-workshop-mpl')
from pathlib import Path
import sys,io,base64,contextlib,traceback,json
import nbformat
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import IPython.display
ROOT=Path(__file__).resolve().parents[1]
os.chdir(ROOT)
sys.path.insert(0,str(ROOT/'src'))
QA=ROOT/'output/qa'
QA.mkdir(parents=True,exist_ok=True)
report=[]
for relative in ['overview/00_overview.ipynb','00_preparation.ipynb','01_circuits_and_ibm.ipynb','02_inference.ipynb','03_mps_mpo_advanced.ipynb']:
    path=ROOT/relative;notebook=nbformat.read(path,4)
    namespace={'__name__':'__main__'};outputs=[];image_count=0;code_count=0
    tag=path.stem
    def display(value):
        global image_count
        if isinstance(value,Figure):
            stream=io.BytesIO();value.savefig(stream,format='png',dpi=110,bbox_inches='tight')
            raw=stream.getvalue()
            (QA/f'{tag}_{image_count:02d}.png').write_bytes(raw)
            image_count+=1
            outputs.append(nbformat.v4.new_output('display_data',data={'image/png':base64.b64encode(raw).decode(),'text/plain':'<matplotlib figure>'},metadata={}))
        else:print(value)
    def show(*args,**kwargs):
        for number in plt.get_fignums():display(plt.figure(number))
        plt.close('all')
    IPython.display.display=display;plt.show=show
    for index,cell in enumerate(notebook.cells):
        if cell.cell_type!='code':continue
        outputs=[];buffer=io.StringIO();code_count+=1
        try:
            with contextlib.redirect_stdout(buffer),contextlib.redirect_stderr(buffer):
                exec(compile(cell.source,f'{relative}:cell-{index}','exec'),namespace)
        except Exception:
            print(relative,'cell',index,'FAILED',buffer.getvalue(),flush=True)
            raise
        text=buffer.getvalue()
        if text:outputs.insert(0,nbformat.v4.new_output('stream',name='stdout',text=text))
        cell.outputs=outputs;cell.execution_count=code_count
    nbformat.validate(notebook)
    notebook.metadata['verification']={'method':'sequential Python exec in one shared process, fresh namespace per notebook; hardware branches disabled',
                                     'python':sys.version.split()[0],'jupyter_kernel_ui_tested':False}
    nbformat.write(notebook,path)
    report.append({'notebook':relative,'executed_code_cells':code_count,'saved_figures':image_count,'status':'passed_local_offline'})
    print(relative,'passed',code_count,'code cells',image_count,'figures',flush=True)
(QA/'execution_report.json').write_text(json.dumps(report,indent=2)+'\n')
