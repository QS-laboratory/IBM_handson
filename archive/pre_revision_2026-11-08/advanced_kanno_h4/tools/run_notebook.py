"""Local offline notebook execution with figure capture, no Jupyter socket/kernel."""
import os
os.environ.setdefault('MPLBACKEND','Agg')
os.environ.setdefault('MPLCONFIGDIR','/tmp/kanno-h4-mpl')
from pathlib import Path
import io,base64,contextlib,json
import nbformat
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
ROOT=Path(__file__).resolve().parents[1]
os.chdir(ROOT)
import IPython.display
path=ROOT/'H4_author_method.ipynb'
notebook=nbformat.read(path,4)
qa=ROOT/'output/qa';qa.mkdir(parents=True,exist_ok=True)
outputs=[];image_count=0

def display(value):
    global image_count
    if isinstance(value,Figure):
        stream=io.BytesIO();value.savefig(stream,format='png',dpi=115,bbox_inches='tight')
        data=stream.getvalue();(qa/f'figure_{image_count:02d}.png').write_bytes(data);image_count+=1
        outputs.append(nbformat.v4.new_output('display_data',data={'image/png':base64.b64encode(data).decode(),'text/plain':'<matplotlib figure>'},metadata={}))
    else:print(value)

def show(*args,**kwargs):
    for number in plt.get_fignums():display(plt.figure(number))
    plt.close('all')
IPython.display.display=display;plt.show=show
namespace={'__name__':'__main__'};count=0
for index,cell in enumerate(notebook.cells):
    if cell.cell_type!='code':continue
    count+=1;outputs=[];stream=io.StringIO()
    try:
        with contextlib.redirect_stdout(stream),contextlib.redirect_stderr(stream):exec(compile(cell.source,f'cell-{index}','exec'),namespace)
    except Exception:
        print('FAILED cell',index,stream.getvalue());raise
    if stream.getvalue():
        outputs.insert(0,nbformat.v4.new_output('stream',name='stdout',text=stream.getvalue()))
        print(stream.getvalue(),flush=True)
    cell.outputs=outputs;cell.execution_count=count
nbformat.validate(notebook)
notebook.metadata['verification']={'local_offline_execution':True,'method':'sequential Python exec with figure capture','code_cells':count,'figures':image_count,'colab_tested':False,'ibm_executed':False}
nbformat.write(notebook,path)
(qa/'execution.json').write_text(json.dumps(dict(notebook.metadata['verification']),indent=2)+'\n')
print('Executed',count,'cells with',image_count,'figures')
