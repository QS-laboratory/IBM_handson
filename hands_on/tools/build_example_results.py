from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from workshop import make_plan,simulate_results,write_json,infer_gap,load_model
if __name__=='__main__':
    for name in ('ethylene','butadiene'):
        for variant in ('compact','reference') if name=='ethylene' else ('compact',):
            for noise in (0.,.15):
                plan=make_plan(name,variant=variant,team='example',shots=1024)
                result=simulate_results(plan,noise_rate=noise)
                result['provenance']='Generated locally from distributed QPY circuits; synthetic counts, NOT IBM hardware results.'
                suffix='ideal' if noise==0 else 'depolarized'
                filename=f'{name}_{variant}_{suffix}.json'
                write_json(ROOT/'data/example_results'/filename,result)
                history=infer_gap(result)
                final=history[-1]
                meta,_=load_model(name)
                print(filename,'reference',meta['gap_hartree'],'posterior',final['mean_hartree'],'CI',final['ci95_hartree'],flush=True)
