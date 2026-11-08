from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from author_analysis import evaluate,run_qpde
if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--geometry',type=int,default=0)
    parser.add_argument('--names',nargs='+',default=['canonical','boys','boys_ordered'])
    args=parser.parse_args()
    metrics=[]
    for name in args.names:
        result=evaluate(name,args.geometry)
        print(json.dumps(result,indent=2),flush=True)
        metrics.append(result)
        run_qpde(name,args.geometry,exact_probabilities=True,reference_control=True)
        run_qpde(name,args.geometry,exact_probabilities=False,reference_control=True)
        for exact in (True,False):
            result=run_qpde(name,args.geometry,exact_probabilities=exact)
            final=result['iterations'][-1]
            print(name,result['source'],result['status'],'last mean',final.get('posterior_mean'),'exact',result['exact_gap'],flush=True)
    (ROOT/f'results/geometry_{args.geometry}/comparison.json').write_text(json.dumps(metrics,indent=2)+'\n')
