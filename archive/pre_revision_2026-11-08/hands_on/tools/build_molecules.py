from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from molecules import build_reference
if __name__=='__main__':
    for molecule in ['ethylene','butadiene']:
        build_reference(molecule,ROOT/'data/reference')
