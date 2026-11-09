"""Recompose existing self-study notebooks into five participant exercises.

Run from any directory. Generated notebooks embed their fixtures and helpers.
Rebuilding clears saved outputs; follow with execute_notebooks.py.
"""
from pathlib import Path
import base64
import copy
import hashlib
import json
import nbformat as nb
import argparse
import numpy as np

parser=argparse.ArgumentParser()
parser.add_argument("--only",choices=["01","02","03","04","05"])
args=parser.parse_args()

ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT.parents[1]/'archive/self_study'

def md(s): return nb.v4.new_markdown_cell(s.strip())
def code(s): return nb.v4.new_code_cell(s.strip())
def read(name): return nb.read(STUDY/name, 4)
def picked(source, ids): return [copy.deepcopy(source.cells[i]) for i in ids]

def add_answer(cell, text):
    cell.source += '\n\n<details><summary>解答例</summary>\n\n'+text.strip()+'\n\n</details>'

def write(name, cells, sources):
    if args.only and not name.startswith(args.only+"_"): return
    notebook = nb.v4.new_notebook(cells=cells, metadata={
        'kernelspec': {'display_name':'Python 3', 'language':'python','name':'python3'},
        'language_info': {'name':'python'},
        'teaching_sources': sources,
    })
    for cell in notebook.cells:
        if cell.cell_type == 'code':
            cell.outputs = []; cell.execution_count = None
    nb.validate(notebook)
    nb.write(notebook, ROOT/name)
    print(name, len(cells), 'cells')

INSTALL = code('''# 初回だけ必要に応じてコメントを外して実行し、カーネルを再起動してください。
# %pip install "numpy>=1.26,<3" "scipy>=1.12,<2" "matplotlib>=3.8,<4" "qiskit==2.4.1" "pylatexenc>=2.10,<3"
# 実機モードのみ追加: %pip install "qiskit-ibm-runtime>=0.43,<0.44"
''')

# 1: self-contained explanation, abstract circuits and B/C Bayesian examples.
from matrix_notebook import build_matrix_cells
write('01_matrix_phase_estimation.ipynb',
      build_matrix_cells(read('phase_estimation_comparison.ipynb'), INSTALL),
      ['archive/self_study/phase_estimation_comparison.ipynb'])

# 2: configuration/spin blocks, exact one-qubit encoding and molecular phase estimation.
from h2_notebook import build_h2_cells
with np.load(STUDY/'molecular_jw/data/h2.npz',allow_pickle=False) as fixture:
    ms0_block=fixture['hamiltonian'][np.ix_([5,10,6,9],[5,10,6,9])].real.tolist()
write('02_h2_one_qubit.ipynb',
      build_h2_cells(read('molecular/molecular_phase_estimation.ipynb'),INSTALL,ms0_block),
      ['archive/self_study/molecular/molecular_phase_estimation.ipynb','archive/self_study/molecular_jw/data/h2.npz'])

# 3: use the audited H2 core and molecular fixture, without author optimization or BPE/QPE sections.
fixture = base64.b64encode((STUDY/'molecular_jw/data/h2.npz').read_bytes()).decode()
init = (STUDY/'molecular_jw/tools/jw_core.py').read_text() + f'''
import base64, io
with np.load(io.BytesIO(base64.b64decode({fixture!r})), allow_pickle=False) as saved:
    DATA = {{k:saved[k].copy() for k in saved.files}}
from IPython.display import display
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.synthesis import LieTrotter, SuzukiTrotter

def draw(qc):
    display(qc.draw(output='mpl', fold=100)); plt.close('all')
'''
from h2_jw_notebook import build_h2_jw_cells
write('03_h2_jw_bpde.ipynb',
      build_h2_jw_cells(read('molecular_jw/h2_jw_phase_estimation.ipynb'),INSTALL,init,
                       (STUDY/'molecular_jw/tools/jw_core.py').read_text()),
      ['archive/self_study/molecular_jw/tools/jw_core.py','archive/self_study/molecular_jw/data/h2.npz',
       'archive/self_study/molecular_jw/h2_jw_phase_estimation.ipynb'])

# Lightweight bootstrap: pinned data and runtime files are cached on demand.
DATA=ROOT/'data/h4_pec_v1'
index=json.loads((DATA/'manifest.json').read_text())
catalog=json.loads((DATA/'R_2p00000000/catalog.json').read_text())
for name in ['preparation','evolution']:
    catalog[name]={int(k) if k!='exact' else k:v for k,v in catalog[name].items()}
h4init=code((ROOT/'tools/h4_bootstrap.py').read_text())
h4init.metadata['jupyter']={'source_hidden':True}
modeltext=md('''## 共通条件と測定方式

直鎖H₄、選択した隣接距離R（既定2 Å）、STO-3G、全4電子・4空間軌道、**Boys局在化＋軌道順最適化**。
系8＋補助1 qubit、モード順α₀…α₃,β₀…β₃、補助q₈です。
S₀と最低三重項T₀のMₛ=0成分（コードT）を使います。Δt=0.1、t=LΔt、β=εt。
準備・発展とも同じ軌道の保存済み回路を内蔵し、符号付きFock順変換も発展回路に含みます。
S₀・Tの両枝に共通の準備Mと時間発展Mを、冒頭の設定セルで独立に選びます（既定6・16）。Mは回路の層数で、時間反復LやMPSのbond dimension χとは別です。
最適化は再実行しません。初期化セルはデータ・補助関数を含むため折り畳んでいます。

| MODE | 測定値 | SHOTSの意味 |
|---|---|---|
| `exact` | 選んだ回路の厳密確率 | 使用しない |
| `shots` | 理想確率からの二項模擬測定 | 各探索点のshots |
| `ibm` | 実回路の全9ビット0の頻度 | 各探索点の実機shots |

`exact`はshot noiseがないという意味です。CASE=`both`で近似部品を指定した場合は、その回路近似が残ります。
`shots`にゲートノイズは含みません。全0頻度は全shotsを分母とします。
実機では選択したMの実回路を使用します。厳密回路との対照はシミュレーションです。
''')
config=code('''MODE = 'exact'     # 'exact', 'shots', 'ibm'
CASE = 'both'      # 'exact', 'prep_only', 'evolution_only', 'both'; ibmはbothのみ
SHOTS = 1024       # 各探索点のshots。MODE='exact'では使わない
SEED = 42
BACKEND_NAME = '' or globals().get('IBM_BACKEND_NAME', '')  # 冒頭の選択を引き継ぐ
INSTANCE = '' or globals().get('IBM_INSTANCE', '')  # 冒頭のCRNを引き継ぐ
RUN_HARDWARE = False  # 実機の費用を確認後、送信する場合にTrue
MAX_HARDWARE_L = 8    # 実機の反復上限。増やす前に回路費用を確認する
service = backend = None
if MODE == 'ibm':
    if SELECTED_EVOLUTION_M == 'exact':
        raise ValueError('厳密発展の密行列回路はシミュレーション参照です。実機にはEVOLUTION_M=8/16/32を選びます')
    if CASE != 'both':
        raise ValueError("MODE='ibm'ではCASE='both'を選んでください（選択した実回路）")
    if not BACKEND_NAME:
        raise ValueError('BACKEND_NAMEを設定してください')
    service, backend = connect_ibm(BACKEND_NAME, INSTANCE)
    hardware_preview(backend, L=1)  # コンパイルのみ。送信しない
meter = Measurement(MODE, SHOTS, SEED, service, backend, RUN_HARDWARE, MAX_HARDWARE_L)
MEASUREMENT_READY = MODE != 'ibm' or RUN_HARDWARE
print('mode:', MODE, 'case:', CASE, 'FCI gap / Ha:', GAP)
if not MEASUREMENT_READY:
    print('Preview only: inspect costs, set RUN_HARDWARE=True, then rerun this setting cell and following cells.')
''')
hardwaretext=md('''### 実機モードの使い方

冒頭の接続準備を実行してから、MODE=`ibm`、CASE=`both`へ変更します。
CRN・backend名・認証は引き継ぎます。接続準備をスキップした場合は先に戻って実行してください。
`hardware_preview(backend,L=...)`で接続制約を含む回路費用を確認してから、
RUN_HARDWAREをTrueにして設定セル以降を実行します。
ジョブIDは送信直後に`outputs/ibm_jobs/`へ保存し、同条件の再実行では再投入せず結果を回収します。
Colab切断に備えて、このフォルダをDrive等へ保存してください。フォルダを失うと自動復旧できません。
未確定の送信マーカーが残った場合はIBM側のジョブを確認してから復旧します。
実機は認証・利用枠・待ち時間が必要です。この教材作成時には送信していません。

[IBM Sampler入出力](https://eu-de.quantum.cloud.ibm.com/docs/en/guides/sampler-input-output)
''')
checks=code('''# 同じ回路の行列確率と、全9 qubit状態ベクトルの全0確率を照合。
for eps in [GAP, GAP+.25]:
    actual=Statevector.from_instruction(experiment_circuit(1,eps)).probabilities()[0]
    np.testing.assert_allclose(actual,probability('both',1,eps),atol=1e-10)
np.testing.assert_allclose(probability('exact',3,np.array([GAP,GAP+.1])),
    (1+np.cos(np.array([0,.1])*3*DT))/2,atol=1e-10)
assert GAP>0
np.testing.assert_allclose(GAP,float(MODEL['energy_T']-MODEL['energy_S0']),atol=1e-12)
print('Circuit / amplitude / exact cosine checks passed')
''')
sources=['exercises/data/h4_pec_v1/manifest.json','exercises/data/h4_pec_v1/exact_preparation.py',
         'exercises/tools/h4_measurement.py','exercises/tools/h4_snapshot.py','exercises/tools/h4_bootstrap.py','exercises/tools/ibm_setup.py']
from h4_assessment_notebook import build_h4_assessment_cells
fit_info={'chi':6,'mpo_microsteps':100,'mpo_cutoff':catalog['evolution'][16]['report']['reference_cutoff'],
          'refinement_cap':catalog['evolution'][16]['report']['refinement_cap'],
          'full_convergence_claimed':catalog['evolution'][16]['report']['full_convergence_claimed']}
write('04_h4_circuit_assessment.ipynb',
      build_h4_assessment_cells(INSTALL,h4init,modeltext,config,hardwaretext,checks,fit_info),
      sources+['exercises/tools/h4_inventory.py','exercises/tools/h4_pec_exercise.py',
               'exercises/tools/h4_pec_scan.py','exercises/data/h4_inventory_reference.json'])

from h4_bayes_notebook import build_h4_bayes_cells
write('05_h4_iterative_bayes.ipynb',
      build_h4_bayes_cells(INSTALL,h4init,modeltext,config,hardwaretext,checks),sources)
