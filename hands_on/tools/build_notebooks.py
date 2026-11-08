from pathlib import Path
import nbformat as nb
ROOT=Path(__file__).resolve().parents[1]

BOOTSTRAP='''# GitHub → Colab / ローカルの共通セットアップ
import sys, subprocess, importlib.util
from pathlib import Path
GITHUB_REPOSITORY = "https://github.com/QS-laboratory/IBM_handson.git"
GITHUB_REF = "main"
try:
    IS_COLAB = importlib.util.find_spec("google.colab") is not None
except ModuleNotFoundError:
    IS_COLAB = False
if IS_COLAB:
    subprocess.run([sys.executable, "-m", "pip", "install", "-q",
                    "qiskit==2.4.1", "qiskit-ibm-runtime==0.43.1", "pylatexenc==2.10"], check=True)
    checkout = Path("/content/IBM_handson")
    if not checkout.exists():
        subprocess.run(["git", "clone", GITHUB_REPOSITORY, str(checkout)], check=True)
    subprocess.run(["git", "-C", str(checkout), "fetch", "origin", GITHUB_REF, "--depth", "1"], check=True)
    subprocess.run(["git", "-C", str(checkout), "checkout", "--detach", "FETCH_HEAD"], check=True)
    WORKSHOP_ROOT = checkout / "hands_on"
else:
    candidates = [Path.cwd(), Path.cwd()/"hands_on", Path.cwd().parent]
    WORKSHOP_ROOT = next(p for p in candidates if (p/"src/workshop.py").exists())
sys.path.insert(0, str(WORKSHOP_ROOT / "src"))
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import display
from workshop import *
import qiskit
print("Qiskit:", qiskit.__version__)
print("教材:", WORKSHOP_ROOT)
print("実機への送信は、明示的に有効化したセルだけが行います。")
OUTPUT_DIR = Path("outputs")
OUTPUT_DIR.mkdir(exist_ok=True)

def show_circuit(qc):
    try:
        figure = qc.draw("mpl", fold=90)
        display(figure)
        plt.close(figure)
    except ImportError:
        print(qc.draw("text", fold=90))
'''

AUTH='''# 接続を試すときだけTrueに変更。APIキーをコードや出力へ書かないでください。
CONNECT_IBM = False
INSTANCE_CRN = ""  # 今回利用できる所属企業の枠、または講師が案内するFlex枠
BACKEND_NAME = ""  # 接続後、利用可能な実機から講師と相談して指定
service = None
if CONNECT_IBM:
    from qiskit_ibm_runtime import QiskitRuntimeService
    if not INSTANCE_CRN:
        raise ValueError("利用権限のあるインスタンスのCRNを指定してください")
    if IS_COLAB:
        from google.colab import userdata
        api_key = userdata.get("IBM_QUANTUM_API_KEY")
    else:
        from getpass import getpass
        api_key = getpass("IBM Quantum API key: ")
    service = QiskitRuntimeService(channel="ibm_quantum_platform", token=api_key, instance=INSTANCE_CRN)
    del api_key
    print("利用可能な実機:", [b.name for b in service.backends(operational=True, simulator=False)])
else:
    print("IBM接続はスキップ。ローカル演習と配布データの解析は続けられます。")
'''


def write(name,parts,title):
    cells=[]
    for kind,source in parts:
        cells.append((nb.v4.new_markdown_cell if kind=='md' else nb.v4.new_code_cell)(source.strip()))
    notebook=nb.v4.new_notebook(cells=cells,metadata={'audience':'participant',
      'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'},
      'language_info':{'name':'python','version':'3.12'},
      'workshop':{'title':title,'date':'2026-11-09','github_repository':'QS-laboratory/IBM_handson',
                  'audience':'participant','hardware_tested':False,'colab_tested':False}})
    nb.write(notebook,ROOT/name)

M='md';C='code'
write('overview/00_overview.ipynb',[
(M,'''# 分子の電子状態を、量子コンピュータで調べる
## IBM Quantum FTQC Hands-on workshop — 概観

対象：材料科学・量子化学・量子計算に関心のある企業人・学生。専門知識を前提にしません。

**今日の問い**：分子の2つの電子状態のエネルギー差を、量子コンピュータの測定からどのように推定するか。

**到達目標**
1. 状態準備・時間発展・測定の役割を説明する。
2. 回路の操作数と精度のトレードオフを比較する。
3. 測定からエネルギー差を推定し、不確かさと限界を述べる。

資料の図はコード実行で表示できます。IBMへの接続は行いません。'''),
(C,BOOTSTRAP),
(M,'''## 1. 材料の性質から、電子の状態へ

材料や分子が光を吸収する、電子を受け渡す、化学反応を起こすとき、電子の状態とエネルギーが関係します。
実材料では温度、環境、構造の変化、分子間相互作用も重要です。今日は対象を絞り、**原子核の位置を固定した1分子の電子状態**を扱います。

- 原子核の配置を決める。
- その配置における電子のエネルギーと状態を計算する。
- 状態間の差を読み、物性理解の一つの手掛かりにする。

**考えてみましょう**：業務や研究で関心があるのは、安定性、光、反応、電荷移動のどれでしょうか。'''),
(M,'''## 2. エチレンからポリエンへ：π電子と共役

炭素の二重結合にはσ結合とπ結合があります。平面に垂直なp軌道が重なり、π電子を分子全体へ広げます。
二重結合と単結合が交互に並ぶ共役鎖を、エチレン → 1,3-ブタジエン → 1,3,5-ヘキサトリエンと延ばして考えます。
エチレンはこの説明の入口となるアルケン、ブタジエンとヘキサトリエンはポリエンです。

今日はブタジエンを主題材とし、実機実行の初期設定は小さなエチレンです。
両分子の幾何は教材用に固定した理想化構造であり、最適化済み構造ではありません。'''),
(C,'''fig, axes = plt.subplots(1, 2, figsize=(10, 3))
for name, ax in zip(["ethylene", "butadiene"], axes): plot_molecule(name, ax)
plt.tight_layout()
plt.show()'''),
(M,'''### 分子軌道を見る

下図は、STO-3G基底で計算したπ空間の分子軌道を、分子平面から0.6 Å上の断面で示したものです。
色は波動関数の符号で、電荷の正負ではありません。符号全体を反転しても同じ軌道です。
軌道は1電子の表現であり、多電子状態全体と区別します。'''),
(C,'''meta, arrays = load_model("butadiene")
fig, axes = plt.subplots(1, 4, figsize=(12, 3))
for i, ax in enumerate(axes):
    values = arrays["orbital_values"][:, :, i]
    bound = np.max(abs(values))
    ax.contourf(arrays["orbital_x"], arrays["orbital_y"], values,
                levels=np.linspace(-bound, bound, 17), cmap="RdBu_r")
    ax.set_title(f"Active pi orbital {i+1}")
    ax.set_aspect("equal"); ax.axis("off")
plt.tight_layout(); plt.show()'''),
(M,'''## 3. 基底状態・励起状態・スピン

電子はフェルミ粒子です。複数の電子を配置すると、軌道の占有だけでなくスピンの組み合わせも状態を区別します。
一重項は全スピンS=0、三重項はS=1です。三重項には異なるスピン射影の3成分があり、この教材では回路用にM_S=0を選びます。

今回の対象は、同じ原子核配置での**最低一重項 S0 と最低三重項 T1 の垂直エネルギー差**です。
三重項は光物性や励起状態の過程を考える上で重要ですが、S0–T1差は通常の許容光吸収のエネルギーとは同一ではありません。
'''),
(C,'''fig, axes = plt.subplots(1, 2, figsize=(8, 3))
for name, ax in zip(["ethylene", "butadiene"], axes):
    meta, _ = load_model(name)
    gap = meta["gap_ev"]
    ax.hlines([0, gap], -.3, .3, colors=["#2563eb", "#c2410c"], linewidth=3)
    ax.annotate("", xy=(.45, gap), xytext=(.45, 0), arrowprops={"arrowstyle":"<->"})
    ax.text(.5, gap/2, f"{gap:.3f} eV", va="center")
    ax.text(-.55, 0, "S0", va="center"); ax.text(-.55, gap, "T1", va="center")
    ax.set(xlim=(-.8, 1.3), ylim=(-.5, 5.5), title=name, ylabel="Energy relative to S0 (eV)", xticks=[])
plt.tight_layout(); plt.show()
print("These are CASCI/STO-3G teaching-model values, not experimental excitation energies.")'''),
(M,'''## 4. 基底関数・電子相関・活性空間

**基底関数**：電子の空間的な広がりを表現する部品です。今日は小さなSTO-3Gを使います。

**電子相関**：電子どうしの影響を、独立した1電子の描像より詳しく扱うことです。
平均場のHartree–Fockだけでは記述しにくい状態があり、複数の電子配置を重ね合わせます。

**活性空間**：すべての軌道を同じ精度で扱う代わりに、今回はπ電子・π軌道を選んで、その中の配置を漏れなく扱うCASCIを使います。
残りの占有軌道は二重占有で凍結します。CASCIでは軌道を固定して配置を混ぜます。CASSCFの軌道最適化は行っていません。

| 分子 | 活性電子・空間軌道 | 系量子ビット | 補助込み |
|---|---|---:|---:|
| エチレン | 2電子・2軌道 | 4 | 5 |
| ブタジエン | 4電子・4軌道 | 8 | 9 |
| ヘキサトリエン（Kanno論文） | 6電子・6軌道 | 12 | 13 |

1空間軌道にα・βの2スピン軌道があり、通常のJordan–Wigner変換ではそれぞれを1量子ビットにします。
量子ビット削減は今回は行いません。基底関数を小さくする近似と、活性空間を絞る近似は別です。

**正解の階層**：教材模型の厳密解 ≠ 基底・活性空間を拡大した解 ≠ 実験値。
今回の基準値は、この固定構造・基底・活性空間の中の厳密解です。'''),
(M,r'''## 5. エネルギーを位相から読む：FTQCとの関係

固有状態の時間発展は $e^{-iHt}|E\rangle=e^{-iEt}|E\rangle$（$\hbar=1$）です。
2状態を重ね合わせると、エネルギー差が相対的な位相の変化として現れます。
適切に干渉させて測定すると、その差を推定できます。測定の1回でエネルギーが直接表示されるわけではありません。

1. **状態準備**：対象の2状態の重ね合わせを作る。
2. **時間発展**：相対的な位相を蓄積する。
3. **逆の状態準備と測定**：干渉を確率へ変換する。
4. **古典計算**：複数の測定結果からエネルギー差を推定する。

FTQCは誤り訂正を使って長い量子計算を信頼できる形で実行する構想です。
位相推定はその重要なアルゴリズムですが、今回使うIBM実機は誤り訂正済みFTQCではありません。
短く工夫した回路で、現在の実機でも原理を試します。'''),
(M,'''## 6. Kanno論文と今日の実験

Kannoらは、基底・励起状態の重ね合わせ準備と短時間発展をテンソルネットワークを使って浅い回路に近似し、位相差推定を実機で実行しました。
Hubbard模型は系32＋補助1量子ビット、分子の本編はヘキサトリエンの系12＋補助1量子ビットです。
分子では軌道の局在化と並べ替えも工夫し、最低一重項–三重項の差を対象にしています。

**論文を読むときの観点**：大きな回路を実行できたこと、信号が得られたこと、高精度で正解したこと、量子優位性を示したことは、それぞれ別の評価です。
論文の分子実機結果には大きな不確かさがあり、今回も誤差を含めて読み解きます。

**今日の教材との違い**
- 同じ分子系列とS0–T1差を扱いますが、論文の数値再現ではありません。
- 状態準備は実数RY/CX回路の事前最適化、時間発展はPauli積の短い回路を使う説明用実装です。
- 本編の回路最適化は密な小規模計算で、論文のMPS/MPOによる大規模最適化とは異なります。
- Fire Opalは基本演習の必須条件にしていません。

**当日の流れ**：準備 → 回路を比較し実機投入 → 結果から推定を更新 → 不確かさと今後の課題を討論。

**出典**
- Kanno et al., PNAS 122, e2425026122 (2025), 本文Fig. 3–5・補足S1–S3： https://doi.org/10.1073/pnas.2425026122
- 活性空間とCASCI： https://pyscf.org/user/mcscf.html
- 電子配置とFCI： https://pyscf.org/user/ci.html

**振り返り**：何の分子の、どの2状態の差を、どんな近似で求めるのか、自分の言葉で説明してください。''')
], '概観：材料・分子・電子状態と量子計算')

write('00_preparation.ipynb',[
(M,'''# 準備：Colabから量子回路を動かす

**できるようになること**：教材を開く／小さな回路を実行する／量子実機へ接続する／対象分子を確認する。

全員が基本操作を行い、実機への試行は各チームの代表者が担当します。
所属企業の利用枠を今回の用途で利用できる方はそちらを使い、それ以外の希望者は講師が案内するFlex枠を使います。

1. GitHubの「Open in Colab」をクリック。
2. Colabの「ファイル → ドライブにコピーを保存」で自分のコピーを作成。
3. 標準のCPUランタイムで、上から実行。GPUは不要です。
4. Google Driveのマウントも不要です。結果ファイルを残す際はダウンロードします。
'''),
(C,BOOTSTRAP),
(M,'''## 小さな回路：0/1の確率を変える
`theta`を0、π/2、πへ変えてみてください。状態ベクトルの各成分の絶対値の2乗が確率です。'''),
(C,'''from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector
theta = np.pi/2
small = QuantumCircuit(1)
small.ry(theta, 0)
show_circuit(small)
probabilities = Statevector.from_instruction(small).probabilities()
print("P(0), P(1) =", probabilities)
fig, ax = plt.subplots(figsize=(4, 2.5))
ax.bar(["0", "1"], probabilities); ax.set_ylim(0,1)
ax.set_ylabel("Probability"); plt.show()'''),
(M,'''## 2量子ビットに相関を作る
Hで重ね合わせを作り、CNOTで相関を作ります。理想計算と、1024回測定したときのばらつきを比較します。
Qiskitのビット列は右端がq0です。これは本編・分子データでも共通です。
第一版の発展教材だけは数式の上からの並びに変換する関数を使うため、変換規則を明記しています。'''),
(C,'''bell = QuantumCircuit(2)
bell.h(0); bell.cx(0, 1)
show_circuit(bell)
state = Statevector.from_instruction(bell)
state.seed(19)
print("Ideal probabilities:", state.probabilities_dict())
print("Simulated counts:", state.sample_counts(1024))'''),
(M,'''## IBM実機への接続

- IBM Quantum Platformで対象インスタンスへの利用権限を確認し、APIキーを用意します。
- Colab左側の鍵アイコン「シークレット」に `IBM_QUANTUM_API_KEY` を登録し、このノートへのアクセスを許可します。
- 利用するインスタンスのCRNを下のセルへ入力します。APIキーはセルへ直接書きません。
- `CONNECT_IBM = True` にして接続します。このセルは接続と一覧取得のみです。

手順： https://quantum.cloud.ibm.com/docs/en/guides/initialize-account

利用可能な実機が表示されない場合は、アカウント・インスタンスの権限・地域を講師と確認してください。'''),
(C,AUTH),
(M,'''### チーム代表者：Bell回路を実機へ送信

以下は実機利用枠を消費します。講師の合図で `RUN_BELL = True` にします。
ジョブIDは送信直後に保存します。Colabから切断されても、再投入せずジョブIDから結果を取得できます。
待ち時間に次の「対象分子」へ進んでください。'''),
(C,'''RUN_BELL = False
bell_record_path = OUTPUT_DIR / "bell_submission.json"
if RUN_BELL:
    if service is None or not BACKEND_NAME:
        raise ValueError("IBM接続とBACKEND_NAMEの指定が必要です")
    if bell_record_path.exists():
        raise FileExistsError("既存ジョブを先に確認してください")
    from qiskit_ibm_runtime import SamplerV2
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
    backend = service.backend(BACKEND_NAME)
    measured_bell = bell.copy(); measured_bell.measure_all()
    isa = generate_preset_pass_manager(backend=backend, optimization_level=3, seed_transpiler=11).run(measured_bell)
    job = SamplerV2(mode=backend).run([isa], shots=1024)
    write_json(bell_record_path, {"job_id":job.job_id(), "backend":BACKEND_NAME})
    print("Job ID:", job.job_id())
else:
    print("Bell回路の実機送信はスキップ")'''),
(C,'''FETCH_BELL = False
if FETCH_BELL:
    if service is None: raise ValueError("IBMへ接続してください")
    bell_record = read_json(bell_record_path)
    saved_job = service.job(bell_record["job_id"])
    print("Status:", saved_job.status())
    if "DONE" in str(saved_job.status()).upper():
        print(saved_job.result()[0].data.meas.get_counts())'''),
(M,'''## 今日の分子を確認する
ブタジエンを共通の主題材にし、実機演習の初期設定は5量子ビットのエチレンとします。
9量子ビットのブタジエンは同じAPIで切り替えられますが、現在の回路は深く、講師による事前検証を条件にします。

計算条件は固定された教材用構造・STO-3G・π活性空間のCASCIです。
これは実験値を予測する精密な分子計算ではありません。'''),
(C,'''for molecule in ["ethylene", "butadiene"]:
    meta, _ = load_model(molecule)
    print(molecule, "CAS", (meta["active_electrons"], meta["active_spatial_orbitals"]),
          "total qubits", meta["total_qpde_qubits"], "S0-T1 gap/eV", round(meta["gap_ev"], 4))
fig, axes = plt.subplots(1,2,figsize=(10,3))
for name, ax in zip(["ethylene","butadiene"],axes): plot_molecule(name, ax)
plt.show()'''),
(M,'''**準備完了の確認**：回路を表示できた／理想確率と測定回数を区別できた／チームの実機利用担当が決まった／対象分子と計算する2状態を説明できた。

次は `01_circuits_and_ibm.ipynb` を開きます。認証は各Colabランタイムで必要です。APIキーはノート間でコピーせず同じシークレットから読み込みます。''')
], '準備：Colab・基本回路・IBM接続')

write('01_circuits_and_ibm.ipynb',[
(M,'''# 演習1：分子の回路を比較し、IBM実機へ送信する

**問い**：状態や時間発展を近似すると、操作数と精度はどう変わるでしょうか。
MPS/MPOの知識は必要ありません。配布済み回路を使います。


回路条件（compact / reference）とチーム名は、当日の案内に従って設定してください。'''),
(C,BOOTSTRAP),
(C,'''MOLECULE = "ethylene"  # 主題材の拡張: "butadiene"
VARIANT = "compact"    # 比較用: "reference"
TEAM = "team01"       # 各チームで変更
SHOTS = 1024
meta, arrays = load_model(MOLECULE)
summary = circuit_summary(MOLECULE)
print("対象:", MOLECULE, "基底:", meta["basis"], "CAS:", (meta["active_electrons"], meta["active_spatial_orbitals"]))
print("一重項・三重項の S(S+1):", meta["s_squared"])
plot_molecule(MOLECULE); plt.show()'''),
(M,r'''## 状態準備：どんな状態を作るか
補助量子ビットと分子の状態を組み合わせ、
$(|0\rangle|S_0\rangle+|1\rangle|T_1,M_S=0\rangle)/\sqrt{2}$ を目標とします。
補助量子ビットは最大の番号、分子の系はq0から並びます。

referenceは小さな模型の厳密状態を汎用的に分解した準備回路です。
compactはRYと近接CNOTの浅い形に合わせた近似回路です。忠実度1なら目標状態と一致します。

**見ること**：操作数だけでなく、目標状態への近さを比較してください。'''),
(C,'''from qiskit.quantum_info import Statevector
joint_target = np.concatenate([arrays["singlet"], arrays["triplet"]])/np.sqrt(2)
for kind in ["reference", "compact"]:
    prep = load_circuit(MOLECULE, f"prep_{kind}")
    fidelity = abs(np.vdot(joint_target, Statevector.from_instruction(prep).data))**2
    cost = summary["circuits"][f"prep_{kind}"]
    print(kind, "fidelity=", round(fidelity,6), "CX=",cost["cx_count"], "depth=",cost["depth"])
# 回路全体は長いため、最初の12操作だけを表示します。
from qiskit import QuantumCircuit
preview = QuantumCircuit(meta["total_qpde_qubits"])
original = load_circuit(MOLECULE,"prep_compact")
for instruction in original.data[:12]:
    # 配布回路内のビット番号を明示的に引き継ぎます。
    indices = [original.find_bit(q).index for q in instruction.qubits]
    preview.append(instruction.operation, indices)
show_circuit(preview)'''),
(M,'''## 時間発展：どこまで短くできるか

referenceは6分割の一次Trotter回路、compactは1分割の短いPauli積回路の時間係数を事前に調整したものです。
ここでいう短縮は論文のMPO回路最適化そのものではなく、分子で精度と操作数の関係を学ぶための例です。

両者を、別途作った密な行列指数の厳密時間発展と比較します。referenceも有限分割の近似です。
表示するCX数と深さは全結合・共通基底での値で、IBM実機への配置後には変わります。'''),
(C,'''for key in ["evol_reference", "evol_compact"]:
    cost = summary["circuits"][key]
    error_key = "reference_relative_operator_error" if key.endswith("reference") else "relative_operator_error"
    print(key, "CX=",cost["cx_count"], "depth=",cost["depth"],
          "relative error=",summary["evolution"][error_key])
fig, axes = plt.subplots(1,2,figsize=(8,3))
keys=["evol_reference","evol_compact"]
axes[0].bar(["reference","compact"],[summary["circuits"][k]["cx_count"] for k in keys])
axes[0].set_ylabel("CX count (all-to-all)")
axes[1].bar(["reference","compact"],[summary["evolution"]["reference_relative_operator_error"],summary["evolution"]["relative_operator_error"]])
axes[1].set_ylabel("Relative operator error")
plt.tight_layout(); plt.show()'''),
(M,'''## 同じ回路を位相差推定へ接続する

準備 → 系の時間発展を反復 → 補助の位相操作 → 準備の逆回路 → 全ビット測定、の順です。
**配布された浅い回路を、そのまま連結します。** 密な行列に戻して汎用ゲートへ再分解することはしません。

まずノイズなしで計算します。ショットノイズと実機ノイズは次の演習で分けて考えます。
基準エネルギー差は正解合わせ用にのみ表示し、測定計画はそれを中心に設定していません。'''),
(C,'''epsilon_values = np.linspace(0,.5,17)
steps = 4
probs = [zero_probability(MOLECULE,eps,steps,VARIANT) for eps in epsilon_values]
fig, ax = plt.subplots(figsize=(6,3))
ax.plot(epsilon_values*HARTREE_TO_EV,probs,"o-")
ax.axvline(meta["gap_ev"],color="gray",ls=":",label="Model reference")
ax.set(xlabel="Trial gap epsilon (eV)",ylabel="P(all zero)")
ax.legend(); plt.show()'''),
(M,'''## 実機投入計画

同じ分子の測定を3ラウンド行います。初期設定は各ラウンド9点、1点1024ショットで、合計27回路・27,648ショットです。
ここでいうラウンドは時間発展の長さが異なる測定セットで、探索点は事前に固定しています。
結果を見て次の探索点を変える適応的な制御は、発展課題です。

当日の実機の状態と講師の指示で、まず1ラウンドだけに短縮できます。
エチレンなら `step_counts=[2]` と指定します。時間の単位は原子単位です。
ブタジエンは既定の計画でも非常に深くなるので、講師の検証前には実機投入しないでください。'''),
(C,'''plan = make_plan(MOLECULE, VARIANT, TEAM, SHOTS)
print("Circuits:",len(plan["rows"]),"Shots:",sum(r["shots"] for r in plan["rows"]))
for iteration in sorted({r["iteration"] for r in plan["rows"]}):
    row=next(r for r in plan["rows"] if r["iteration"]==iteration)
    native_cx=2*summary["circuits"][f"prep_{VARIANT}"]["cx_count"]+row["steps"]*summary["circuits"][f"evol_{VARIANT}"]["cx_count"]
    print("Round",iteration,"time=",row["time_atomic_units"],"steps=",row["steps"],"CX before routing ~",native_cx)
write_json(OUTPUT_DIR/f"{TEAM}_plan.json",plan)'''),
(C,AUTH),
(M,'''### 配置後の深さを確認（送信しない）

講師と接続先・利用枠を確認し、代表回路をIBM実機の命令セットへ変換します。
長い回路で信号が得られるとは限りません。ここでの数値を見て測定計画を短縮します。'''),
(C,'''PREVIEW_ISA = False
if PREVIEW_ISA:
    if service is None or not BACKEND_NAME: raise ValueError("IBM接続と実機名が必要です")
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
    backend=service.backend(BACKEND_NAME)
    pm=generate_preset_pass_manager(backend=backend,optimization_level=3,seed_transpiler=11)
    for iteration in sorted({r["iteration"] for r in plan["rows"]}):
        row=next(r for r in plan["rows"] if r["iteration"]==iteration)
        isa=pm.run(qpde_circuit(MOLECULE,row["epsilon_hartree"],row["steps"],VARIANT))
        print("round",iteration,"ISA depth",isa.depth(),"2-qubit operations",sum(len(op.qubits)==2 for op in isa.data))'''),
(M,'''### チーム代表者：送信する

以下は実機利用枠を消費します。講師の合図で `RUN_IBM = True` に変更します。
各ラウンドの送信直後にジョブIDを記録します。途中で止まった場合は、保存済み記録を先に確認してください。
再実行による重複送信を防ぐため、同じ記録ファイルがあると送信を停止します。'''),
(C,'''RUN_IBM = False
record_path = OUTPUT_DIR / f"{TEAM}_submission.json"
if RUN_IBM:
    if service is None or not BACKEND_NAME: raise ValueError("IBM接続と実機名が必要です")
    if MOLECULE == "butadiene":
        raise ValueError("ブタジエン実機版は講師の事前ベンチマーク後に有効化してください")
    submission = submit_ibm(plan,service,BACKEND_NAME,record_path)
    print("Saved:",record_path)
    print("Job IDs:",[job["job_id"] for job in submission["jobs"]])
else:
    print("実機送信なし。次の演習は配布の模擬測定データでも進められます。")'''),
(M,'''### 記録ファイルを持ち出す
Colabのファイル欄から `outputs/..._submission.json` をダウンロードしてください。
APIキーはこの記録に含みません。次の演習ノートでアップロードし、再認証して結果を取得します。
ジョブID・実機名・回路条件を残すことは、再現性のための作業です。

**討論（休憩前）**：ゲート数が小さい回路は常に優れているでしょうか。準備と時間発展のどちらの誤差が結果へ影響するでしょうか。''')
], '演習1：分子の回路比較・IBM実機投入')

write('02_inference.ipynb',[
(M,'''# 演習2：測定から分子のエネルギー差を推定する

**できるようになること**：測定結果を取得する／推定を更新する／不確かさと限界を説明する。


配布データは**回路から生成した模擬測定データ**で、IBM実機の測定結果ではありません。グラフでも出所を明記します。
実機結果が弱い・間に合わない場合も、このデータで解析を練習できます。'''),
(C,BOOTSTRAP),
(M,'''## A. 配布データで始める
`compact_ideal`は回路近似と有限ショットを含み、実機ノイズは含みません。
`compact_depolarized`はさらに教材用の一様な脱分極ノイズを加えたもので、IBM実機の詳細なノイズ模型ではありません。'''),
(C,'''EXAMPLE_FILE = "ethylene_compact_depolarized.json"
data = validate_counts(read_json(WORKSHOP_ROOT / "data/example_results" / EXAMPLE_FILE))
meta, arrays = load_model(data["molecule"])
print("DATA SOURCE:",data["source"])
print("Molecule:",data["molecule"],"variant:",data["variant"])
print("Noise model:",data.get("noise_model", "Actual raw counts; not a simulated noise model"))'''),
(M,'''## B. 自分のIBM実機結果を使う場合（任意）
配布データで進める場合は、この節のフラグをFalseのままにします。
演習1で保存した `..._submission.json` をアップロードし、下の `RECORD_FILE` にそのパスを指定します。
Colab左側のファイル欄からアップロードできます。実機で計算が完了していない場合は、待ってから結果取得だけを再実行します。'''),
(C,AUTH),
(C,'''USE_IBM_RESULTS = False
RECORD_FILE = "outputs/team01_submission.json"
if USE_IBM_RESULTS:
    if service is None: raise ValueError("IBMへ再接続してください")
    record=read_json(RECORD_FILE)
    data=fetch_ibm(record,service)
    write_json(OUTPUT_DIR/f"{data['team']}_results.json",data)
    meta, arrays=load_model(data["molecule"])
    print("DATA SOURCE:",data["source"],"backend:",data["backend"])'''),
(M,'''## 測定したものを確認する

各条件で、全ビットが0だった回数kと総ショット数nを取り出します。観測確率はk/nです。
成功・失敗の集計は実際の整数回数を使います。誤差緩和後の擬似確率を整数回数として扱いません。
Fire Opalの結果を取り込むには、その出力の統計的な扱いを別途設計する必要があります。'''),
(C,'''fig, ax=plt.subplots(figsize=(7,3.5))
for iteration in sorted({r["iteration"] for r in data["rows"]}):
    rows=[r for r in data["rows"] if r["iteration"]==iteration]
    x=np.array([r["epsilon_hartree"] for r in rows])*HARTREE_TO_EV
    p=np.array([r["successes"]/r["shots"] for r in rows])
    ax.plot(x,p,"o-",label=f"round {iteration}")
ax.set(xlabel="Trial gap epsilon (eV)",ylabel="Observed P(all zero)",title="SOURCE: "+data["source"],ylim=(-.02,1.02))
ax.legend(); plt.tight_layout(); plt.show()'''),
(M,r'''## 推定を更新する

求めたい値を $\Delta E$ とします。理想的には
$p=(1+\cos[(\Delta E-\varepsilon)t])/2$。
仮の $\Delta E$ ごとに、実際の回数kがどのくらい起こりやすいかを二項分布で計算します。

**事後分布 ∝ 事前分布 × 測定結果の尤度**。
今回は0〜0.5 Hartreeの一様な事前分布を使います。この範囲には別の周期の解が入りにくい時間設定を使います。
参照値を中心に事前分布を狭めることはしません。

実機の信号低下に対応するため、背景とコントラストも不明な量として候補を置き、データと共に更新してから平均化します。
モデルは $p=b+v(1+\cos[(\Delta E-\varepsilon)t])/2$ です。
背景bとコントラストvが全ラウンドで共通という簡略化を含みます。実機では回路の深さで変わるため、適合を確認する必要があります。

表示区間は**この事前分布と尤度モデルを仮定した95%ベイズ信用区間**です。模型・基底・回路近似の誤差をすべて含む区間ではありません。'''),
(C,'''history=infer_gap(data,noise_aware=True)
fig=plot_posterior(history,reference=meta["gap_hartree"])
fig.axes[0].set_title("SOURCE: "+data["source"])
plt.show()
for row in history:
    low,high=np.array(row["ci95_hartree"])*HARTREE_TO_EV
    print(f"Round {row['iteration']}: {row['mean_hartree']*HARTREE_TO_EV:.3f} eV, conditional 95% interval [{low:.3f}, {high:.3f}]")
    if row["edge_mass"]>.1: print("探索範囲の端に質量があります。探索範囲とモデルを確認してください。")
print("Classical reference for this CASCI model:",meta["gap_ev"],"eV")
diagnostics=inference_diagnostics(data,history)
print("Scan variation detected:",diagnostics["detectable_scan_variation"])
print("Exploratory model check:",diagnostics["approx_model_tail_probability"])
if not diagnostics["detectable_scan_variation"]:
    print("信号を確認できません。狭い推定区間を成功と判断せず、測定データを再点検してください。")
if diagnostics["model_mismatch_warning"]:
    print("共通の背景・コントラストという尤度モデルとのずれが疑われます。区間を過信しないでください。")'''),
(M,'''## 誤差を3つに分ける

| 要因 | 比較の仕方 | ショットを増やせば消えるか |
|---|---|---|
| 有限ショットの統計変動 | 同じ回路・同じノイズ模型で測定回数だけ変える | 小さくなる |
| 回路近似による偏り | ノイズなしのreference回路とcompact回路を比較 | 消えない |
| 実機ノイズ | 同じ回路の理想計算と実機の測定を比較 | 一般には消えない |

以下では配布の模擬測定を使ってノイズを無視した推定と比較します。
**区間が狭いことと、正しいことは別です。** 信号が消えた結果を、無理に成功した推定へ変換しないでください。'''),
(C,'''simple=infer_gap(data,noise_aware=False)[-1]
aware=history[-1]
for title,row in [("Ideal likelihood",simple),("Background/contrast marginalized",aware)]:
    print(title,"mean/eV",round(row["mean_hartree"]*HARTREE_TO_EV,3),
          "interval/eV",np.round(np.array(row["ci95_hartree"])*HARTREE_TO_EV,3))
# 生データを見てから区間を解釈します。複雑なノイズや近似誤差まで補正できたことにはなりません。
summary_result={"source":data["source"],"molecule":data["molecule"],"variant":data["variant"],
                "mean_hartree":aware["mean_hartree"],"conditional_ci95_hartree":aware["ci95_hartree"],
                "likelihood":"binomial with shared background/visibility nuisance grid",
                "limitations":"conditional on likelihood and scan range; excludes molecular model error"}
write_json(OUTPUT_DIR/"inference_summary.json",summary_result)'''),
(M,'''### ショット数を変える（時間があるチーム向け）
教材の同じ回路とノイズ模型で、ショット数だけを変えます。
同じ一回の模擬実験で、区間幅が必ず単調に減ると保証されるわけではありません。'''),
(C,'''RUN_SHOT_COMPARISON = False
if RUN_SHOT_COMPARISON:
    for shots in [256,1024,4096]:
        simulated=simulate_results(make_plan("ethylene",shots=shots),seed=19,noise_rate=.15)
        last=infer_gap(simulated)[-1]
        print("Shots per circuit:",shots,"interval width/eV:",np.diff(last["ci95_hartree"])[0]*HARTREE_TO_EV)'''),
(M,'''## チーム討論

1. 何の分子の、どの2状態の差を求めたか。実験値と比べる前に何を確認するか。
2. 短い回路と長い回路で、どの指標が良くなり、どの指標が悪くなったか。
3. 事後分布が狭くなったら成功と言えるか。基準値からのずれと整合するか。
4. 実機の信号が弱かった場合、回路を浅くする・時間を短くする・ショットを増やす、どれを試すか。
5. 活性空間と基底関数を拡大したとき、今日の結果の意味はどう変わるか。

各組は「対象」「実行条件」「分かったこと」「残る課題」を1分で共有してください。

次に詳しく学びたい方は `03_mps_mpo_advanced.ipynb` へ進んでください。''')
], '演習2：結果取得・推定更新・誤差評価')

# Preserve the working v1 numerical cells, add an explicit molecular bridge.
legacy=nb.read(ROOT.parent/'versions/v1_2026-11-04/01_mps_mpo_circuits.ipynb',4)
advanced=[nb.v4.new_markdown_cell('''# 発展：MPS・MPOと分子の量子回路を理解する

本編を終えた方への任意教材です。SVD・ボンド次元・回路生成の数学を扱います。
小さな量子系の例で原理を確かめた後、分子データに適用します。
本編の進行に、このノートの理解は必須ではありません。

**配列順の注意**：第一版の関数は回路図の上からの順へ変換します。
新版の分子配列はQiskitの標準順です。分子をMPSへ渡すときは下記で明示的に軸を反転します。
また、最終の第一版QPDE例は密なStatePreparation/UnitaryGateによる原理確認で、本編の実機回路ではありません。'''),nb.v4.new_code_cell(BOOTSTRAP)]
for c in legacy.cells:
    if c.cell_type=='code':
        c.source=c.source.replace('from demo_core import *','from tensor_demo import *')
        c.execution_count=None;c.outputs=[]
    if c.cell_type == 'markdown':
        from participant_text import clean_legacy_markdown
        c.source = clean_legacy_markdown(c.source)
    c.metadata['audience'] = 'participant'
    advanced.append(c)
advanced.extend([
nb.v4.new_markdown_cell('''## 今回の分子とMPS/MPOを対応させる

本編の対象分子について、厳密な分子状態からMPSを作り、打ち切りの影響を確認します。
これは密な配列を使う教育用の処理で、大規模化できる実装ではありません。

本編の配布回路はRY/CXとPauli積の最適化で作ったもので、このMPS/MPOを使うテンソル最適化で生成したものではありません。
論文ではテンソルネットワークを使い、軌道局在化・並べ替えも工夫しています。'''),
nb.v4.new_code_cell('''import workshop as molecular_workshop
molecular_meta, molecular_arrays=molecular_workshop.load_model("ethylene")
n=molecular_meta["system_qubits"]
# Qiskit配列(q[n-1]...q0)を、MPSで扱うq0,q1,...の軸順に変換。
state_topdown=molecular_arrays["singlet"].reshape((2,)*n).transpose(tuple(reversed(range(n)))).reshape(-1)
for chi in [1,2,4]:
    mps,_=state_to_mps(state_topdown,n,max_bond=chi)
    reconstruction=mps_to_state(mps)
    reconstruction/=np.linalg.norm(reconstruction)
    print("chi",chi,"fidelity",abs(np.vdot(state_topdown,reconstruction))**2,"shapes",[a.shape for a in mps])
# 分子の小さな厳密時間発展もMPOへ分解できる。
U_native=expm(-1j*molecular_arrays["hamiltonian"]*.5)
axes=list(reversed(range(n)))+list(reversed(range(n,2*n)))
U_topdown=U_native.reshape((2,)*(2*n)).transpose(axes).reshape(2**n,2**n)
for chi in [1,2,4,16]:
    mpo,_=operator_to_mpo(U_topdown,n,max_bond=chi)
    reconstructed=mpo_to_operator(mpo)
    print("MPO chi",chi,"relative error",relative_error(reconstructed,U_topdown),"unitarity error",unitarity_error(reconstructed))''')])
adv=nb.v4.new_notebook(cells=advanced,metadata={'audience':'participant','kernelspec':{'display_name':'Python 3','language':'python','name':'python3'},'language_info':{'name':'python'}})
nb.write(adv,ROOT/'03_mps_mpo_advanced.ipynb')
print('Created overview and four participant notebooks.')
