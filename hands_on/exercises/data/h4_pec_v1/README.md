# H4/JW 全20構造の教材用スナップショット v1

直鎖H4/STO-3G、4電子・4空間軌道、削減なしJW 8 qubit。
Rは隣接H原子間距離（Å）、20点は約0.56164〜3.0 Å。Boys局在化＋軌道順最適化。
モード順はq0..q7=alpha0..3,beta0..3、ketはq7..q0、BPDEの補助はq8です。
電子Hamiltonianで Δt=0.1 原子単位、U=exp(-iHΔt)。全エネルギーには核反発を加えます。

研究用ディレクトリから明示的に書き出した固定版です。通常のNotebook再生成・実行では最適化も研究用ファイルの読み込みも行いません。
各点の軌道・Hamiltonian・準備状態・発展を揃え、発展には符号付きFock順変換を含めています。
BPEでも使えるように回路の全体位相を保持しています。

## 内容

- `manifest.json`：20構造のR・参照エネルギー・ファイル一覧・SHA-256。
- `R_*/model.npz`：電子H、軌道係数、積分、核反発、FCI S0/S1/T0と電子エネルギー。
- `R_*/prep_{4,6,8,exact}_{S0,T,A}.qpy`：単独準備とBPDE両枝の9 qubit準備A。逆準備は同じAのinverse()。
- `R_*/evolution_{8,16,32,exact}.qpy`：一歩の発展。L回の反復は教材で組み立てます。
- `R_*/catalog.json`：fidelity、CX/depth、L=1信号の参照値。MPSのχは6、Mは回路層数です。
- `R_*/provenance.json`：元の研究用ファイルの相対パスとハッシュ。これらの元ファイルは受講者の実行には不要です。
- `exact_preparation.py`：厳密な実振幅準備・スピン演算子等の補助関数。

`exact`は数値精度での厳密参照です。厳密発展QPYは256×256の密行列unitaryであり、ネイティブゲートへ分解した回路ではありません。
そのCX数・depthは未評価（null）で、0ではありません。実機には近似発展M=8/16/32を選びます。
近似発展はcutoff=1e-18・二次Trotter100分割の参照MPOへの有限予算のフィットです。完全収束・大域最適性は保証していません。
CX値は全結合・rz/sx/x/cxでの部品単位の値です。実機配置と全回路の最適化は含みません。
保存済みの発展basic回路と元のブロック列の微小な合成誤差は、元のreportのbasic_operator_errorと区別して扱います。

## 受講者の使い方

このフォルダを含めたリポジトリを取得し、04/05を上から実行します。01〜03の実行状態は不要です。
初期化後、`list_geometries()`で構造IDを確認し、冒頭のセルで指定します。

```python
PREP_M = 'exact'       # 'exact', 4, 6, 8
EVOLUTION_M = 32       # 'exact', 8, 16, 32
GEOMETRY = 'R_3p00000000'
select_geometry(GEOMETRY)
select_circuits(PREP_M, EVOLUTION_M)
```

以降を上から再実行します。CASE='both'で指定部品を使います。CASE='exact'は選択構造の厳密準備＋厳密発展の対照です。
MODE='exact'はshot noiseなし、MODE='shots'は選んだshotsの二項模擬測定、MODE='ibm'は実機です。
通常は作業場所から自動検索します。別の場所へ置いた場合は、選択前に `os.environ['H4_DATA_DIR']` にこのフォルダの絶対パスを設定します。
既定R=2 ÅだけはNotebook内蔵データでも実行できます。別構造を選んでファイルが見つからない場合は停止します。

QPYはQiskit 2.4.1で書き出しています。教材のrequirements.txtを使用してください。PySCF・Juliaは受講者側では不要です。
個々の回路は `qiskit.qpy.load()`、模型は `numpy.load(..., allow_pickle=False)` でも読み込めます。

## GitHubで共有するファイル

`hands_on/exercises/data/h4_pec_v1/` をフォルダごと、04/05とtoolsの補助・生成コード、requirements.txtと一緒に登録します。
研究用の全ログ・最適化環境・アカウント情報を配布する必要はありません。実機の認証キー・ジョブ記録は含めていません。
授業ではデータとNotebookを含む同じコミットを使用します。今後の研究成果を反映するときはv2等の別版へ書き出し、v1を上書きしません。

教師用の初回書き出しは `tools/freeze_h4_pec.py` です。既存manifestがある場合は上書きせず停止します。
`manifest.json` は配布データのハッシュ、`tools/validate_h4_snapshot.py` は内容・回路・信号の検証に使用します。
出所は[Kanno et al., PNAS (2025)](https://doi.org/10.1073/pnas.2425026122)と固定版の著者手法に基づくこのプロジェクトの保存計算です。
著者コードのライセンスは同梱の `UPSTREAM_LICENSE` を参照してください。
