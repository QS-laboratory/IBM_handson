# 新本編：5つの演習

既存の自習教材を使い、行列→H₂の2状態部分空間→H₂/JW→H₄/JW→反復的な推定へ進みます。
旧本編は `archive/pre_revision_2026-11-08/` に保存し、元の自習教材・再生成用データは `archive/self_study/` に保持しています。授業の時間割はまだ割り当てていません。

| 順番 | Notebook | 内容 |
|---|---|---|
| 1 | [2×2行列の位相推定](01_matrix_phase_estimation.ipynb) | diag(1/4,3/4)、IQFT型QPE・BPE・BPDE |
| 2 | [H₂・1 qubit](02_h2_one_qubit.ipynb) | 化学への入口、S₀/S₁/S₂/T₀、S₀–S₂部分空間 |
| 3 | [H₂・JW・4 qubit](03_h2_jw_bpde.ipynb) | 厳密状態準備、BPDE、Trotter分解 |
| 4 | [H₄の回路評価](04_h4_circuit_assessment.ipynb) | 厳密準備/Trotterの費用、全20構造・厳密/近似選択、fidelity、費用、L=1曲線、実機接続 |
| 5 | [H₄の反復的Bayesian推定](05_h4_iterative_bayes.ipynb) | 測定条件の選択、事前→測定→事後を繰り返す |

各Notebookは前のNotebookの実行状態に依存しません。01〜03と04・05の既定R=2 Åはデータ・補助関数を内蔵します。
04・05の全20構造を使うときは、同じコミットの [固定データ集](data/h4_pec_v1/README.md) を含めてリポジトリを取得してください。
カーネルを再起動して上から実行します。必要なインストール例は冒頭にあります。
予想・課題・確認問題には、その場で開ける解答例を付けています。
演習1は模式回路、式、測定例、BPE／BPDE B／BPDE Cのベイズ更新を単独で確認できます。

冊をまたぐときは、位相差とエネルギー差、基本時間の変更に注意してください。

| 教材 | 基本操作・基本時間 | コードの差の変数 |
|---|---|---|
| 01 | `U=exp(+2πiK)` | `GAP`は位相差、1周単位 |
| 02 | `U=exp(-iH DT)`、`DT=1`原子単位 | `GAP`は周期1の位相差、`ENERGY_GAP`はHa単位 |
| 03 | 同じH₂、`DT=0.1` | `GAPS`は各状態への正のエネルギー差、Ha単位 |
| 04・05 | 同じH₄と保存回路、`DT=0.1` | `GAP`は正のT₀−S₀差、Ha単位 |

物理時間は`t=L DT`、分子の固有位相は`φ=-E_el DT/(2π) mod 1`です。
位相候補では`β=-2πL ε`、電子エネルギー／差の候補では`β=ε_E t`へ換算します。
QPEの補助数m、Trotter分割数sと、05のGaussianフィット中心m・幅sは文脈が異なります。

```sh
python -m pip install -r hands_on/exercises/requirements.txt
jupyter lab hands_on/exercises
```

ColabではルートREADMEのリンク、またはNotebookのアップロードから開けます。
検証済みコミットと条件は [5冊の初回検証](../../docs/colab_validation_1203538.json) と [追加改訂後の04・05の検証](../../docs/colab_validation_900224b.json) を参照してください。以降の改訂は記録ごとに検証範囲を区別します。IBM実機への送信はこの検証に含めません。
本編にPySCF・Julia・Aerは不要です。分子積分の再生成や最適化は元の自習教材を参照してください。

## 演習4・5の測定モード

構造ID `GEOMETRY` と、準備 `PREP_M='exact'/4/6/8`、発展 `EVOLUTION_M='exact'/8/16/32` を選べます。既定はR=2 Å・6/16です。
`list_geometries()`で全20構造を確認できます。設定後は選択セル以降を再実行してください。
厳密発展は密行列を持つシミュレーション参照回路で、ネイティブCX数は未評価です。IBMには発展M=8/16/32を選びます。
演習4・節5.2では9通りのstate fidelity・operator fidelity・回路費用・L=1信号を比較します。
演習5・節8.1では同じ初期事前から9通りを反復し、推定平均のずれ・SD・L・フィット残差を比較します。更新上限は `M_COMPARE_ROUNDS` で指定します。
これらの追加比較は常にshot noiseなしで、主実験の設定・履歴を保持し、実機には送信しません。Mは回路の層数で、MPS/MPOのボンド次元とは区別します。

- `MODE='exact'`：選択した回路の厳密確率。shot noiseなし。
- `MODE='shots'`：理想確率から二項分布で模擬測定。各探索点の`SHOTS`と`SEED`を指定。
- `MODE='ibm'`：IBM SamplerV2で全9ビット0の頻度を取得。`CASE='both'`を使用。

`CASE='exact'/'prep_only'/'evolution_only'/'both'`で、準備・発展の近似を切り替えます。
CASE=bothは指定した両部品を使い、CASE=exactは選択構造の厳密準備・厳密発展です。
MODEとCASEは独立です。MODE=exactでも、指定部品が近似ならその誤差は残ります。
模擬測定にゲートノイズや読み出しノイズは加えていません。

実機を使う場合だけ追加でインストールします。

```sh
python -m pip install 'qiskit-ibm-runtime>=0.43,<0.44'
```

APIキーはColabシークレット`IBM_QUANTUM_API_KEY`、環境変数、または非表示入力から取得します。
backend名と必要ならインスタンスCRNを指定し、`hardware_preview`で実機配置後の費用を確認してから
`RUN_HARDWARE=True`にします。初期値では送信しません。

演習4はL=1・9探索点。演習5は1ラウンド21探索点を初期値とし、結果を回収・更新してから次を送ります。
各ラウンドの事後を次の事前とするため、後続の実機ジョブは前の結果に依存します。
`MAX_HARDWARE_L=8`を初期上限とし、超過時は送信前に停止します。上限を増やす前に
対象Lの費用を確認してください。長い回路でも実機で信号が残るとは限りません。
演習5は上限までの更新を表示・保存して終了します。`RUN_HARDWARE=False`ならプレビューだけを行い、測定はスキップします。
MODE・CASE・SHOTSを変更したら設定セル以降を再実行し、測定器と乱数を初期化してください。

送信IDと生countsは`outputs/ibm_jobs/`へ保存し、同条件の再実行では既存ジョブを復旧します。
実行中断後は同じ設定で最初から再実行できます。ColabではこのフォルダをDrive等へ保存してください。
送信直後の中断でIDなしマーカーが残った場合は、自動再送せず、IBM側を確認するよう停止します。
別の実機試行を行う場合は、元の記録を保持したまま`Measurement(directory=...)`に別の保存先を指定します。

[IBM Sampler入出力](https://eu-de.quantum.cloud.ibm.com/docs/en/guides/sampler-input-output)

## 反復更新の範囲

演習5はラウンド数を2に固定しません。`ROUNDS`、初期分布、shots、探索点数、最大L、目標SDを変更できます。
演習4と5は同じ模型・保存回路・全0確率を使います。演習4はL=1の診断、演習5は独立した初期事前から測定条件を選ぶ反復です。
著者コードのGaussianフィット＋Gaussian積の更新を参照し、時間と探索幅もコードの分散ベースの規約に合わせます。
二項尤度によるグリッド事後分布と同じ手順ではありません。GaussianのSDは全誤差の保証ではありません。

既存のH₄自習Notebookは正準軌道ですが、新演習4・5は一致したBoys＋orderingの模型とM6/M16回路を使います。
S₀とTの両枝を近似し、実際の準備回路の逆を使います。符号付きFock順変換も発展回路に含めます。
全0確率は全shotsを分母とし、成功したショットだけで再正規化しません。
Notebook内の`PROVENANCE`に入力ファイルのSHA-256・元の最適化記録を保持します。

## 開発・検証

```sh
python hands_on/exercises/tools/build_notebooks.py
python hands_on/exercises/tools/execute_notebooks.py
python -m pytest hands_on/exercises/tests -q
```

生成元は`tools/build_notebooks.py`、各冊の本文と例は`tools/matrix_notebook.py`・`tools/h2_notebook.py`・`tools/h2_jw_notebook.py`・`tools/h4_assessment_notebook.py`・`tools/h4_bayes_notebook.py`、演習4・5の共通測定処理は`tools/h4_measurement.py`です。
演習2だけを再生成する場合は`python hands_on/exercises/tools/build_notebooks.py --only 02`を使えます。
再生成するとNotebookの保存出力は消えます。続けて実行してください。
演習5は初回を段階ごとに実行してから反復を初期値から再実行します。同条件の初回IBMジョブは再利用します。
各ラウンドの4図は計画・フィット・残差・Gaussian積を示します。末尾の回路近似比較は常にshot noiseなしのシミュレーションで、主実験の履歴を上書きしません。
実行器は各冊を別のPython/IPythonプロセス・一時作業ディレクトリで実行し、出力をNotebookへ保存します。
実行記録・QA画像・数値結果はGit管理外の`outputs/qa/`に保存します。
この検証はJupyterカーネル・Colab UI・実機の検証を代替しません。

## 出所

元の教材は[自習教材](../../archive/self_study/phase_estimation_comparison.ipynb)と、そのH₂/JW・H₄/JWの各冊です。
M6/M16回路は`archive/self_study/h4_mps_preparation_study/results/both_states_M6_evolution_M16_L1/`を使用します。
参考：[Kanno et al., PNAS (2025)](https://doi.org/10.1073/pnas.2425026122)、
[著者コード（固定コミット）](https://github.com/sk888ks/TQPDE_open/tree/9cd88464c5c50dda783150661813609092d7f2b6)、
[原著コードのLGPL-3.0ライセンス](../../archive/self_study/molecular_h4_jw/UPSTREAM_LICENSE)。
論文の分子・実機実験全体の再現ではありません。

## 全20構造の固定版データ

[データ集と読み込み方](data/h4_pec_v1/README.md) に、S₀・T₀の厳密/近似準備、BPDE用A、厳密/近似発展、模型・参照エネルギー・ハッシュをまとめています。
研究側の最適化を再実行せずに利用できます。ファイルを改訂するときは別バージョンへ書き出します。
全データの独立検証は `python hands_on/exercises/tools/validate_h4_snapshot.py` で実行できます。
