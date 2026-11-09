# 演習01〜09

行列 → H₂の1 qubit模型 → H₂の4 qubit JW表現 → H₄のPEC → Gaussian更新 → 実機デモ → 回路長の実機比較の順に進みます。
各冊は独立したNotebookです。前の冊の変数・出力ファイルは必要ありません。

| 順番 | 教材 | 内容 | Google Colab |
|---|---|---|---|
| 01 | [行列の位相推定](01_matrix_phase_estimation.ipynb) | QPE・BPE・BPDE、Gaussian更新 | [開く](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/exercises/01_matrix_phase_estimation.ipynb) |
| 02 | [H₂・1 qubit](02_h2_one_qubit.ipynb) | 電子配置とスピン、S₀–S₂の位相差 | [開く](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/exercises/02_h2_one_qubit.ipynb) |
| 03 | [H₂・Jordan–Wigner・4 qubit](03_h2_jw_bpde.ipynb) | 02と同じ流れでS₀–T₀の差を測定 | [開く](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/exercises/03_h2_jw_bpde.ipynb) |
| 04 | [H₄のポテンシャル曲線](04_h4_circuit_assessment.ipynb) | 保存回路の表、期待値のPEC、Lを増やす位相実験 | [開く](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/exercises/04_h4_circuit_assessment.ipynb) |
| 05 | [H₄のGaussian更新](05_h4_iterative_bayes.ipynb) | シミュレーションで測定・フィット・事後更新を反復 | [開く](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/exercises/05_h4_iterative_bayes.ipynb) |
| 06 | [H₄の実機デモ](06_h4_hardware_demo.ipynb) | 同じ9点で理想確率・有限shots・IBM実機を比較 | [開く](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/exercises/06_h4_hardware_demo.ipynb) |
| 07 | [H₄の実機信号と回路長](07_h4_hardware_depth_scan.ipynb) | L=1・2・4の位相走査で信号の減衰を調べる | [開く](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/exercises/07_h4_hardware_depth_scan.ipynb) |
| 08 | [H₄のjoint-state圧縮と誤り抑制](08_h4_joint_compression_suppression.ipynb) | joint準備・較正配置・twirling・DDをfakeと実機で比較 | [開く](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/exercises/08_h4_joint_compression_suppression.ipynb) |
| 09 | [H₄のpaired/joint準備比較](09_h4_joint_preparation_demo.ipynb) | 06のpaired準備とjoint準備をL=1・位相一周期で比較 | [開く](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/exercises/09_h4_joint_preparation_demo.ipynb) |

## 実行方法

Colabでは上の「開く」から開始し、自分のコピーを保存してください。
初回は冒頭の `%pip` 行のコメントを外して実行し、再起動を求められたらランタイムを再起動して上から実行します。

ローカルでは次のように起動します。

```sh
python -m pip install -r hands_on/exercises/requirements.txt
jupyter lab hands_on/exercises
```

PySCF・Juliaは不要です。08・09および06のfakeモードでは `qiskit-ibm-runtime` と `qiskit-aer` を導入します。
08・09は冒頭のインストール行を使ってください。06・07の実機利用時は接続準備セルでRuntimeも導入します。

## 冊ごとの変数と時間

| 演習 | 系 | 時間発展 | 推定する差 |
|---|---|---|---|
| 01 | 1 qubit | U=exp(+2πiK) | 位相差、1周単位 |
| 02 | H₂の閉殻2配置、1 qubit | 厳密発展、Δt=1 a.u. | S₂−S₀、コードGAPは位相差 |
| 03 | H₂/JW、4 qubit | 二次Trotter・分割数1、Δt=0.1 a.u. | T₀−S₀、コードGAPは位相差 |
| 04 | H₄/JW、8 qubit | 保存発展または厳密発展、Δt=0.1 a.u. | S₀・T₀のPECと正のエネルギー差 |
| 05・06 | H₄/JW、8＋補助1 qubit | Δt=0.1 a.u. | 正のT₀−S₀差、Ha単位 |
| 07 | H₄/JW、8＋補助1 qubit | Δt=0.1 a.u.、L=1・2・4 | 位相走査の信号検出、ギャップ推定なし |

物理時間は t=LΔt、分子の固有位相は φ=−E_elΔt/(2π) mod 1です。
位相候補ではβ=−2πLε、電子エネルギー／エネルギー差の候補ではβ=ε_E tを使います。
01〜03のGaussian更新の関数・表示形式は共通です。03は02のJW版で、事前分布とQPE分解能をΔtに合わせています。
Gaussianの幅は、その近似更新の幅であり、全誤差を覆う保証ではありません。

## 04：期待値とL依存のPEC

- `EXPECT_PREP_M`：4・6・8・`'exact'`。Rの20点で保存準備状態の〈H〉を計算します。
- `PHASE_PREP_M`：4・6・8・`'exact'`。
- `PHASE_EVOLUTION_M`：8・16・32・`'exact'`。
- `L_MAX`：L=1,2,3,…,L_MAXを計算。既定64、最大1024。
- `L_SHOW`：PECの図に重ねるL。計算するLを間引く設定ではありません。

S₀のBPE位相とBPDEの差から、S₀・T₀の全エネルギーを求めます。核反発は両PECに加え、差では相殺します。
20点のPEC、FCIとの差、全LのRMSE、最終Lの数値表を表示します。
長時間の位相の枝は補助時間を使うシミュレーションで確認します。FCIや〈H〉で正解の枝を選びません。
末尾に、短時間では位相が期待値を読み、目的状態が優勢なら長時間の誤差上限が1/tで減る説明があります。
結果は `outputs/04_expectation_pec.json` と `outputs/04_phase_pec.json` に保存します。

## 05：Gaussian更新の反復

構造 `GEOMETRY`、準備 `PREP_M`、発展 `EVOLUTION_M` を選びます。
`MODE='exact'`はshot noiseなし、`MODE='shots'`は有限shotsの模擬測定です。
`CASE='exact'/'prep_only'/'evolution_only'/'both'`で準備・発展の近似を切り替えます。
MODEとCASEは別の設定で、MODE=exactでも回路近似は残り得ます。

`ROUNDS`、初期事前、`POINTS`、`MAX_L`を指定し、測定条件の選択 → 測定 → Gaussianフィット → 事後更新を繰り返します。
左に初期・途中・最終の分布、右に最終分布の拡大、下に更新表を表示します。
末尾のCASE比較とM比較はshot noiseなしのシミュレーションです。
結果は `outputs/exercise5_updates.json`・`outputs/exercise5_measurements.json` に保存します。

## 06：同じ9点のシミュレータ・実機比較

R=2 Å、準備M6、発展M16、L=4、候補9点、初期事前μ=0・v=4に固定し、1回だけGaussian更新します。
各点のshots数は `SHOTS` で指定できます。
理想確率、有限shotsの模擬測定、実機を同じ条件で比較します。

実機を使う場合は以下の順です。

1. [登録](../../docs/01_ibm_quantum_registration.md)と[接続確認](../../docs/02_colab_connection.md)を済ませる。
2. 06のIBM接続準備セルを `USE_IBM=True` にし、本人のAPIキーを非表示入力してCRN・backendを選ぶ。
3. L=4の回路プレビューを実行し、配置後の回路費用を確認する。
4. 送信する場合だけ `RUN_HARDWARE=True` にし、実機測定セルを実行する。

初期設定では接続・送信を行いません。ジョブID・counts・更新結果は `outputs/exercise6/` に保存します。
同条件の再実行では保存済みジョブを利用します。Colabの切断に備え、必要な記録をDrive等へ保存してください。

06の第3a節では、`qiskit-aer` とRuntimeを導入して `USE_FAKE_KAWASAKI=True` にすると、
APIキーなしでFakeKawasakiのノイズ付き模擬測定を実行できます。この間は後のIBM接続・実機送信をスキップします。
127 qubit・ECRゲートの過去のスナップショットを使うため、現在のCZ装置の較正とは異なります。
全countsとスナップショット情報は `outputs/exercise6/fake_kawasaki_measurements.json` に保存します。

## 07：回路長と実機信号の比較

06と同じH₄・準備M6・発展M16で、L=1・2・4を比較します。位相角を同じ一周期の12点に揃え、
短い時間による理想曲線の平坦化と、実機ノイズによる干渉の減衰を区別します。
各点2回、既定1024 shots、合計72回路・73,728 shotsを1ジョブにまとめます。

`USE_IBM=True`で接続し、同じ初期配置で各Lの回路費用をプレビューした後、`RUN_HARDWARE=True`で送信します。
既定ではシミュレーションのみです。生counts・ジョブID・送信したISA回路を `outputs/exercise7/ibm_jobs/` に保存し、
同条件の再実行では保存結果を利用します。独立に測定し直す場合は `RUN_LABEL` を変えてください。

全0頻度の95%区間と、一定確率に対する位相依存性の検査を表示します。Gaussian更新やギャップ推定は行いません。
保存JSONを `RESULT_FILE` に指定すると、再接続せずに解析できます。
[07のローカル実行・模擬ジョブ検証](../../docs/validation_20261109_depth_scan.json)では、実機へ送信していません。

## 08：joint-state圧縮とhardware-awareな誤り抑制

06のjoint-stateをMPSへ変換し、分岐制御を使わない9-qubitの局所回路へ直接圧縮します。
paired準備、joint準備、joint＋較正配置・level 3・Pauli twirling・XX-DDの3条件を比較します。
FireOpalの再実装ではなく、読み出し誤り緩和や独自パルス最適化は含みません。

既定はL=1、各armの各位相1024 shots、12位相・2回で、1 backendあたり合計36,864 shotsです。
L=4へ変えると06と同じ発展時間になります。joint標的とのfidelityが0.99未満なら送信を止めます。

FakeKawasakiではゲート・読出し・明示的DelayのT1/T2ノイズを含めて、既定でローカル実行します。
実機では `USE_IBM=True` で接続し、3条件の費用をプレビューした後 `RUN_HARDWARE=True` で実行します。
送信回路・ジョブID・全countsは `outputs/exercise8/` に保存します。実機接続・送信は既定で無効です。
[08の検証記録](../../docs/validation_20261109_joint_compression.json)にはfake実行と模擬実機ジョブの結果を記録しています。

## 09：paired準備とjoint準備を同条件で比較

06のM6 paired準備と、その出力を直接圧縮したjoint準備を比較します。
既定はL=1・位相一周期β∈[-π,π)の12点です。旧09の±0.4 radではピーク付近しか測れないため、ピークと谷を含めて信号を調べます。M16発展・初期物理qubit配置・level 1・seed=7・shotsを揃えます。
各条件・各点1024 shots、合計24回路・24,576 shotsです。Lを変えるときも両条件に共通です。候補はε=β/(LΔt)で計算し、一周期を保ちます。これは干渉信号の診断で、ギャップ精度を示す走査ではありません。
較正配置探索・level 3・twirling・DDは追加しません。

理想確率、shot noiseのみ、FakeKawasaki、実機で両条件のcounts・位相依存性・回路費用を比較します。
各条件のp値は2条件分を補正します。これは条件間の優劣の検定やギャップ精度の保証ではありません。
既定はFakeKawasakiです。実機では `USE_FAKE_KAWASAKI=False`、`USE_IBM=True` で接続し、
両条件のプレビュー後に `RUN_HARDWARE=True` にします。両条件を交互に並べて1ジョブへ送信します。
保存先は `outputs/exercise9/` です。同条件の保存ジョブを再利用し、06のジョブと混同しません。
[09の検証記録](../../docs/validation_20261109_full_phase_comparison.json)にはfake実行と模擬実機ジョブの検証を記録しています。

## H₄の固定データ

[全20構造の固定データ集](data/h4_pec_v1/README.md)を使います。近似準備はχ=6・M=4/6/8、近似発展はM=8/16/32です。
厳密発展はシミュレーション参照です。Mは回路の層数で、MPS/MPOのボンド次元とは別です。

04〜09は固定コミット `4ab6ef3d06f4256b5c818feb8c3fa2af4381dd43` を使い、`H4_DATA_DIR` → `H4_CACHE_DIR` → 取得済みリポジトリ → `~/.cache/ftqc_handson/h4_pec_v1` の順で保存先を決めます。
指定した場所に不足ファイルをダウンロードし、同じ場所を指定すれば再起動後も再利用します。04の従来の保存先を使う場合は `H4_CACHE_DIR = 'outputs/h4_pec_v1'` と指定できます。
各ファイルはSHA-256で検証して再利用します。Google Driveを使う場合は先にmountし、`H4_CACHE_DIR`へパスを指定してください。
05〜09でオフライン使用する範囲を先に取得するには `download_h4_data()`、取得禁止には `H4_OFFLINE=True` を使います。

## 開発・検証・旧版

正式版はNotebook本体を編集します。差し替え前の生成スクリプトはarchiveへ移しました。

```sh
python hands_on/exercises/tools/execute_notebooks.py
python -m pytest hands_on/exercises/tests -q
```

実行スクリプトは各冊を別のPython/IPythonプロセス・一時フォルダで実行し、出力をNotebookに保存します。
06も初期設定のまま実行するため、IBM実機へのジョブ送信を行いません。
実行記録と図はGit管理外の `outputs/qa/`、今回の記録は [正式版01〜06の検証](../../docs/validation_20261109_promotion.json) にあります。
これはローカル実行の検証です。ColabやIBM実機の検証とは区別します。

[差し替え前の01〜05と旧生成コード・テスト](../../archive/pre_promotion_2026-11-09/README.md)／[それ以前の教材](../../archive/README.md)。
元の自習教材は `archive/self_study/` にあります。固定データの入力一覧は [SOURCE_FILES.json](SOURCE_FILES.json) に記録しています。

## 参考研究

[Kanno et al., PNAS (2025)](https://doi.org/10.1073/pnas.2425026122)と[著者コード](https://github.com/sk888ks/TQPDE_open)を参照しています。
[論文PDF・補足資料・参考資料](../../docs/README.md)／[原著コードのLGPL-3.0ライセンス](../../archive/self_study/molecular_h4_jw/UPSTREAM_LICENSE)。
本教材は論文の分子・実機実験全体の再現ではありません。
