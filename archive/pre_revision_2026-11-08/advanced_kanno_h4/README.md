# 発展演習：著者のMPS/MPO回路最適化でH₄を計算する

既存の本編・第一版を変更せず、独立したフォルダに用意した発展演習です。

[H4_author_method.ipynb](H4_author_method.ipynb) を上から読んでください。

GitHub配置後のColabリンク：
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/advanced_kanno_h4/H4_author_method.ipynb)

ローカルの成果物をGitHubへ配置してから使えます。Colabの実環境での動作は未検証です。

## 今回「同じやり方」にした部分

Kannoらの [TQPDE_open](https://github.com/sk888ks/TQPDE_open) のJuliaコアをバイト同一で使用し、以下を実際に実行します。

1. DMRGで基底・第一励起状態を求め、補助付き重ね合わせMPSを作る。
2. Δt=0.1を100分割した二次Trotterから参照MPOを作る。
3. 環境テンソルとSVDで隣接2量子ビットbrick-wall回路を最適化する。
4. 得た4×4ゲートを、そのままQiskitの位相差推定回路へ連結する。
5. 著者と同じ時間選択・Gaussianフィット・Gaussian積の更新式で推定する。

前の教材のRY/CX角度フィットや、短いPauli積の時間係数フィットには置き換えていません。
ただし実行環境、DMRGの収束スイープ数、実行バックエンドには差があります。詳細は [UPSTREAM.md](UPSTREAM.md) に分けて明記しています。
IBM実機・Fire Opalを含むハードウェア実験の完全再現ではありません。

## 比較対象

- 4電子4空間軌道のH₄/STO-3G。系8＋補助1量子ビット。
- 補足Fig. S6の直鎖→正方形の5構造を復元し、全ての厳密なエネルギー差を表S1と照合。
- 全5構造について、正準軌道・Boys局在化・局在化後のGA順序最適化のHamiltonianを配布。
- 重いMPS/MPO回路最適化は、まず直鎖（ID 0）の3表現を同じ深さ8・1000スイープで比較。
- 他の構造の回路は再生成できます。正方形は対称性や初期値でBoys局在化が停滞する場合があり、必ず改善するとは仮定しません。

厳密解は軌道変換に対して不変ですが、固定深さの回路近似は変わります。
Pauli項数・交換積分コスト・MPSの重なり・回路近似誤差を別々に評価します。
推定平均が厳密値からずれる場合も、そのまま記録して議論します。

## 参加者向け：保存済み結果を比較する

NotebookをColabまたはJupyterで開き、上から実行してください。保存済みの数値と図も読めます。
重いJulia最適化は最初から再実行せず、配布されたNumPyゲート・最適化履歴・推定履歴を読み込みます。
受講者版には再生成セルを含めません。APIキーも不要で、QPU時間は使いません。
進行・解説は [講師用ガイド](instructor/GUIDE.md)、再生成は [講師用ノート](instructor/H4_regeneration.ipynb) に分離しています。

本編の追加30〜60分程度の読み解き、または別日の上級演習を想定しています。再生成の時間は含みません。

## 講師向け：全て再生成する

Python依存は本編の環境と分離して準備することを推奨します。

```sh
cd advanced_kanno_h4
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

著者のソースは `upstream/` に取得済みです。取得元コミットとライセンスは `UPSTREAM.md` を参照してください。
`runtime/Gate_generation_from_MPS_MPO.jl` はその未変更コピー、`runtime/input.jl` は保存先を与えるアダプタです。

分子前処理：

```sh
python src/h4_models.py 0 1 2 3 4
```

Julia 1.11.1を用意したうえで、専用のJulia環境を構築します。
`JULIA_DEPOT_PATH`はこのフォルダ内に置き、既存のJuliaプロジェクトを変更しません。

```sh
bash tools/setup_julia.sh
bash tools/run_author.sh
```

既定値はID 0、準備・時間発展とも深さ8、最適化1000スイープ、DMRG 200スイープです。
DMRGの収束は独立した厳密対角化と比較し、不一致なら圧縮前に停止します。

次に回路の整合性検証と位相差推定を行います。

```sh
python tools/analyze_all.py
python tools/build_notebook.py
python tools/run_notebook.py
python -m pytest tests -q
python tools/check_preservation.py
```

`tools/check_preservation.py`は教材本体78ファイルとJuliaコアを確認します。
保存時の80ファイルのハッシュを保持し、今回の対象読者の分離で変更した文書・生成スクリプトは`audience_revision_sha256.json`の改訂ハッシュで照合します。
数値計算のソース・分子データ・第一版は元のハッシュで照合します。公開用のルート`README.md`と`.gitignore`は対象外です。
単独配布のため元教材を置かない場合は、元教材の照合部分だけ実行できません。これは生成した計算の検証とは別の確認です。
`tools/build_notebook.py`は本文を再生成して出力を消すので、修正後は`run_notebook.py`で保存済み出力を更新してください。

他の構造を生成する例：

```sh
H4_GEOMETRY=4 bash tools/run_author.sh canonical boys_ordered boys
python tools/analyze_all.py --geometry 4
```

Juliaの初回JITコンパイル・参照MPO生成・1000スイープの最適化には数十分以上を要する場合があります。参加者用には生成済み結果を配布してください。
JuliaのSerializationファイルは環境依存なので、参加者向けにはNumPyゲート・JSON・QPYを使います。

## 数値と解釈

- 比較する基準値は、この固定構造・STO-3G模型の厳密解で、実験値ではありません。
- 元のDMRG 20スイープが未収束だった例は `logs/diagnostic_dmrg20_canonical.log` に記録しています。
- 本版では200スイープへ増やして参照状態を検証します。DMRGアルゴリズムや回路最適化コアは変えていません。
- `f`はMPSとの重なりの絶対値、状態忠実度はその2乗。`delta`は論文の演算子誤差指標です。
- 参照MPO自体にもSVD収縮の誤差があるので、密な100分割二次Trotter行列・厳密行列指数との比較を別々に記録しています。
- 同じDMRG参照状態と厳密時間発展でGaussian更新だけを行う対照計算も付け、回路近似の偏りと区別します。
- 量子ビット変換内部のPauli簡約の数値閾値により、FCIのエネルギーと量子ビットHamiltonianの固有値には約1e-8 Haの差が残る場合があり、両方をメタデータに記録しています。
- 推定履歴の分散は著者のGaussian更新のパラメータで、全誤差を含む信用・信頼区間ではありません。
- ショット版は理想測定の全0イベントの二項周辺分布から生成したもので、IBM実機結果ではありません。
- 最終反復の時間発展回数が大きくなった回路は、実機で実行可能な深さとは限りません。この発展演習は実機へ自動投入せず、資源の増加も比較します。
- 元のPythonコードは古いIBM認証をローカルシミュレーター経路でも初期化するため、その接続処理を分離し、同じ回路・同じ確率分布・同じ推定更新を使う実行アダプタを用意しています。

## 出典

Kanno et al., PNAS 122, e2425026122 (2025), [論文](https://doi.org/10.1073/pnas.2425026122)、補足S1/S2/S7。
著者コード固定コミット：`9cd88464c5c50dda783150661813609092d7f2b6`。
著者コードのLGPL-3.0ライセンスと出所を保持しています。

<!-- GENERATED_RESULTS -->
## 実行済みの比較（直鎖ID 0）

| 軌道 | 準備の重なり f | 状態忠実度 f² | 時間発展 δ | ショット版の最終推定 / Ha |
|---|---:|---:|---:|---:|
| canonical | 0.963732 | 0.928780 | 0.014266 | -0.271082 |
| boys | 0.993985 | 0.988007 | 0.009396 | 0.005619 |
| boys_ordered | 0.996375 | 0.992762 | 0.007726 | -0.039616 |

厳密ギャップは約0.015905 Haです。推定値のずれを隠さず、回路近似とGaussian推定の限界を議論します。
