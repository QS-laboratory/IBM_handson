# 分子を題材としたFTQCハンズオン — 新版

参加者向け入口は [リポジトリのREADME](../README.md) です。
第一版は `../versions/v1_2026-11-04/` に保存し、本編5教材へ再構成しました。

## ローカルで実行する

検証環境はPython 3.12.7、Qiskit 2.4.1、qiskit-ibm-runtime 0.43.1、NumPy 2.4.6、SciPy 1.18.0、Matplotlib 3.10.0、PySCF 2.13.1です。
新規venvへのインストールとColab標準環境での検証は未実施です。

```sh
cd IBM_handson/hands_on
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-instructor.txt
python -m ipykernel install --prefix .venv --name ftqc-hands-on --display-name 'FTQC Hands-on'
jupyter lab
```

Jupyterで「FTQC Hands-on」のカーネルを選び、各ノートを上から実行します。
配布された実行済み出力を読んで進めることもできます。
全ノートのIBM接続・送信フラグはFalseで、通常の通し実行ではQPUを使いません。

## 構成

- `overview/00_overview.ipynb`：材料・分子・電子状態からFTQCへの導入。講師用の時間配分と問いを含みます。
- `00_preparation.ipynb`：Colab起動、基本回路、IBM認証、Bell回路の試行。
- `01_circuits_and_ibm.ipynb`：分子の準備済み回路比較、測定計画、ISA確認、実機送信。
- `02_inference.ipynb`：実機結果取得または模擬測定の読込み、ベイズ更新、誤差・モデル診断。
- `03_mps_mpo_advanced.ipynb`：第一版の数学的説明と、分子のMPS/MPOへの橋渡し。
- `src/workshop.py`：教材共通API。インポート時には接続・送信しません。
- `src/molecules.py`：分子の生成・積分・活性空間・厳密参照計算。
- `src/circuit_building.py`：準備回路と短時間回路の事前生成。
- `src/tensor_demo.py`：第一版の小規模テンソル例。
- `data/reference/`：分子条件JSONと数値NPZ。
- `data/circuits/`：QPY回路、最適化条件・指標JSON、準備回路パラメータ。
- `data/example_results/`：出所を明記した模擬測定JSON。IBM実機データはまだありません。

## 生成・検証

リポジトリのルートで実行します。

```sh
python hands_on/tools/build_molecules.py
python hands_on/tools/build_circuits.py
python hands_on/tools/build_example_results.py
python hands_on/tools/build_notebooks.py
python -m pytest hands_on/tests -q
python hands_on/tools/execute_notebooks.py
python hands_on/tools/build_manifest.py
```

`build_notebooks.py`は説明を含めてノートを再生成し、保存した出力を消します。ノート本文を変更する場合はこの生成元も更新してください。
`execute_notebooks.py`はJupyterカーネルを使わず、各ノートを別の名前空間でPythonとして順次実行し、図と出力を保存します。
IPythonマジックは使っていません。これはJupyter/ColabのUIでの通し実行を代替して確認したものではありません。
ローカル図の確認用PNGと実行記録は `output/qa/`、実行時の解析結果は `outputs/` に保存し、Git管理から除外します。

## 分子・回路の条件

幾何は固定した教材用の平面構造。基底はSTO-3G、π空間のCASCIです。σ軌道等の占有軌道は二重占有で凍結しています。
PySCFのRHF分子軌道のうち、炭素2pz空間への射影が大きい軌道を選びます。
基底一重項と最低三重項を粒子数と全スピンで分離し、回路にはともにM_S=0の状態を使います。
エネルギーは独立したPySCFの一重項・M_S=1の三重項FCI計算と照合しています。
配列は新版全体でQiskit標準のlittle endian。補助量子ビットは最大番号です。

本編の短い回路は、実数RY/CX状態フィットとPauli積回路の時間係数の最適化による教育用実装です。
論文のMPS/MPOスイープではありません。短時間回路の最適化はエチレンでは全演算子、ブタジエンでは固定ランダムプローブ12本に対して行い、最終評価は全演算子で行います。
厳密行列・全状態・QPYの読み込みと最適化は小規模系向けで、指数的な資源を使います。

推定区間は、有限範囲の事前分布と簡略化した二項尤度に条件付けた95%ベイズ信用区間です。
分子模型・基底・活性空間・回路近似の全誤差を含む保証ではありません。
推定のラウンドは事前に固定した測定セットで、適応的な実機再投入ループはまだ実装していません。

## GitHub/Colab

配置先は `QS-laboratory/IBM_handson`、既定ブランチは`main`を想定しています。
GitHub側にはまだ教材をアップロードしていません。リポジトリの公開状態・配置後のリンクは講師が確認してください。
このColab導線はGitHubリポジトリを未認証でcloneできることを前提とします。非公開運用には別の取得手順が必要です。
各ノート冒頭で同リポジトリのコードとデータを取得します。配布時には `GITHUB_REF` とColabリンクを検証済みタグまたはコミットへ固定してください。
第一版・論文PDFはColab実行には不要です。PDFは独立した元資料として扱い、GitHubへ自動配布する対象に含めていません。
