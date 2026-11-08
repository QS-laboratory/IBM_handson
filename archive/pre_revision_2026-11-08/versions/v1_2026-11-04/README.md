# MPS・MPO・量子回路の入門デモ

2026年11月9日開催予定のハンズオンに向けた、60–90分の導入教材です。
`01_mps_mpo_circuits.ipynb` に日本語の説明、演習、講師向け回答の要点、回路図、実行結果を収録しています。
補助実装は `demo_core.py` にあります。2つのファイルを同じフォルダに置いてください。

## 起動

動作確認環境はPython 3.12.7です。新しい環境では次のように準備します。

```sh
cd /Users/sato/project/FTQC_HandsOn/hands_on
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m ipykernel install --prefix .venv --name ftqc-hands-on --display-name 'FTQC Hands-on'
jupyter lab 01_mps_mpo_circuits.ipynb
```

Jupyterでは「FTQC Hands-on」のカーネルを選択し、上から実行します。
依存バージョンは既存環境で確認したものです。新しいvenvへのインストールは未検証です。
実行済みの出力も保存しているので、最初はパラメータを変えずに図を読み進められます。

## 構成と所要時間の目安

| 内容 | 目安 |
|---|---:|
| GHZ状態、明示的なMPS、テンソル図 | 15分 |
| SVDとボンド次元による近似 | 10分 |
| MPSで与えた基底状態に状態準備回路を合わせる | 10分 |
| ZZ時間発展の行列・MPO・回路、MPOをMPSに作用させる | 15分 |
| 非可換Ising模型、MPO圧縮、浅いユニタリ回路の最適化 | 20分 |
| 短時間回路の反復と位相差推定への接続 | 15分 |

参加者はθ、時間幅、ボンド次元、回路の層数、ショット数を変更できます。
4時間全体のうち、この教材を基礎パートとし、著者コードやIBM実機での演習を別パートにする構成を想定しています。

## 実装と検証範囲

全処理はローカルの理想計算です。IBM実機への送信、Fire Opalの呼び出し、アカウント情報の読取りはありません。
PythonのNumPy/SciPy/Qiskitを使用しています。

- 全13コードセルを同一Pythonプロセスで順に実行し、出力と12枚の図を保存しました。
- Jupyterカーネルを起動するnbclientの方式は、この環境のローカルソケット制限で利用できませんでした。Jupyter UIでの通し実行は未確認です。
- MPS/MPOからの復元、回路のビット順、ZZ回路、MPOのMPSへの作用、位相差推定の確率式を数値検証しています。
- `python demo_core.py` でも復元と確率式の基本検証を実行できます。

密な配列による教材用の実装です。一般の状態や演算子を配列に展開するため、大規模系には使えません。
MPS/MPOのボンド切詰めと、ユニタリ回路の最適化を区別して説明しています。
最終QPDEの重ね合わせ準備には、小規模系用の汎用StatePreparationを使用します。
論文のHubbard模型、テンソルネットワークによる大規模最適化、Fire Opalを含む実機実験の再現は対象外です。

## 参照

- Kanno et al., PNAS 122, e2425026122 (2025): https://doi.org/10.1073/pnas.2425026122
- 補足資料 S1・S2：親フォルダ内のPDF。
- 著者の公開コード: https://github.com/sk888ks/TQPDE_open
- Qiskit Quantum Information: https://quantum.cloud.ibm.com/docs/en/api/qiskit/quantum_info

教材コードは説明用に独立実装したもので、著者のコードを転載していません。
