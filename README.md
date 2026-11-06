# IBM Quantum FTQC Hands-on workshop

2026年11月9日（月）14:00–18:00、東京大学 本郷キャンパス。
材料科学・量子化学を入口に、分子の状態準備・時間発展・測定からエネルギー差を推定する流れを体験します。

**主題材は1,3-ブタジエン、実機演習の初期設定は軽量なエチレン**です。
どちらも教材用固定構造、STO-3G、π活性空間のCASCIを用い、最低一重項–最低三重項の垂直エネルギー差を対象にします。
小規模模型の厳密解を基準にしており、実験の吸収エネルギーを予測する教材ではありません。

## 教材を開く

GitHubへ配置後、以下のボタンからGoogle Colabで開きます。現時点のローカル成果物だけではColabリンクはまだ利用できません。
各ノートで「ファイル → ドライブにコピーを保存」し、標準CPUランタイムで実行してください。

| 時間 | 教材 | Colab |
|---|---|---|
| 14:05–15:00 | [概観：材料・分子・電子状態と量子計算](hands_on/overview/00_overview.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/overview/00_overview.ipynb) |
| 15:10–16:00 | [準備：基本操作・IBM接続](hands_on/00_preparation.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/00_preparation.ipynb) |
| 16:10–17:00 | [演習1：回路比較・実機投入](hands_on/01_circuits_and_ibm.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/01_circuits_and_ibm.ipynb) |
| 17:10–17:55 | [演習2：推定更新・誤差評価](hands_on/02_inference.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/02_inference.ipynb) |
| 任意・復習 | [発展：MPS/MPO・SVD・回路生成](hands_on/03_mps_mpo_advanced.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/03_mps_mpo_advanced.ipynb) |

開会14:00–14:05、休憩15:00–15:10・16:00–16:10・17:00–17:10、まとめ17:55–18:00。懇談会18:15–20:00。

## 事前準備

- GoogleアカウントでColabへログインできること、会社PC・社内規則でGitHub/Colabの利用が可能なことを確認してください。
- IBM実機へのアクセスは各自のAPIキーを使用します。今回利用可能な所属企業の枠、または講師が案内するFlex枠の権限が必要です。
- APIキーはColabのシークレット `IBM_QUANTUM_API_KEY` に登録します。教材セル・出力・GitHubへ記載しません。
- シミュレーション・配布データの解析は、IBMアカウントなしでも実行できます。

全員がColabで操作し、実機送信はチーム代表者が行います。実機送信セルは初期状態では無効です。
実機待ち時間に備え、**模擬測定データ（実機結果ではない）**を配布しています。

## 開発・検証状況

分子データ・準備済み回路・模擬測定データを生成し、5教材をローカルで順次実行しています。
IBM実機・Colab実環境での動作と授業全体の所要時間は未検証です。実機性能・利用枠の所要時間は事前ベンチマークで確定します。

詳しい起動・再生成手順は [開発用README](hands_on/README.md)、当日の運営と未完了の確認項目は [講師用ガイド](hands_on/INSTRUCTOR.md) を参照してください。
第一版は [versions/v1_2026-11-04](versions/v1_2026-11-04/) に独立保存しています。

## 参考研究

S. Kanno et al., *Tensor-based quantum phase difference estimation for large-scale demonstration*, PNAS 122, e2425026122 (2025).
[論文](https://doi.org/10.1073/pnas.2425026122) / [著者コード](https://github.com/sk888ks/TQPDE_open)

本教材は独立した説明用実装です。論文のテンソルネットワークによる最適化やFire Opalを含む実機実験の再現ではありません。
