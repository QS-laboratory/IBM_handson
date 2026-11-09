# IBM Quantum FTQC Hands-on workshop

2026年11月9日（月）14:00–18:00、東京大学 本郷キャンパス。

小さな行列からH₂、H₄へ進み、状態準備・時間発展・測定と、エネルギー差の推定を学びます。

## 演習を開く

各冊は単独で実行できます。Colabで自分のコピーを保存し、冒頭のインストールセルを有効にして上から実行してください。

| 順番 | 教材 | 内容 | Google Colab |
|---|---|---|---|
| 01 | [行列の位相推定](hands_on/exercises/01_matrix_phase_estimation.ipynb) | QPE・BPE・BPDE、Gaussian更新 | [開く](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/exercises/01_matrix_phase_estimation.ipynb) |
| 02 | [H₂・1 qubit](hands_on/exercises/02_h2_one_qubit.ipynb) | 電子配置とスピン、S₀–S₂の位相差 | [開く](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/exercises/02_h2_one_qubit.ipynb) |
| 03 | [H₂・Jordan–Wigner・4 qubit](hands_on/exercises/03_h2_jw_bpde.ipynb) | 02と同じ流れでS₀–T₀の差を測定 | [開く](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/exercises/03_h2_jw_bpde.ipynb) |
| 04 | [H₄のポテンシャル曲線](hands_on/exercises/04_h4_circuit_assessment.ipynb) | 保存回路の表、期待値のPEC、Lを増やす位相実験 | [開く](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/exercises/04_h4_circuit_assessment.ipynb) |
| 05 | [H₄のGaussian更新](hands_on/exercises/05_h4_iterative_bayes.ipynb) | シミュレーションで測定・フィット・事後更新を反復 | [開く](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/exercises/05_h4_iterative_bayes.ipynb) |
| 06 | [H₄の実機デモ](hands_on/exercises/06_h4_hardware_demo.ipynb) | 同じ9点で理想確率・有限shots・IBM実機を比較 | [開く](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/exercises/06_h4_hardware_demo.ipynb) |
| 07 | [H₄の実機信号と回路長](hands_on/exercises/07_h4_hardware_depth_scan.ipynb) | L=1・2・4の位相走査で信号の減衰を調べる | [開く](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/exercises/07_h4_hardware_depth_scan.ipynb) |
| 08 | [H₄のjoint-state圧縮と誤り抑制](hands_on/exercises/08_h4_joint_compression_suppression.ipynb) | joint準備・較正配置・twirling・DDをfakeと実機で比較 | [開く](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/exercises/08_h4_joint_compression_suppression.ipynb) |
| 09 | [H₄のjoint準備だけの実機デモ](hands_on/exercises/09_h4_joint_preparation_demo.ipynb) | 06をjoint準備・L=1・位相一周期へ最小変更 | [開く](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/exercises/09_h4_joint_preparation_demo.ipynb) |

01〜05はシミュレーション、06〜09はシミュレータと実機の比較です。06も初期設定ではIBM接続・ジョブ送信を行いません。
[実行方法・記号・設定](hands_on/exercises/README.md)を参照してください。

## 事前準備と資料

- [① IBM Quantumの利用登録](docs/01_ibm_quantum_registration.md)
- [② Google Colabからの接続確認](docs/02_colab_connection.md)
- [Kanno論文・補足資料・参考資料](docs/README.md)

シミュレーションにはIBMアカウントは不要です。06〜09で実機を使う方は本人のAPIキーと利用可能なインスタンスCRNを用意してください。
APIキーは非表示入力し、コード・出力・GitHubへ記載しないでください。
PySCF・Juliaは不要です。08・09と06のfakeモードではQiskit Aerも使用します。

## 参考研究

S. Kanno et al., *Tensor-based quantum phase difference estimation for large-scale demonstration*, PNAS 122, e2425026122 (2025).
[論文](https://doi.org/10.1073/pnas.2425026122) / [著者コード](https://github.com/sk888ks/TQPDE_open)

本教材は論文の分子・実機実験全体の再現ではありません。

<details><summary>検証記録・旧版・開発資料</summary>

- [正式版01〜06のローカル実行記録](docs/validation_20261109_promotion.json)
- [07のローカル実行・模擬ジョブ検証](docs/validation_20261109_depth_scan.json)
- [08のjoint圧縮・fake実行・模擬実機ジョブ検証](docs/validation_20261109_joint_compression.json)
- [09のjoint準備のみ・fake実行検証](docs/validation_20261109_joint_only.json)
- [差し替え前の01〜05・生成コード・テスト](archive/pre_promotion_2026-11-09/README.md)
- [それ以前の教材と自習資料](archive/README.md)

旧版のColab実行記録は、[1203538](docs/colab_validation_1203538.json)・[900224b](docs/colab_validation_900224b.json)の各コミットに対応します。

</details>
