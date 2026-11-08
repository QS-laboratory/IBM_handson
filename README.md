# IBM Quantum FTQC Hands-on workshop

2026年11月9日（月）14:00–18:00、東京大学 本郷キャンパス。

小さな行列からH₂、H₄へ進み、状態準備・時間発展・測定と、エネルギー差の推定を学びます。

## 演習を開く

各冊はデータと補助関数を内蔵し、前の冊の実行状態に依存しません。Colabで自分のコピーを保存し、上から実行してください。

| 順番 | 内容 | Colab |
|---|---|---|
| 01 | [2×2行列：QPE・BPE・BPDE](hands_on/exercises/01_matrix_phase_estimation.ipynb) | [開く](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/exercises/01_matrix_phase_estimation.ipynb) |
| 02 | [H₂・1 qubit：電子状態と化学への入口](hands_on/exercises/02_h2_one_qubit.ipynb) | [開く](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/exercises/02_h2_one_qubit.ipynb) |
| 03 | [H₂・JW・4 qubit：厳密準備とTrotter分解](hands_on/exercises/03_h2_jw_bpde.ipynb) | [開く](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/exercises/03_h2_jw_bpde.ipynb) |
| 04 | [H₄：回路近似の評価・実機接続](hands_on/exercises/04_h4_circuit_assessment.ipynb) | [開く](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/exercises/04_h4_circuit_assessment.ipynb) |
| 05 | [H₄：反復的Bayesian推定](hands_on/exercises/05_h4_iterative_bayes.ipynb) | [開く](https://colab.research.google.com/github/QS-laboratory/IBM_handson/blob/main/hands_on/exercises/05_h4_iterative_bayes.ipynb) |

[実行方法・記号・測定モード](hands_on/exercises/README.md)を参照してください。全ての問いに解答例があります。
演習4・5は、厳密確率・有限shotsの模擬測定・IBM実機を切り替えられます。実機送信は初期状態では無効です。
Colab UIとIBM実機での動作検証は、ローカル実行の検証と区別しています。

## 事前準備

- [① IBM Quantumの利用登録](docs/01_ibm_quantum_registration.md) / [PDF](output/pdf/01_ibm_quantum_registration.pdf)
- [② Google Colabからの接続確認](docs/02_colab_connection.md) / [PDF](output/pdf/02_colab_connection.pdf)

本編の実行にPySCF・Julia・Aerは不要です。APIキーをコード・出力・GitHubへ記載しないでください。

## 参考研究

S. Kanno et al., *Tensor-based quantum phase difference estimation for large-scale demonstration*, PNAS 122, e2425026122 (2025).
[論文](https://doi.org/10.1073/pnas.2425026122) / [著者コード](https://github.com/sk888ks/TQPDE_open)

本教材は論文の分子・実機実験全体の再現ではありません。

<details><summary>旧版・開発資料</summary>

[旧教材の保存場所](archive/README.md)。新本編の元教材・再生成用データは `hands_on/self_study/` に保持しています。

</details>
