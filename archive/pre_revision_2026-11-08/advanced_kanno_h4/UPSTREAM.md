# 著者実装の出所と変更範囲

- 論文：S. Kanno et al., *Tensor-based quantum phase difference estimation for large-scale demonstration*, PNAS 122, e2425026122 (2025), https://doi.org/10.1073/pnas.2425026122
- 公開コード：https://github.com/sk888ks/TQPDE_open
- 取得コミット：`9cd88464c5c50dda783150661813609092d7f2b6`
- 著者コードのライセンス：`upstream/LICENCE`（LGPL-3.0）。原文・著者表記を保持しています。

## そのまま実行する箇所

`runtime/Gate_generation_from_MPS_MPO.jl` は `upstream/src/Gate_generation_from_MPS_MPO.jl` のバイト同一コピーです。
`tools/check_preservation.py`で一致を検査します。変更しているのはデータ保存先を与える`runtime/input.jl`と実行ドライバです。

- `get_psi_sp`：DMRGで基底・第一励起MPSを作り、補助量子ビットを付けて和を取る。
- `get_gate_second_trotter`：二次Trotterの細分化とMPO収縮・SVDを使って参照MPOを作る。
- `get_and_save_prep_MPSBW_2d_list`：重ね合わせMPSに合わせた隣接2量子ビットbrick-wall回路を、環境テンソルとSVDによって最適化。
- `get_and_save_evol_BW_2d_list`：参照MPOに合わせた隣接2量子ビットbrick-wall回路の最適化。
- 保存された4×4ゲートは、著者のPythonコードと同じ層順・量子ビット位置でQiskitに読み込む。

## 分子前処理

`src/h4_models.py`は著者のhexatriene用スクリプトと同じ流れをH4へ適用したものです。
RHF → Boys局在化 → 交換積分によるDEAPの遺伝的アルゴリズム → ElectronicEnergy → InterleavedQubitMapper(JordanWignerMapper)。
H4/STO-3Gは4電子4空間軌道を全て使い、凍結軌道はありません。
正準／Boys／Boys+順序最適化の3表現で電子積分を作り直します。

構造は補足Fig. S6の2 Å直鎖・2 Å正方形と、その座標の線形補間を復元しました。
著者リポジトリにH4の入力XYZはありませんが、5構造の厳密ギャップは補足Table S1の表示桁と一致します。
これは構造の妥当性の強い照合であり、未公開の元入力とバイト同一という主張ではありません。

## 明記する設定差と実行アダプタ

- H4の設定は補足S7に合わせ`d_prep=d_evol=8`。最適化は1000スイープ、seed=1、Δt=0.1、参照MPOは100分割二次Trotter。
- 著者サンプルのDMRG 20スイープでは正準H4の基底状態が未収束だったため、200スイープへ増やした。アルゴリズムは変更せず、得た状態のエネルギーを厳密対角化と照合してから圧縮する。
- Juliaはローカルの1.11.1を使い、ITensors 0.6.16・PyCall 1.96.4・Infiltrator 1.7.0は著者と同じ版。論文READMEのJulia 1.10.2とは異なる。Manifestを保存。
- Python 3.12.7/Qiskit 2.4.1/PySCF 2.13.1で周辺処理を実行。著者のPython 3.11.8/Qiskit 1.2.4環境そのものではない。
- PythonのGaussian関数は著者ファイルからASTで直接読み込む。時間選択、探索点、curve_fit、Gaussianの積による更新は著者コードと同じ式。
- 古いIBM認証・Session初期化をローカル推定から分離した。全0振幅を線形代数で評価し、ショット実験はその二項周辺分布から生成する。理想的な測定分布は同じだが、Aerの乱数列を再現するものではない。
- Fire Opal/IBM実機への接続・誤差抑制は実行していない。

したがって本演習は「著者のテンソル回路最適化を実行するH4発展演習」であり、論文のハードウェア実験・全ソフトウェア環境・表S1の推定値まで完全同一に再現したという意味ではありません。
