# 直鎖H4のFull-CIポテンシャルエネルギー曲線

[PEC図（PNG）](h4_linear_pec.png) / [ベクトル図（SVG）](h4_linear_pec.svg) / [全69点のCSV](linear_stretch.csv)

- 対象：中性H4、最低一重項S0と最低三重項T1。
- 配置：Hを一直線に等間隔で置き、隣接距離Rを0.60〜4.00 Åまで0.05 Å刻みで変える。
- 座標：`(0,1.5R,0), (0,0.5R,0), (0,-0.5R,0), (0,-1.5R,0)` Å。
- 基底：STO-3G。全4電子・全4空間軌道のFull-CIで、凍結軌道・活性空間の切詰めはない。
- 縦軸：電子エネルギーと核間反発エネルギーを合計した全エネルギー。
- 同じ等間隔直鎖という拘束の下でのPEC。曲線の最小点は、全原子座標を最適化した平衡構造を意味しない。

軌道基底には対称直交化AOを用いる。全軌道Full-CIなので、正準軌道や局在化軌道へユニタリ変換しても厳密なエネルギーは変わらない。
SCF解を追跡する必要がなく、SCFの未収束や解の乗り換えによる曲線の不連続を避けられる。

全M_S=0の36行列式空間からS²を対角化し、S=0の20次元空間、S=1の15次元空間をそれぞれ厳密対角化する。
全69点でスピン、固有状態残差、HamiltonianとS²の可換性を確認し、標準PySCF direct_spin0/direct_spin1(M_S=1)のエネルギーとも1e-9 Ha以内で照合した。
長距離側の縮退でも、S0とT1をスピンで分離する。

R=2.00 Åの結果は既存の発展演習の参照値と一致：

| 状態 | 全エネルギー / Ha |
|---|---:|
| S0 | -1.897780645990 |
| T1 | -1.881875688074 |

差は0.015904957916 Ha = 0.432795953084 eV。
破線は同じSTO-3G基底で計算した4個の孤立H原子のエネルギー和。

## 再生成

```sh
OPENBLAS_NUM_THREADS=2 python compute_pec.py
MPLBACKEND=Agg python plot_pec.py
```

PySCF 2.13.1、NumPy、SciPy、Matplotlibを使用。計算条件と数値はmetadata.json、核間反発・スピン・残差はCSVにも保存。
元の教材・発展演習ノート・既存データは変更せず、このpecフォルダへ追加した。

参考： [PySCF FCI](https://pyscf.org/user/ci.html)、[一重項FCIソルバー](https://pyscf.org/_modules/pyscf/fci/direct_spin0.html)。

## DMRGとの比較

[Full-CIとDMRGの重ね描き（PNG）](h4_linear_fci_dmrg_pec.png) / [SVG](h4_linear_fci_dmrg_pec.svg) / [DMRG結果CSV](linear_dmrg.csv)

`pyqmd2/src/pyqmd2/experimental/dmrg_benchmark.py` のblock2ベンチマークを参考に、PEC専用の高速な全電子DMRGを追加した。
参照元のSHA-256はdmrg_metadata.jsonに記録。pyqmd2のコードや本編・論文再現のコードは変更していない。

- block2 0.5.3の電子SU(2)・two-site DMRG。DOCI／seniority-zeroのペア制限は使わない。
- 電子数4、全スピンS=0／S=1をそれぞれ固定。block2のspin引数は2Sなので0／2。
- H4は空間軌道4、局所次元4で、保守的な全空間Schmidt上限は16。ボンド次元を16で打ち止めにする。
- 独立2シード20261024／20261025を各点・各スピンで実行。距離間でMPSを使い回さない。
- 初めの4スイープに1e-5,1e-5,1e-6,1e-6のnoiseを与え、その後はnoise=0。
- 最後の2スイープのエネルギー変化1e-11 Ha未満、discarded weight 1e-12未満で停止。上限は16スイープ。
- DMRGにはFCI行列、FCIベクトル、FCI初期値、FCIエネルギーを与えない。FCI比較は計算後の検証だけ。
- 核間反発を含む全エネルギー、STO-3G、Full-CIと同じ対称直交化AO積分を使う。

今回のローカル実行では全69点×2状態×2シード＝276回が約4.43秒で終了した（ライブラリのimport時間を除く実行区間）。
Full-CIとの差の最大絶対値はS0で1.07e-14 Ha、T1で1.24e-14 Ha、独立シード間の最大差は1.69e-14 Ha。
この小系ではボンド次元上限に到達しており、差は浮動小数点丸め程度。図のDMRG点はFCI曲線に重なる。
全スピンの期待値と全276回の収束履歴はdmrg_metadata.jsonに保存した。

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python compute_dmrg.py
MPLBACKEND=Agg python plot_dmrg.py
```

必要な追加依存はblock2（pyblock2を含む）。今回確認した版は0.5.3。

前回長時間を要したのは、主として参照MPOから浅い量子回路を作る1000スイープのSVD最適化であり、DMRGエネルギー計算そのものではない。
PEC計算ではその回路最適化を実行しない。
