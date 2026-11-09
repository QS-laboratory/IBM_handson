"""Build static Markdown inventories from the frozen R=2 reference metadata."""
from pathlib import Path
import json


def circuit_inventory():
    data=Path(__file__).resolve().parents[1]/'data'
    cat=json.loads((data/'h4_pec_v1/R_2p00000000/catalog.json').read_text())
    reference=json.loads((data/'h4_inventory_reference.json').read_text())
    assert reference['reference_sha256']==cat['evolution']['16']['report']['reference_sha256']
    lines=[r'''## 2.1．保存回路の一覧：χ・M・精度・CXを対応させる

以下の数値は**固定データ集のR=2 Å**の例です。選択構造を変えたときの値は後続セルで計算します。
χはMPSのボンド次元、Mは交互の隣接2-qubitゲート層数で、同じものではありません。

| 準備の保存範囲 | χ | 保存されたM | 状態・構造 | この演習で選択 |
|---|---|---|---|---|
| 固定データ集v1 | 6 | 4, 6, 8 | S₀・T₀、全20構造 | 可 |
| 厳密準備 | 使用しない | exact | S₀・T₀、全20構造 | 可 |
| 研究アーカイブの別系列 | 4 | 1, 2, 4, 8 | R=2 ÅのS₀のみ | v1未収録 |
| χ=5、χ=4/M=6 | — | — | 該当する保存回路なし | 不可 |

χ=4/5/6 × M=4/6/8の全組合せが揃っているわけではありません。
この演習の20構造走査では、S₀とT₀が揃ったχ=6の保存系列を使います。最適化は再実行しません。

「回路化後fidelity」はFCI状態との重なりです。MPS目標への重なりとは区別します。
下表の「実機変換前CX」は、保存済みのrz/sx/x/cx回路の値で、IBM backendへの配置・再トランスパイル前です。
**U(4)ブロックを基本ゲートへ分解する前の回路には、明示的CXはありません。費用0という意味ではなく未分解です。**
厳密準備は初めからRy/CXで合成しています。

| 状態 | χ | M | MPS化後F_state | 回路化後F_state | 未分解2-qubitブロック数 | 実機変換前CX(Pk) | 実機変換前CX(A) |
|---|---|---|---|---|---|---|---|''']
    for m in ['4','6','8','exact']:
        p=cat['preparation'][m]
        for state,label in [('S0','S₀'),('T','T₀')]:
            chi='—' if m=='exact' else '6'
            mps='—（経由せず）' if m=='exact' else f"{p['fidelities'][state+'_MPS']:.10f}"
            blocks='Ry/CX合成' if m=='exact' else str(p['reports'][state]['two_qubit_blocks'])
            lines.append(f"| {label} | {chi} | {m} | {mps} | {p['fidelities'][state]:.10f} | {blocks} | {p['resources'][state]['cx']} | {p['resources']['A']['cx']} |")
    lines.append(r'''
AはS₀/T₀両枝を一緒に準備する9 qubit回路です。各行のCX(A)は同じ両枝回路の値を再掲しており、S₀とT₀で足しません。
全BPDEにはA†と時間発展も必要です。MPS化後がほぼ1でも、回路化後の一致はMに依存します。

### 保存時間発展の一覧：MPOとゲート列を分ける

演算子は**MPO化**します。保存系列は発展M=8・16・32と厳密参照です。
三つのMは同じ参照MPO（二次Trotter100分割、SVD cutoff=10⁻¹⁸）を別々の長さへ回路化したものです。
MPOは非ユニタリになり得るため、表のMPO fidelityは正規化Hilbert–Schmidt重なり
$$F_{\rm MPO}=\frac{|\mathrm{Tr}(U_{\rm exact}^\dagger V_{\rm MPO})|^2}
{\mathrm{Tr}(U_{\rm exact}^\dagger U_{\rm exact})\mathrm{Tr}(V_{\rm MPO}^\dagger V_{\rm MPO})}$$
を使います。回路はユニタリなので、この定義は前述の$|\mathrm{Tr}(U^\dagger U_{app})/256|^2$と一致します。
MPO fidelityの高さだけではユニタリ性は保証しません。参照MPO行列自体を実機へ投入しない理由です。
表の回路化後fidelityも、MPO目標ではなく共通のHから得た厳密Uと比較しています。

| 発展M | MPO化後F_op | 回路化後F_op | 未分解U(4)数（順変換を除く） | 実機変換前CX（変換なし） | 実機変換前CX（符号付き順変換込み） |
|---|---|---|---|---|---|''')
    for m in ['8','16','32']:
        e=cat['evolution'][m]
        lines.append(f"| {m} | {reference['mpo_operator_fidelity_normalized_HS']:.16f} | {e['operator_fidelity']:.12f} | {e['report']['two_qubit_blocks']} | {e['report']['CX_evolution_only']} | {e['resources']['cx']} |")
    lines.append(r'''| exact | —（MPOを経由せず） | 1（数値精度） | 256×256密行列unitary | 未分解・未評価 | 未分解・未評価 |

実際のBPDEは右端の順変換込みの回路を使います。未分解の発展U(4)列も明示的CX数は0ですが、基本ゲートの費用0を意味しません。
厳密発展はexp(−iHΔt)を密行列として持つシミュレーション参照で、短いネイティブ回路ではありません。
有限予算の回路フィットを保存しており、Mを増やしたものが完全収束したとは主張しません。
MPO fidelityの算出条件・参照行列のハッシュは `data/h4_inventory_reference.json` に記録しています。
''')
    legacy=['### 参考：研究側に保存されたχ=4のS₀（R=2 Å）',
            'この別系列の参照FCIに対する保存評価値です。T₀・全20構造が揃っていないため、下の走査では選択しません。',
            '| χ | M | MPS化後F_state | 回路化後F_state | 未分解U(4)数 | 実機変換前CX |',
            '|---|---|---|---|---|---|']
    for row in reference['archived_chi4_reference']['rows']:
        legacy.append(f"| 4 | {row['M']} | {row['mps_fidelity']:.10f} | {row['circuit_fidelity']:.10f} | {7*row['M']//2} | {row['native_cx']} |")
    return '\n'.join(lines).replace('### 保存時間発展の一覧', '\n'.join(legacy)+'\n\n### 保存時間発展の一覧')
