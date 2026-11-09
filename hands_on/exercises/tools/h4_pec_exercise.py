"""Small PEC scan section appended to exercise 04 before the handoff."""
from pathlib import Path
import nbformat as nb


def pec_exercise_cells():
    md=nb.v4.new_markdown_cell;code=nb.v4.new_code_cell
    return [md(r'''## 6.3．全20構造を走査する：BPE・BPDEのピークと〈H〉

`PEC_PREP_M`（exact/4/6/8）、`PEC_EVOLUTION_M`（exact/8/16/32）、`PEC_L`（1〜1024の整数）を選びます。
これは**常にshot noiseなし**のシミュレーションです。MODE=ibmでも、この節から実機へ送信しません。
全20構造の固定データを使い、不足ファイルは同じキャッシュへ自動取得します。既定は準備M8・発展M32・L64で、主実験とは独立に設定します。
変更してこの節の設定・計算セル以降を再実行すれば、厳密準備／厳密発展との違いを比較できます。
比較前の構造・M・主実験のmeterと保存ファイルは保持し、この走査は別ファイルへ保存します。

各構造で実際の準備状態と保存発展から$a_k=\langle\widetilde\psi_k|U^L|\widetilde\psi_k\rangle$を計算し、
補助の位相βを9点で掃引します。BPEとBPDE Cの確率は
$$P_{\rm BPE}(\beta)=\tfrac12[1+\mathrm{Re}(e^{i\beta}a_0)],\qquad
P_{\rm BPDE,C}(\beta)=\tfrac14|a_0+e^{i\beta}a_T|^2.$$
どちらも定数＋cosβ＋sinβへ線形フィットでき、ピークを読み取れます。
枝を選んだピークβ₀、β_gapから$E_{0,el}=\beta_0/(L\Delta t)$、$\Delta E=\beta_{gap}/(L\Delta t)$、
$E_{T,total}=E_{0,el}+\Delta E+E_{nuc}$を得ます。T₀はBPE(S₀)＋BPDEの差から再構成します。
〈H〉は同じ準備状態で評価します。厳密準備なら〈H〉はFCIと一致しますが、近似発展の推定には偏りが残り得ます。

### 位相の枝：長時間のピークには追加情報が必要

βのピークだけではエネルギーは$2\pi/t$周期でしか分かりません。
この演習では、t=0から細かい補助時間で生存振幅の位相を追跡して枝を選びます。**FCIや〈H〉を使って正解の枝へ寄せません。**
$U=V\mathrm{diag}(e^{i\theta_n})V^\dagger$から、補助時間では$e^{i\theta_n\ell}$を使います。
近似回路の分数反復は保存Uの主値logによる連続補間で、実機でM回路を分数回実行したという意味ではありません。
これはシミュレーション上の追加情報です。この枝情報を、L固定の9測定だけから得たと解釈してはいけません。
各枝の重みは $w_n^{(k)}=|\langle v_n|\widetilde\psi_k\rangle|^2$ で、補助時刻の振幅は $a_k(\ell)=\sum_nw_n^{(k)}e^{i\theta_n\ell}$ です。
細かい刻みと倍の刻みの一致に加え、各枝の $D_k=\sum_nw_n^{(k)}|\theta_n|$ を使った導関数境界で未解像の巻き付きがないことを確認します。
検証できない場合は停止します。この節はピークの診断で、05のGaussian更新そのものではありません。
'''),code(Path(__file__).with_name('h4_pec_scan.py').read_text()),code('''RUN_PEC_SCAN = True  # Falseにすると全20構造の走査をスキップ
PEC_PREP_M = 8       # 'exact', 4, 6, 8（χ=6の保存回路）
PEC_EVOLUTION_M = 32 # 'exact', 8, 16, 32
PEC_L = 64          # 1..1024。総時間t=PEC_L*DT
PEC_RESULTS=[]
if RUN_PEC_SCAN:
    PEC_RESULTS=run_pec_scan(PEC_PREP_M,PEC_EVOLUTION_M,PEC_L)
    radii=[r['R'] for r in PEC_RESULTS]
    fig,axes=plt.subplots(2,3,figsize=(14,7))
    for j,(key,title) in enumerate([('S0','S0 total energy'),('T','T0 total energy'),('gap','T0 - S0 gap')]):
        peak=key+'_peak'+('_total' if key!='gap' else '')
        fci=key+'_FCI'+('_total' if key!='gap' else '')
        mean=key+'_expect'+('_total' if key!='gap' else '')
        for field,label,style in [(fci,'FCI','k-'),(mean,'<H>','--'),(peak,'BPE / BPDE peaks','o-')]:
            axes[0,j].plot(radii,[r[field] for r in PEC_RESULTS],style,label=label,ms=3)
        axes[1,j].axhline(0,color='black',lw=.7)
        for field,label in [(mean,'<H> - FCI'),(peak,'Peak - FCI')]:
            axes[1,j].plot(radii,[1000*(r[field]-r[fci]) for r in PEC_RESULTS],'o-',ms=3,label=label)
        axes[0,j].set(title=title,ylabel='Energy / Ha')
        axes[1,j].set(ylabel='Signed error / mHa')
        for ax in axes[:,j]:ax.set_xlabel('Adjacent H-H distance / Angstrom');ax.legend(fontsize=8)
    fig.suptitle(f'Preparation={PEC_PREP_M}, evolution={PEC_EVOLUTION_M}, L={PEC_L}; no shots')
    plt.tight_layout();plt.show();plt.close(fig)
    print('Quantity    peak RMSE / mHa    <H> RMSE / mHa')
    for key in ['S0','T','gap']:
        suffix='_total' if key!='gap' else ''
        error=lambda kind:np.array([r[key+'_'+kind+suffix]-r[key+'_FCI'+suffix] for r in PEC_RESULTS])
        print(key,1000*np.sqrt(np.mean(error('peak')**2)),1000*np.sqrt(np.mean(error('expect')**2)))
    print('Phase certificates:',sum(r['phase_certificate']['certified'] for r in PEC_RESULTS),'/20')
    print('Maximum cosine-fit residual:',max(r['max_fit_residual'] for r in PEC_RESULTS))
    print('Saved: outputs/exercise4_pec_scan.json')
'''),md('''### 課題6補足：ポテンシャル曲線と長時間の推定

1. 両方exactなら、FCI・〈H〉・ピークの曲線はどうなるか。
2. 近似準備のままLを増やすと、ピークは必ず〈H〉より良くなるか。発展を近似にするとどうか。
3. L=64の9点の位相掃引だけで、図の絶対エネルギーの枝が決まったといえるか。

<details><summary>解答例</summary>

1. 数値精度で一致します。E_nucは両状態へ加え、ギャップでは相殺します。
2. 必ずではありません。短時間の位相は〈H〉に近く、優勢な固有成分があれば長時間でFCIへ近づく場合がありますが、各Lで単調には改善しません。
   近似発展では別のユニタリの位相を追うため、準備を改善しても発展由来の偏りが残り得ます。図の符号付き誤差とRMSEを別々に読みます。
   準備exact＋発展近似では〈H〉がFCIと一致しても、ピークがずれる例を確認できます。
3. いえません。固定Lのピークは2π/tの同値類です。この演習は補助時間の連続位相追跡も使っています。
   FCIは比較にのみ使い、枝選択には渡しません。実際の測定で枝を決めるには追加の測定時間や事前情報が必要です。

</details>''')]
