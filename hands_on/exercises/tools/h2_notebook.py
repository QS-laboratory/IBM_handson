"""Exercise 2: connect configurations, spin, invariant blocks and one-qubit circuits."""
import copy
import nbformat as nb


def build_h2_cells(source, install, ms0_block):
    cells=[]
    def md(text):cells.append(nb.v4.new_markdown_cell(text.strip()))
    def code(text):cells.append(nb.v4.new_code_cell(text.strip()))
    def answer(text):
        cells[-1].source+='\n\n<details><summary>解答例</summary>\n\n'+text.strip()+'\n\n</details>'

    md(r'''# 演習2：H₂の電子配置から1 qubitの位相推定へ

演習1の2×2行列を、分子の電子Hamiltonianへ置き換えます。
**電子配置→スピンと固有状態→4×4のブロック→2×2の符号化→回路→測定**の順に進みます。
核間距離0.735 Å、STO-3G、核固定、中性H₂（2電子・2空間軌道）の小さな模型です。

| 本冊で確認するもの | 対応する計算 |
|---|---|
| S₀・S₁・S₂・T₀の違い | 配置の振幅、全スピン、エネルギー |
| S₀・S₂だけなら1 qubitで扱える理由 | 二重占有配置が張る不変な2次元空間 |
| 電子配置から固有状態を作る方法 | 1 qubitのRy回転 |
| QPE・BPE・BPDEの化学での意味 | 絶対エネルギーとS₂−S₀差の推定 |

全数値を内蔵しているため、このNotebookだけで実行できます。IBMアカウント・PySCFは不要です。
4×4行列は状態の説明に使い、位相推定の系は1 qubitに保ちます。
S₁・T₀の位相推定は、次冊の削減なしJW・系4 qubit版へ進んでから扱います。
''')
    basic=copy.deepcopy(install)
    basic.source='\n'.join(line for line in basic.source.splitlines() if '実機モード' not in line)
    cells.append(basic)
    md('''### 実行準備

次セルはNumPy、Qiskit、描画関数を読み込みます。
`draw`は回路図、`draw_measured`は測定を追加した回路図です。
`Statevector`で厳密確率を計算し、shotsの模擬測定とは区別します。
''')
    cells.append(copy.deepcopy(source.cells[2]))
    md(r'''## 0．二つの空間軌道と四つのスピン軌道

STO-3Gでは各H原子に1個ずつ基底関数があり、RHF計算から二つの空間軌道を作ります。
結合性軌道を$g=\sigma_g$、反結合性軌道を$u=\sigma_u$と書きます。
それぞれにα（↑）、β（↓）のスピンがあり、計4スピン軌道です。
占有ビットの順序を$(q_0,q_1,q_2,q_3)=(g\alpha,u\alpha,g\beta,u\beta)$、
表示ケットを$|q_3q_2q_1q_0\rangle$とします。

```text
                 Dg            Du            O1            O2
u (反結合性)    [    ]        [↑↓ ]         [ ↑ ]         [ ↓ ]
g (結合性)      [↑↓ ]         [    ]        [ ↓ ]         [ ↑ ]
```

| 配置名 | 意味 | 占有ケット |
|---|---|---|
| $D_g$ | gを二重占有（RHF配置） | $|0101\rangle$ |
| $D_u$ | uを二重占有（二重励起配置） | $|1010\rangle$ |
| $O_1$ | gβとuαを一つずつ占有 | $|0110\rangle$ |
| $O_2$ | gαとuβを一つずつ占有 | $|1001\rangle$ |

いずれも電子数2、$M_S=(N_\alpha-N_\beta)/2=0$です。
この4配置のほかにgα・uα、gβ・uβの配置があり、電子数2の空間全体は$\binom42=6$次元です。
各ケットは反対称化済みのSlater行列式です。電子を区別できる粒子として並べた積ではありません。
「配置」は占有の基底、「分子の固有状態」はその線形結合です。
''')
    md(r'''## 1．S₀・S₁・S₂・T₀は何が違うか

全スピン$S$に対し$\hat S^2$の固有値は$S(S+1)$（$\hbar=1$）です。
一重項はS=0、三重項はS=1。S₀/S₁/S₂の添字は**一重項の中でのエネルギー順位**です。
T₀は本教材で最低三重項を指し、計算ではその$M_S=0$成分を使います（コード名は`T`）。
最低三重項をT₁と呼ぶ流儀もあります。ここでは「三重項の順位」と「スピン射影M_S」を分けて考えます。

$$\begin{aligned}
|S_0\rangle&=\cos(\theta/2)|D_g\rangle+\sin(\theta/2)|D_u\rangle,\\
|S_2\rangle&=-\sin(\theta/2)|D_g\rangle+\cos(\theta/2)|D_u\rangle,\\
|S_1\rangle&=(|O_1\rangle+|O_2\rangle)/\sqrt2,\\
|T_0,M_S=0\rangle&=(|O_1\rangle-|O_2\rangle)/\sqrt2.
\end{aligned}$$

ここでθは後でHamiltonianから決まる混合角です。
最後の二つの±は、上で固定した行列式・軌道の位相規約によるものです。
電子1・2のスピン関数で書く一重項$(\alpha\beta-\beta\alpha)/\sqrt2$の符号を、
反対称化を含む行列式の係数へそのまま移してはいけません。後で$S^2=0,2$を数値照合します。

| 状態 | 配置の特徴 | 全スピン | この距離でのエネルギー順位 |
|---|---|---:|---|
| S₀ | 主にDg、少量のDu | 0 | 最低 |
| T₀ | O1とO2、相対符号− | 1 | S₀より上、S₁より下 |
| S₁ | O1とO2、相対符号＋ | 0 | 三重項とは別の状態 |
| S₂ | 主にDu、少量のDg | 0 | 本模型の最高一重項 |

三重項には$M_S=-1,0,+1$の3成分があります。スピンに依存しないこのHでは同じエネルギーです。
$M_S=0$の三重項を「二電子がどちらも↑の配置」と説明することはできません。

### 課題1：占有とスピンを区別する

1. O1、O2はそれぞれS₁またはT₀そのものか。
2. S₁とT₀を占有ビットの測定確率だけで区別できるか。
3. 同じM_S=0なら同じSか。
''')
    answer('''1. それぞれは行列式であり、この場合はS₁とT₀の重ね合わせです。固有状態はO1±O2です。
2. 両状態ともO1、O2が各1/2で、占有確率だけでは区別できません。相対符号とS²の測定に違いがあります。
3. M_S=0でもS=0とS=1を含みます。S₁はS²=0、T₀はS²=2です。''')
    md(r'''## 2．配置を基底にした4×4行列で違いを見る

基底を$(D_g,D_u,O_1,O_2)$の順に並べます。電子Hamiltonianは、丸めて
$$H_{M_S=0}=\begin{pmatrix}
-1.836968&0.180931&0&0\\
0.180931&-0.245218&0&0\\
0&0&-1.063653&0.180931\\
0&0&0.180931&-1.063653
\end{pmatrix}\quad\text{Hartree}.$$

左上が二重占有配置、右下が開殻配置のブロックです。
左上の非対角項がHF配置DgとDuを混ぜ、S₀・S₂を作ります。
右下の対角要素は等しいので、固有ベクトルはO1±O2です。
対称結合のS₁は$-1.063653+0.180931$、反対称結合のT₀は$-1.063653-0.180931$の電子エネルギーです。

この基底で$S^2$は左上が0、右下が$\begin{pmatrix}1&-1\\-1&1\end{pmatrix}$です。
次セルでは二つのブロックを対角化し、エネルギー固有方程式とS²を照合します。
図の**符号付き振幅**と表のエネルギーを見てください。全エネルギーには核間反発0.7199689944 Haを加えます。
''')
    code(f'''H_MS0 = np.array({ms0_block!r}, dtype=float)
# 独立したJW教材から取り出した中性・Ms=0ブロック。微小なゼロ要素は保持する。
S2_MS0 = np.zeros((4,4))
S2_MS0[2:,2:] = [[1,-1],[-1,1]]
energies_closed, vectors_closed = np.linalg.eigh(H_MS0[:2,:2])
for k in range(2):
    if vectors_closed[np.argmax(abs(vectors_closed[:,k])),k] < 0:
        vectors_closed[:,k] *= -1
STATES_MS0 = {{'S0':np.r_[vectors_closed[:,0],0.,0.],
              'S1':np.array([0,0,1,1])/np.sqrt(2),
              'S2':np.r_[vectors_closed[:,1],0.,0.],
              'T':np.array([0,0,1,-1])/np.sqrt(2)}}
ENERGIES_MS0 = {{name:float(psi@H_MS0@psi) for name,psi in STATES_MS0.items()}}
print('State   total energy / Ha   S^2    amplitudes (Dg,Du,O1,O2)')
for name,psi in STATES_MS0.items():
    energy=ENERGIES_MS0[name]
    spin=float(psi@S2_MS0@psi)
    np.testing.assert_allclose(H_MS0@psi,energy*psi,atol=1e-12)
    np.testing.assert_allclose(S2_MS0@psi,(2 if name=='T' else 0)*psi,atol=1e-12)
    print(f'{{name:4s}} {{energy+.7199689944489797:18.10f}} {{spin:5.1f}}',psi)
fig,axes=plt.subplots(1,2,figsize=(11,3.5))
names=['S0','S1','S2','T']
amplitudes=np.array([STATES_MS0[name] for name in names])
im=axes[0].imshow(amplitudes,cmap='coolwarm',vmin=-1,vmax=1)
axes[0].set_xticks(range(4),['Dg','Du','O1','O2']);axes[0].set_yticks(range(4),names)
axes[0].set_title('Signed configuration amplitudes')
for i in range(4):
    for j in range(4):axes[0].text(j,i,f'{{amplitudes[i,j]:+.3f}}',ha='center',va='center',fontsize=9)
fig.colorbar(im,ax=axes[0],shrink=.8)
for i,name in enumerate(['S0','T','S1','S2']):
    energy=ENERGIES_MS0[name]+.7199689944489797
    axes[1].plot([i-.25,i+.25],[energy,energy],lw=3)
    axes[1].text(i,energy+.04,f'{{energy:.4f}}',ha='center',fontsize=9)
axes[1].set_xticks(range(4),['S0','T0','S1','S2'])
axes[1].set(ylabel='Total energy / Ha',title='Energy order at R = 0.735 A',ylim=(-1.3,.7))
plt.tight_layout();plt.show();plt.close(fig)
''')
    md(r'''## 3．S₀・S₂だけなら、なぜ1 qubitで足りるか

H₂の反転対称性により、Dg・Duは偶、O1・O2は奇の空間対称性を持ちます。
g×gとu×uは偶、g×uは奇なので、この固定模型のHは両ブロックを結びません。
S₁も一重項ですが**別の空間対称性のブロック**にあります。一重項という条件だけでは除外できません。

二重占有空間$\mathcal V=\mathrm{span}\{D_g,D_u\}$は、Hを作用させても外へ出ない不変部分空間です。
その2次元だけを、論理qubitの2次元へ対応させます。
$$V|0\rangle_L=|D_g\rangle,\quad V|1\rangle_L=|D_u\rangle,\quad
H_L=V^\dagger H_{M_S=0}V.$$
$V=\begin{pmatrix}1&0\\0&1\\0&0\\0&0\end{pmatrix}$とすれば
$H_{M_S=0}V=VH_L$であり、$e^{-iH_{M_S=0}t}V=Ve^{-iH_Lt}$です。
このブロック内では厳密な符号化です。電子間相互作用を捨てたり、片方の電子を消したりしているわけではありません。

```text
4スピン軌道のFock空間: 16次元   (削減なしJWでは4 qubit)
             ↓ 電子数2: 6次元
             ↓ Ms=0: 4次元 = {Dg,Du} ⊕ {O1,O2}
             ↓ 二重占有の不変ブロックを選択
          2次元 = {Dg,Du}  ↔  {|0>_L, |1>_L}: 1 qubit
```

演習3への接続に使う4占有ビットでは、$|0\rangle_L\leftrightarrow|0101\rangle$、
$|1\rangle_L\leftrightarrow|1010\rangle$です。元の自習教材はParity変換を介して同じブロックを得ます。
本冊では上の等長写像Vで直接対応を示します。

### 課題2：表現できる状態の範囲

1. 論理$\alpha|0\rangle_L+\beta|1\rangle_L$は、分子ではどの状態か。
2. 論理|1⟩を準備すればS₁またはS₂になるか。
3. この1 qubitでS₁・T₀も含む四つの直交状態を同時に表現できるか。
''')
    answer('''1. αDg+βDuです。任意の二電子状態ではなく、選んだ2配置内の状態です。
2. |1⟩Lは配置Duです。S₁は別ブロック、S₂は主にDuですがDgも混ざるため、どちらそのものでもありません。
3. 1 qubitの次元は2なので、四つの直交状態を同時には保持できません。
本冊のS₀/S₂の符号化を使ってS₁/T₀を計算することはできません。
4状態の抽象的な符号化なら2 qubitで可能ですが、次冊では全Fock空間を残すJWの4 qubit表現を使います。''')
    md(r'''### 3.1．2×2 Hamiltonianと位相の変数を定義する

以下のHはVで取り出した電子Hamiltonianです。
$$H=cI-zZ+xX,\quad c=\tfrac12(H_{00}+H_{11}),\quad
z=\tfrac12(H_{11}-H_{00}),\quad x=H_{01}.$$
$\rho=\sqrt{x^2+z^2}$、$\theta=-\operatorname{atan2}(x,z)$とすると、
$R_y(\theta)^T H R_y(\theta)=\operatorname{diag}(c-\rho,c+\rho)$です。

| コード | 意味 |
|---|---|
| `zero`, `one` | 配置Dg、Duの論理表現 |
| `g`, `e` | S₀、S₂の論理振幅 |
| `E`, `ENUC` | 電子固有エネルギー、核間反発 |
| `DT`, `L` | 基本時間Δt=1原子単位、反復数 |
| `PHASES` | −EΔt/(2π)で得る位相 |
| `GAP` | 位相差を周期1で表した値 |
| `ENERGY_GAP` | 正のS₂−S₀差、Hartree |

位相とエネルギーのGAPは同じ数ではありません。符号と周期の換算は第5節で確認します。
次セルではHを対角化し、HFの重なりと参照値を表示します。
''')
    definition=source.cells[4].source
    code(definition[:definition.index('prep_g =')])
    md(r'''## 4．配置からS₀・S₂を作る回路

$R_y(\theta)=\begin{pmatrix}\cos(\theta/2)&-\sin(\theta/2)\\
\sin(\theta/2)&\cos(\theta/2)\end{pmatrix}$なので、二つの列が固有状態です。

```text
S0: |0>_L ─────────Ry(θ)──→ cos(θ/2)|0> + sin(θ/2)|1>
S2: |0>_L ──X──────Ry(θ)──→ -sin(θ/2)|0> + cos(θ/2)|1>
```

θ≈−0.223537で、S₀は約$0.99376D_g-0.11154D_u$です。
HFはDg一配置で、S₀と同じではありません。非対角項があるためH|Dg⟩はDu成分を持ちます。
重なり$|\langle S_0|HF\rangle|^2\approx0.98756$、残り約0.01244はS₂です。

次セルで回路振幅、4配置への埋め込み、固有方程式を照合します。
1 qubitの回転は**二つの電子配置の係数を変える**操作で、実在の片方の電子スピンを回すという意味ではありません。
''')
    code(definition[definition.index('prep_g ='):]+'''
V=np.eye(4)[:,:2]
np.testing.assert_allclose(V.T@H_MS0@V,H,atol=1e-12)
np.testing.assert_allclose(H_MS0@V,V@H,atol=1e-12)
for name,qc,psi in [('S0',prep_g,g),('S2',prep_e,e)]:
    actual=Statevector.from_instruction(qc).data
    np.testing.assert_allclose(actual,psi,atol=ATOL)
    np.testing.assert_allclose(V@actual,STATES_MS0[name],atol=ATOL)
    energy=float(np.vdot(psi,H@psi).real)
    np.testing.assert_allclose(H@psi,energy*psi,atol=ATOL)
    print(name,'logical amplitudes:',actual,'configuration amplitudes:',V@actual)
print('HF variance:',float(np.vdot(zero,H@H@zero).real-(H[0,0])**2),'Ha^2')
''')
    md(r'''### 課題3：配置・固有状態・回路を結ぶ

1. S₂を作る回路にXだけを残すと、どの配置を作るか。
2. HFのエネルギー期待値とS₀の固有エネルギーは同じか。
3. 二つの準備回路の出力が直交することを、Ryの列または内積で確認する。
''')
    answer('''1. 論理|1⟩、つまりDuを作ります。RyがないのでS₂のDg成分を作れません。
2. HF全エネルギーは約−1.116999 Ha、S₀は約−1.137306 Haで、異なります。
HFは固有状態ではなく分散も非零です。混合によって模型内の基底エネルギーが下がります。
3. Ryはユニタリなので列が直交し、⟨S₀|S₂⟩=0です。二つを4配置へ写してもV†V=Iにより内積は変わりません。''')
    md(r'''## 5．時間発展と干渉：行列版との符号の違い

演習1は$U=e^{+2\pi iK}$、分子では$U=e^{-iH\Delta t}$です。
$$\phi_j=-E_j^{\rm el}\Delta t/(2\pi)\pmod1,\qquad t=L\Delta t.$$
本模型は既知の電子エネルギー探索区間−2<E_el<0 Haにあり、Δt=1では
$E_{\rm total}=-2\pi\phi/\Delta t+E_{\rm nuc}$で周期をほどけます。

| 状態 | 全エネルギー / Ha | 位相（1周単位） |
|---|---:|---:|
| S₀ | −1.1373060358 | 約0.29559450 |
| S₂ | ＋0.4950577416 | 約0.03579574 |

S₂の方が高エネルギーでも位相は小さくなります。
この1 qubitでは発展を厳密に
$$U^L=e^{-ict}R_y(\theta)R_z(-2\rho t)R_y(-\theta)$$
と作れます。回路は右端の演算子から実行します。

```text
s: ──Ry(-θ)──Rz(-2ρt)──Ry(θ)──   (全体位相 exp(-ict) も保持)

a: |0> ──H──●──────P(β)──H──M
            │
s: |Sj> ───[U^L]─────────────
```

干渉確率は$[1+\cos(2\pi L\phi_j+\beta)]/2$。
位相候補εならβ=−2πLε、**電子エネルギー候補**ε_Eならβ=ε_E tです。
制御Uではexp(−ict)が補助の一方の枝だけに付き、相対位相になるため保持します。
次セルは回転回路と厳密行列、干渉確率と余弦式を照合します。
''')
    cells.append(copy.deepcopy(source.cells[7]))
    md(r'''### 課題4：発展時間と定数項

1. 候補位相ε∈[0,1)を掃引し、L=1、2、4にするとピークはいくつか。
2. 制御UからcIを省いた場合、BPEで得る値は正しい分子エネルギーか。
3. U⁵は位相の第5桁を読む操作か。
''')
    answer('''1. 順に1、2、4個です。候補の周期は1/Lで、長時間だけでは周期の別候補が残ります。
2. H−cIの固有エネルギーを読むことになり、そのままでは正しい絶対エネルギーではありません。
両枝へ共通の時間発展をかけるBPDEでは共通位相は差から消えます。
3. L=5回、t=5Δtの時間発展であり、第5桁ではありません。''')
    md(r'''## 6．QPE：同じ回路から、分子の固有エネルギーを読む

```text
補助m個: |0...0> ──H⊗m──制御 U^(2^r)──IQFT──測定 y → φ=y/2^m → E_total
系1個  : |Sj>   ─────────Uの作用───────────
```

演習1と同じ分布$P(y|\phi)=|2^{-m}\sum_{k=0}^{2^m-1}e^{2\pi ik(\phi-y/2^m)}|^2$です。
分子の位相は有限二進数で正確に表せないため、厳密準備でも出力が周囲の格子へ広がります。
HF入力なら$P(y)=w_0P(y|\phi_0)+w_2P(y|\phi_2)$で、HFの期待値が一つ出るわけではありません。

次セルの`QPE_BITS`を変え、厳密S₀・HF・S₂成分25%の入力を比較します。
系は1 qubitのまま、補助の数と制御Uのべきが増えます。
図の赤い点線はFCIの模型内の参照エネルギーで、実験値ではありません。
''')
    qpe=source.cells[10].source
    code(qpe[:qpe.index('M=8;')]+'''
QPE_BITS=6  # 3, 6, 8を比較
M=2**QPE_BITS;xx=np.arange(M)
np.testing.assert_allclose(Operator(inverse_qft(QPE_BITS)).data,
    np.exp(-1j*TAU*np.outer(xx,xx)/M)/np.sqrt(M),atol=ATOL)
draw_measured(qpe_circuit(g,m=3),[0,1,2])  # 図は読みやすい3補助版
fig,axes=plt.subplots(1,3,figsize=(13,3.2))
for ax,state,label in zip(axes,[g,zero,psi_mixed],['Exact S0','HF configuration Dg','25% S2 component']):
    p=qpe_distribution(state,m=QPE_BITS)
    w=abs(R.T@state)**2
    expected=sum(wj*qpe_formula(phi,QPE_BITS) for wj,phi in zip(w,PHASES))
    np.testing.assert_allclose(p,expected,atol=ATOL)
    energies=phase_to_total(np.arange(M)/M)
    ax.bar(energies,p,width=TAU/(M*DT)*.8)
    for Ej in E+ENUC:ax.axvline(Ej,color='C3',ls=':',label='Exact energy' if Ej==E[0]+ENUC else None)
    ax.axvline(H[0,0]+ENUC,color='black',ls='--',label='HF expectation')
    ax.set(xlim=(-1.6,1),xlabel='Reported total energy / Ha',ylabel='Probability',title=label)
axes[0].legend(fontsize=7);plt.tight_layout();plt.show();plt.close(fig)
print('HF weights S0,S2:',HF_WEIGHTS)
print('Energy grid spacing / Ha:',TAU/(M*DT))
''')
    md('''### 課題5：QPEの確率を解釈する

1. HF入力の約1.24%の成分はS₁か。
2. QPE_BITSを増やせばHFが厳密S₀へ変わるか。
3. 厳密S₀なのに有限mで参照値ぴったりの出力にならない理由は何か。
''')
    answer('''1. 残りは同じブロックのS₂です。S₁・T₀の振幅はこの符号化にはありません。
2. 元のHFの成分の重みは変わりません。増えるのは位相・エネルギー格子の分解能です。
十分な分解能で補助の結果に条件付ければ対応固有状態に近づきますが、無条件の全入力がS₀になる操作ではありません。
3. 分子の位相がy/2^mの格子に載らないためです。時間発展の近似誤差や準備誤差とは別の有限ビット効果です。''')
    md(r'''## 7．BPE：S₀の絶対エネルギーを事後分布から推定する

回路は第5節の1補助干渉測定です。各ショットでS₀を準備し直します。
$$p_{\rm new}(\phi)\propto p_{\rm old}(\phi)
 P_0(\phi)^{n_0}[1-P_0(\phi)]^{R-n_0},\quad
P_0=[1+\cos(2\pi L\phi+\beta)]/2.$$

```text
S0準備 → 干渉測定(L,β)をR shots → n0 → 事前×二項尤度 → 正規化 → 次の事前
```

次セルは演習1と同じ更新関数です。周期を観察するためグリッドは[0,1)全体に取り、
物理的な探索区間で事前を切っていません。`phase_to_total`による換算は上で選んだエネルギーの枝に従います。
`data`は(L,β,R,n₀,厳密確率)、`posteriors`は各設定後の確率質量です。
''')
    bpe=source.cells[16].source
    code(bpe[:bpe.index('settings=')])
    md('''### 7.1．測定countsと更新履歴を見る

SHOTSは各(L,β)での測定回数です。二つのβを使い、短時間から長時間へ進めます。
左図は位相全体、右図は正解付近の拡大です。表のcountsが尤度へ入る値です。
真値は模擬データの生成と照合のみに使い、推定グリッドの正解を固定しません。
''')
    code('''SHOTS=256
SEED=20261108
L_VALUES=[1,2,4,8,16,32]
settings=[(L,beta) for L in L_VALUES for beta in [0.,-np.pi/2]]

def plot_updates(grid,posteriors,truth,title):
    fig,axes=plt.subplots(1,2,figsize=(11,3))
    scale=1/(grid[1]-grid[0])
    axes[0].plot(grid,np.ones_like(grid),':',label='Uniform prior')
    for i in sorted(set([0,len(posteriors)//2,len(posteriors)-1])):
        axes[0].plot(grid,posteriors[i]*scale,label=f'After {i+1} settings')
    axes[1].plot(grid,posteriors[-1]*scale,label='Final posterior')
    for ax in axes:
        ax.axvline(truth,color='black',ls='--',label='Exact phase')
        ax.set(xlabel='Phase or gap / cycles',ylabel='Posterior density');ax.legend(fontsize=8)
    axes[0].set_xlim(0,1);axes[1].set_xlim(truth-.02,truth+.02)
    fig.suptitle(title);plt.tight_layout();plt.show();plt.close(fig)

def print_data(data):
    print(' L  beta/pi  shots  observed n0  exact P0')
    for L,beta,R,n0,p0 in data:print(f'{L:2d} {beta/np.pi:8.2f} {R:6d} {n0:12d} {p0:9.4f}')

post_exact,data_exact=bayes_update(GRID,settings,lambda L,b:bpe_probability([1,0],L,b),shots=SHOTS,seed=SEED)
plot_updates(GRID,post_exact,PHASES[0],'BPE: exact molecular S0')
print_data(data_exact)
summary=posterior_summary(GRID,post_exact[-1])
print(summary,'MAP total energy / Ha:',phase_to_total(summary['MAP']))
''')
    md('''### 課題6：BPEの対象と統計を変える

1. SHOTSを64、256、1024へ変えると、何が変わるか。
2. 生成確率の重みを`[1,0]`から`[0,1]`へ変えると、何のエネルギーを推定するか。
3. `[1,0]`から`HF_WEIGHTS`へ変えたデータを、同じ単一位相尤度で正しく扱えるか。
''')
    answer('''1. 測定頻度の統計変動は一般に1/√SHOTSのスケールで減り、単一位相モデルが正しければ事後は狭くなります。
2. 入力をS₂へ変えた計算になり、位相約0.03579574、全エネルギー約+0.49505774 Haを推定します。
3. HFはS₀・S₂の二周波数の信号です。一つの固有位相だけの尤度とは一般に一致しません。
shotsを増やしてもモデルの不一致は除けません。比較後は重みを[1,0]へ戻してください。''')
    md(r'''## 8．BPDE：S₀とS₂の差を直接測る

論理0/1を分子固有状態へ回す前に補助と相関させます。

```text
準備W_B:
a: |0> ──H──●──────────
            │
s: |0> ─────X──Ry(θ)───

実験全体:
a: |0> ──┌─────┐────────P(β)──┌──────┐──M_a
         │ W_B │              │ W_B† │
s: |0> ──└─────┘──[U^L]───────└──────┘─────
```

準備状態は$(|0\rangle_a|S_0\rangle_s+|1\rangle_a|S_2\rangle_s)/\sqrt2$です。
制御なしの系発展と逆準備で、
$$P_B(a=0)=[1+\cos(2\pi L\Delta+\beta)]/2
=[1+\cos((\varepsilon_E-\Delta E)t)]/2\quad(\beta=\varepsilon_Et).$$
位相差と正のエネルギー差の関係は
$$\Delta=(\phi_{S_2}-\phi_{S_0})\bmod1\approx0.74020124,\qquad
\Delta E=2\pi(1-\Delta)/\Delta t\approx1.63236378\ {\rm Ha}.$$
この換算は0<ΔE<2 Haの既知の探索枝を使います。ΔとΔEを混同しないでください。
cIや核間反発は両枝に共通で、差からは消えます。

次セルで準備と読出しを定義します。既定eta=0が厳密準備で、CNOTとcontrolled-Ry(π)の両方を照合します。
eta>0は励起側にS₀成分を混ぜる診断です。その場合、単一ギャップ尤度が正しいとは限りません。
''')
    bpde=source.cells[23].source
    code(bpde[:bpde.index('post_gap,')])
    md('''### 8.1．BPDEの二項更新を、エネルギー差へ換算する

BPEと同じSHOTS・L_VALUESと一様事前を使い、変数をφからΔへ置き換えます。
次セルはcounts、事後分布、MAPの正のエネルギー差を表示します。
原点からの絶対エネルギーを二つ別々に推定して引く計算とは異なります。
''')
    code('''post_gap,data_gap=bayes_update(GRID,settings,lambda L,b:bpde_probability(0,L,b),shots=SHOTS,seed=83)
plot_updates(GRID,post_gap,GAP,'BPDE B: S2-S0 phase difference')
print_data(data_gap)
print('Phase gap summary:',posterior_summary(GRID,post_gap[-1]))
print('MAP energy gap / Ha:',phase_to_gap(posterior_summary(GRID,post_gap[-1])['MAP']))
''')
    md(r'''### 課題7：BPDEの測定対象

1. この回路でS₁−S₀やT₀−S₀差を得られるか。
2. H→H+dIでBPEとBPDEはそれぞれどう変わるか。
3. eta=1では何の差を測るか。
''')
    answer('''1. 得られません。二つの枝に準備しているのはS₀とS₂であり、S₁/T₀は論理空間に含まれません。
2. BPEの位相は−dΔt/(2π)だけずれ、BPDEでは両枝の共通位相が消えて差は不変です。
3. 両枝がS₀となり、差0の信号です。近似状態の混入は差の直接測定でも自動的には消えません。''')
    md(r'''## 9．Fig. 1C型：全0への戻りを測る

演習1と同じ説明用$W_C=W_B\,\mathrm{SWAP}$を使います。
初期00への作用は同じでも、逆準備の最後のSWAPで情報の出る場所が変わります。
厳密準備eta=0のこの例は補助だけでは常に0、**全0の確率**は位相差に応じて変わります。

```text
|00> → W_C → [P_a(β) と U_s^L] → W_C† → M_a,M_s → 全0の回数 / 全shots
```

$A_j=\langle\widetilde S_j|U^L|\widetilde S_j\rangle$とすると
$$P_C(00)=\tfrac14|A_0+e^{i\beta}A_2|^2.$$
厳密準備ならBと同じ余弦、近似準備なら複数の固有位相の影響が残り得ます。
全0読出しは成功した部分だけで再正規化する操作ではありません。
次セルはB/Cの回路確率、補助のみと全0の読出し、近似準備の影響を比較します。
SWAPは原著の回路そのものではなく、読出しを示す例です。
''')
    cells.append(copy.deepcopy(source.cells[28]))
    md('''### 9.1．全0countsでも同じ差を推定する

次セルはC回路の厳密な全0確率から模擬測定を生成します。
更新関数のn₀はここではn₀₀（全0の回数）です。分母は全shotsです。
B/Cで同じ理想信号を測ることを照合し、Cの事後も正のS₂−S₀差へ換算します。
''')
    code('''for L,beta in settings:
    np.testing.assert_allclose(readout(return_circuit(L=L,beta=beta,kind='C'))[1],
        bpde_probability(0,L,beta),atol=ATOL)
post_c,data_c=bayes_update(GRID,settings,lambda L,b:readout(return_circuit(L=L,beta=b,kind='C'))[1],shots=SHOTS,seed=91)
plot_updates(GRID,post_c,GAP,'BPDE C: all-zero counts, S2-S0')
print_data(data_c)
print('C MAP energy gap / Ha:',phase_to_gap(posterior_summary(GRID,post_c[-1])['MAP']))
''')
    md('''### 課題8：同じ状態と同じ読出しを区別する

1. W_BとW_Cの準備状態が同じなら、逆準備後の補助0確率も同じか。
2. この例で全0を補助0の回数に置き換えてベイズ更新してよいか。
3. 全0測定は、近似状態を厳密固有状態へ直す操作か。
''')
    answer('''1. 同じではありません。時間発展後の状態に対する逆回路の作用は、準備回路全体に依存します。
2. eta=0ではCの補助は常に0で、差の情報がありません。全0という事象に対応するcountsと尤度を揃える必要があります。
3. 近似準備状態への戻りを測るので、厳密状態への浄化ではありません。eta>0では信号と単一ギャップモデルの不一致が残ります。''')
    md(r'''## 10．演習3への橋渡し

| 本冊で使った対応 | 次冊で変わる点 |
|---|---|
| Dg/Duの2配置 ↔ 論理1 qubit | 4スピン軌道の占有を削減なしJWで表す系4 qubit |
| RyでS₀/S₂を厳密準備 | S₁/T₀も厳密準備し、別のギャップを測る |
| Ry・Rzで厳密な1 qubit発展 | Pauli文字列の指数関数、Trotter分解を導入 |
| 補助1個と逆準備、全0測定 | 同じBPDE原理を系4＋補助1 qubitへ拡大 |

次冊では候補変数を位相差Δから正のエネルギー差ΔEへ切り替えます。
本冊の`GAP`は周期1で表す位相差、`ENERGY_GAP`はHa単位の差です。次冊の`GAPS`はHa単位です。
補正位相はβ=−2πLε（位相差候補）からβ=ε_E LΔt（エネルギー差候補）へ換算します。
さらに基本時間は本冊のΔt=1から次冊の0.1へ変わるので、Lだけでなくt=LΔtを揃えて比較します。

### 課題9：1 qubit化で残したもの・残さなかったもの

1. 1 qubit化しても電子数2とDg/Du間の電子相関は保たれるか。
2. S₁/T₀もこの1 qubit空間に含まれるか。
3. 何が厳密で、どの模型上の制約が残るか。
''')
    cells[-1].source+='''

### 解答例（課題9：演習3への橋渡し）

1. **保たれます。** 論理0/1はどちらも2電子配置です。DgとDuの混合を残すので、相関も含みます。
2. **含まれません。** S₁/T₀は開殻配置O1/O2の別ブロックです。次冊では4スピン軌道の占有を4 qubitで表し、これらも準備します。
3. **選んだブロック内の固有状態・時間発展・1 qubitへの写像が厳密です。** 分子全体の全状態を表す符号化ではなく、核固定・STO-3G・固定構造という模型上の制約は残ります。
'''
    return cells
