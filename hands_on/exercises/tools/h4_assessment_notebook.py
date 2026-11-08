"""Exercise 4: preparation/evolution costs, tensor-derived circuits and C-readout."""
import copy
import nbformat as nb


def build_h4_assessment_cells(install, init, modeltext, config, hardwaretext, checks, fit_info):
    cells=[]
    def md(text):cells.append(nb.v4.new_markdown_cell(text.strip()))
    def code(text):cells.append(nb.v4.new_code_cell(text.strip()))
    def answer(text):cells[-1].source+='\n\n<details><summary>解答例</summary>\n\n'+text.strip()+'\n\n</details>'

    md(r'''# 演習4：H₄の回路費用・fidelity・L=1信号を結ぶ

演習3と同じBPDEを、H₄/STO-3Gの**系8＋補助1 qubit**へ拡大します。
まず厳密準備とTrotter回路の費用を確認し、次にMPS/MPOから得た短い回路を評価します。

| 手順 | 確認する対応 |
|---|---|
| 模型と参照状態 | 8スピン軌道、4電子、S₀と最低三重項T₀ |
| 厳密準備・Trotter | 振幅やPauli項数の増加 → 単独部品と全BPDEの費用 |
| MPS/MPO由来の回路 | テンソル表現 → 有限層数の回路へのフィット → 基底ゲート |
| fidelityとL=1曲線 | 状態の重なり、演算子の平均的一致、信号の高さと位相 |
| 測定モード | 厳密確率／指定shots／IBM実機で同じ全0事象を測る |
| 次冊へ | L=1の診断から、Lと探索点を変える反復推定へ |

保存済みの準備M=4・6・8と発展M=8・16・32から選び、重い最適化は再実行しません。
本冊だけで実行できます。前のNotebookの変数・出力ファイルを読み込む必要はありません。
各課題の直後に解答を付けています。
''')
    cells.append(copy.deepcopy(install))
    cells.append(copy.deepcopy(init))
    md('### 回路層数Mを選ぶ\n\n準備はS₀・T₀の両枝に共通の `PREP_M`（4・6・8）、時間発展は独立な `EVOLUTION_M`（8・16・32）を選びます。\n同じ模型で保存済みの9通りを利用します。任意の整数の回路を新規最適化するセルではありません。\n既定値は6と16です。変更したら、この設定セル以降を上から再実行してください。\nCASE=`exact`はMによらない対照、他のCASEは指定した部品を使います。Mと時間反復Lは別です。\n実機モードも同じ選択回路を使います。Mの異なるジョブを再利用せず、保存履歴にもMを記録します。\n演習4と5を比較するときは両冊で同じMを指定します。')
    code("PREP_M = 6         # 4, 6, 8（S0/T0共通）\nEVOLUTION_M = 16   # 8, 16, 32\nselect_circuits(PREP_M, EVOLUTION_M)")
    molecule_intro=copy.deepcopy(modeltext)
    molecule_intro.source=molecule_intro.source.split('| MODE |')[0].rstrip()
    molecule_intro.source=molecule_intro.source.replace('## 共通条件と測定方式','## 共通の分子・回路条件')
    cells.append(molecule_intro)
    md(r'''## 0．H₂からH₄へ：模型の大きさと計算する差

STO-3GのH₄では4空間軌道×2スピンで8スピン軌道です。
削減なしJWは8占有qubit、Fock空間は256次元です。4電子の空間は70次元、Nα=Nβ=2は36次元です。
補助q₈を加え、回路は9 qubitで動かします。

本冊は等間隔直鎖、隣接距離2 Åで、Boys局在化と軌道順の最適化を施した軌道を使います。
H₂の正準軌道や別のH₄データの振幅・回路を混ぜず、**同じ軌道のH・準備・発展**を組み合わせます。
全軌道を活性に含めた同じ模型内の軌道回転はFCIエネルギーを変えませんが、配置振幅・Pauli係数・回路費用は変わります。
符号付きのFock順変換も保存発展回路に含まれています。普通のビットSWAPだけではフェルミオンの符号を表せません。

対象はS₀と最低三重項T₀のM_S=0成分です。Hは電子Hamiltonian、全エネルギーはE+E_nucです。
次セルで基準エネルギー、ギャップ、粒子数・スピンを確認します。
''')
    code('''NA,NB,SPIN2=spin_operators()
print('State  total energy / Ha   Nalpha Nbeta  S^2')
for name,psi,key in [('S0',EXACT_STATES[0],'energy_S0'),('T',EXACT_STATES[1],'energy_T')]:
    energy=float(MODEL[key])
    np.testing.assert_allclose(H@psi,energy*psi,atol=1e-11)
    na,nb,ss=[float(np.vdot(psi,op@psi).real) for op in [NA,NB,SPIN2]]
    np.testing.assert_allclose([na,nb,ss],[2,2,0 if name=='S0' else 2],atol=1e-10)
    print(f'{name:4s} {energy+float(MODEL["enuc"]):18.10f} {na:6.1f} {nb:5.1f} {ss:4.1f}')
print('FCI gap / Ha:',GAP,'DT:',DT,'system dimension:',len(H))
''')
    md('''### 課題1：大きさと対象を確認する

1. H₄が系8 qubitになるのは電子が8個だからか。
2. 本冊のS₀とT₀の差を、H₂のS₀とS₂の差と同じ値だと考えてよいか。
''')
    answer('''1. 電子は4個です。qubit数は8スピン軌道の占有を残すことから決まります。
2. 分子・構造・二つの状態が異なるので同じ値ではありません。本模型の正のT₀−S₀差は約0.01590496 Haです。
状態・軌道・Hamiltonianを揃えた上で、その模型内の参照値と比較します。''')
    md(r'''## 1．厳密準備とTrotter発展は、どこが重くなるか

H₂の固有状態は二配置の重ね合わせで短く準備できました。
H₄では各状態を256成分の振幅として与え、多重化Ry回転で$P_k|0^8\rangle=|\psi_k\rangle$を作ります。
振幅から角度を直接求める**厳密な合成**で、状態探索やVQEではありません。
一般の実振幅準備で扱う係数や回転数はqubit数とともに指数的に増えます。

```text
既知のFCI振幅 ψ(q7,...,q0)
       ↓ 部分ノルムから角度を計算
上位ビットのRy → 条件付きRyの列 → 下位ビットのRy → Pk|0^8> = ψk
```

発展は演習3と同じ、$H=\sum_jc_jP_j$からPauli回転を並べる二次Trotterです。
Pauli項数は次セルで表示します（I項を含む）。一歩Δt=0.1を分割数sで近似します。
小さい係数を既定の相対許容誤差で落とさないよう、展開時はrtol=0を指定し、Hへ戻ることを照合します。
Trotterは厳密なexp(−iHΔt)そのものではありません。分割数を増やすと誤差を減らせますが回路は長くなります。

費用を比較する際は、**個別状態の準備P₀**と、**二つの枝を作る9 qubitのA**を分けます。
同様に、発展一歩と、準備・逆準備・位相補正まで含む全BPDEも分けます。
次セルは同じ基底ゲート・全結合条件でこれらをコンパイルします。
''')
    code('''BASIS=['rz','sx','x','cx']
COSTS={}
COMPILED_CIRCUITS={}

def record_cost(name,qc):
    native=transpile(qc,basis_gates=BASIS,optimization_level=1,seed_transpiler=7)
    row=dict(qubits=qc.num_qubits,cx=int(native.count_ops().get('cx',0)),depth=native.depth())
    COSTS[name]=row
    COMPILED_CIRCUITS[name]=native
    print(f'{name:34s} qubits={row["qubits"]} CX={row["cx"]:5d} depth={row["depth"]:6d}')
    return native

EXACT_PREPS=[prepare_real(psi) for psi in EXACT_STATES]
pauli=SparsePauliOp.from_operator(Operator(H),atol=1e-12,rtol=0)
np.testing.assert_allclose(pauli.to_matrix(),H,atol=1e-12)
TROTTER_SLICES=1  # 1, 2, 4で一歩DTの費用を比較
qc_trotter=QuantumCircuit(8)
qc_trotter.append(PauliEvolutionGate(pauli,time=DT,
    synthesis=SuzukiTrotter(order=2,reps=TROTTER_SLICES)),range(8))
qc_trotter=qc_trotter.decompose()
for name,qc in [('Exact S0 preparation',EXACT_PREPS[0]),('Exact T preparation',EXACT_PREPS[1]),
                ('Exact coherent A',A_EXACT),('Trotter step',qc_trotter)]:
    record_cost(name,qc)
for name,qc,psi in zip(['S0','T'],EXACT_PREPS,EXACT_STATES):
    np.testing.assert_allclose(Statevector.from_instruction(qc).data,psi,atol=1e-11)
print('Pauli terms:',len(pauli),'Trotter slices:',TROTTER_SLICES,'time step:',DT)
''')
    md('''### 課題2：費用表の読み方

1. 表のTrotter stepは厳密な時間発展か。
2. 個別S₀の準備とcoherent AのCX数が違う理由は何か。
3. TROTTER_SLICESを増やすこととLを増やすことは同じか。
''')
    answer('''1. 近似です。厳密準備とTrotter発展を組み合わせても「全て厳密」にはなりません。
2. Aは補助に応じてS₀/Tの両枝を作る回路なので、個別準備とは制御と回転の数が異なります。
3. 分割数は一歩DTの内部、Lは一歩の回路の反復です。前者は総時間DTを保ち、後者は総時間LDTを増やします。
既定の全結合条件では個別厳密準備は各254 CX、coherent Aは510 CXです。
Trotterの費用は表示されたPauli項数・分割数で確認します。係数を落とした別のHと混同しないことが必要です。''')
    md(r'''## 2．MPS/MPOから短い回路へ：表現と回路化を分ける

MPSは状態の振幅を、物理添字qᵢを持つ小さなテンソルの連鎖で表します。
$$\psi_{q_7\cdots q_0}=A^{[7]q_7}A^{[6]q_6}\cdots A^{[0]q_0}.$$
MPOは各サイトに出力・入力の二つの物理添字を持ち、演算子を同様に表します。
$$U_{q'_7\cdots q'_0,q_7\cdots q_0}=W^{[7]q'_7q_7}\cdots W^{[0]q'_0q_0}.$$
隣のテンソルをつなぐ内部添字の最大次元がbond dimensionです。
MPSのχ、回路の層数M、時間反復L、Trotter分割数sはそれぞれ別の量です。

```text
FCI状態 → χ=6のMPS目標 → 指定Mの2-qubitゲート列へフィット → 準備回路Pk
時間発展 → TrotterとMPO積・SVDによる参照MPO → 指定Mゲート列へフィット → 発展回路Uapp
```

テンソルが状態や演算子を正確に表せても、限られた層数の回路が同じ精度になるとは限りません。
特に、MPS/MPOの打切りと、回路へのフィットは別段階の近似です。
参照MPOは圧縮で非ユニタリになり得るため、その行列をそのまま測定回路として使いません。
測定には保存された**ユニタリなゲート列**を使います。

本冊の保存計算では、χ=6 MPSは両状態を数値精度で保持しています。
準備回路は状態fidelityへのフィット、発展回路は参照MPOへの演算子誤差へのフィットです。
発展回路は二次Trotter100分割とSVD cutoff=10⁻¹⁸から得たMPOを目標とする有限予算の保存結果で、完全収束や大域最適解は主張していません。
参照MPOを厳密なexp(−iHDT)と同一視せず、最終的なfidelityと信号は共通のHから得たU_exactで再評価します。
このNotebookでは最適化を再実行せず、実際の保存回路を検証します。

保存回路は隣接2-qubitゲートを交互に並べる層構造です。1個のGは1個の2-qubitユニタリで、CX 1個とは限りません。

```text
層1: G(0,1)  G(2,3)  G(4,5)  G(6,7)
層2:     G(1,2)  G(3,4)  G(5,6)
層3: G(0,1)  G(2,3)  G(4,5)  G(6,7) ...
```

準備・発展ともM層には交互の4・3ブロックがあります（既定M6は21、M16は56ブロック）。
M4/M8の枝回路は小さい成分も保つ厳密な局所分解で合成しています。CX数はMだけで決めず、次の費用表で確認します。
発展のゲート列はinterleavedスピン順で作られており、準備のblocked順へ接続するときの変換も含めて費用を数えます。
''')
    code(f'FIT_INFO={fit_info!r}\n'+'''
print('Preparation M, chi:',SELECTED_PREP_M,PREP_INFO['reports']['S0']['chi'],'evolution M:',SELECTED_EVOLUTION_M)
print('Ground-state MPS bonds:',FIT_INFO['S0_MPS_bonds'])
print('MPO target microsteps/cutoff:',FIT_INFO['mpo_microsteps'],FIT_INFO['mpo_cutoff'])
print('Evolution fit iteration cap:',EVOLUTION_INFO['report']['refinement_cap'],
      'full convergence claimed:',EVOLUTION_INFO['report']['full_convergence_claimed'])
''')
    md(r'''## 3．state fidelityとoperator fidelityは、何を比べるか

状態は一つのベクトルの重なり、演算子は全256次元での一致を評価します。
$$F_{\rm state}^{(k)}=|\langle\psi_k^{\rm FCI}|P_k^{\rm app}|0^8\rangle|^2,\qquad
F_{\rm op}=\left|\frac{\mathrm{Tr}(U_{\rm exact}^\dagger U_{\rm app})}{256}\right|^2.$$
state fidelityはS₀とTを別々に、operator fidelityは一歩DT=0.1に対して計算します。
F_opは全Fock空間の指標で、目的の二状態だけの位相差誤差ではありません。
F_op=1なら大域位相を除いて一致しますが、一般に1に近いだけではギャップの精度を保証しません。

次セルでは**MPS目標のfidelity**と**実際の回路状態のfidelity**を並べます。
準備回路のNα,Nβは2を保ちますが、S²を厳密に固定する最適化ではないことも確認します。
''')
    code('''STATE_FIDELITIES={}
print('State  MPS fidelity         circuit fidelity      energy error/Ha    <S^2>')
for name,exact,approx,key in zip(['S0','T'],EXACT_STATES,APPROX_STATES,['energy_S0','energy_T']):
    fidelity=float(abs(np.vdot(exact,approx))**2)
    STATE_FIDELITIES[name]=fidelity
    target_fid=SELECTED_RESULTS['fidelities'][name+'_MPS']
    energy_error=float(np.vdot(approx,H@approx).real-float(MODEL[key]))
    spin=float(np.vdot(approx,SPIN2@approx).real)
    np.testing.assert_allclose([np.vdot(approx,NA@approx),np.vdot(approx,NB@approx)],[2,2],atol=1e-10)
    np.testing.assert_allclose(fidelity,SELECTED_RESULTS['fidelities'][name],atol=1e-10)
    print(f'{name:4s} {target_fid:.14f} {fidelity:20.12f} {energy_error:18.8f} {spin:10.6f}')
OPERATOR_FIDELITY=float(abs(np.trace(U_EXACT.conj().T@U_APPROX)/256)**2)
print('One-step operator fidelity:',OPERATOR_FIDELITY)
print('One-step relative operator error:',np.linalg.norm(U_APPROX-U_EXACT)/np.linalg.norm(U_EXACT))
np.testing.assert_allclose(U_APPROX.conj().T@U_APPROX,np.eye(256),atol=1e-10)
''')
    md('''### 課題3：テンソルと回路の精度を混同しない

1. χ=6 MPSのfidelityがほぼ1なら、M6準備回路もほぼ1といえるか。
2. operator fidelityがS₀のstate fidelityより高ければ、BPDEのギャップ誤差も小さいといえるか。
3. 準備回路がユニタリで、粒子数が正しければ、参照S₀そのものを作れているか。
''')
    answer('''1. いえません。MPSから限られた層数の回路へフィットする段階が別にあります。
既定M6/M16の例ではMPS目標は両方ほぼ1ですが、回路はS₀≈0.93860、T≈0.999564です。
2. 指標の対象が違うため、その数値だけでは順位付けできません。既定M16のF_op≈0.999999858でも、目的の相対位相を別に調べます。
3. ユニタリ性と粒子数は必要な検証ですが、特定の固有状態への一致を保証しません。重なり・S²・エネルギーなどを確認します。''')
    md(r'''## 4．近似回路の費用：単独準備、枝制御、逆準備を数える

個別の保存準備P₀、P_Tは、配置準備を含めて$P_k|0^8\rangle=|\widetilde\psi_k\rangle$を作ります。
BPDE用には$A=\mathrm{diag}(P_0,P_T)(H_a\otimes I)$を実装し、
$$(|0\rangle_a|\widetilde\psi_0\rangle+|1\rangle_a|\widetilde\psi_T\rangle)/\sqrt2$$
を準備します。同じ入力配置とゲート配置を共有し、対応する局所ゲートG₀の後へ
補助1で制御した$G_TG_0^\dagger$を置くことで、二つの枝を選びます。
実際のAの逆を逆準備に使い、近似状態を作った後だけ厳密な逆回路へ取り替えることはしません。

```text
a=q8: |0> ──┌─────┐────────P(εt)──┌─────┐──M──┐
            │  A  │               │ A†  │     ├─ 全9ビット0を数える
s=q0..7: ───└─────┘──[Uapp]^L──────└─────┘──M──┘
```

既定の発展M16の基底ゲート列だけは168 CXですが、モード順変換を前後に含めると216 CXです。
次セルは実際に本冊で使う部品を同じ条件で再コンパイルし、**同じ役割同士**を比較します。
最後の表では、厳密準備＋Trotterと、指定Mの準備＋発展の全回路を比べます。
前者の発展もTrotter近似で、完全に厳密なexp(−iHt)回路の資源とは区別します。
''')
    code('''for name,qc in [('Approx S0 preparation',load_circuit_asset('S0')),
                ('Approx T preparation',load_circuit_asset('T')),
                ('Approx coherent A',A_APPROX),('Approx A inverse',A_APPROX.inverse()),
                ('Approx evolution with conversion',EVOLUTION_APPROX)]:
    record_cost(name,qc)
expected_pair=np.r_[APPROX_STATES[0],APPROX_STATES[1]]/np.sqrt(2)
np.testing.assert_allclose(Statevector.from_instruction(A_APPROX).data,expected_pair,atol=1e-10)
exact_prep_trotter=QuantumCircuit(9)
exact_prep_trotter.compose(A_EXACT,inplace=True)
exact_prep_trotter.compose(qc_trotter,range(8),inplace=True)
exact_prep_trotter.p(GAP*DT,8)
exact_prep_trotter.compose(A_EXACT.inverse(),inplace=True)
record_cost('Exact prep + Trotter L=1',exact_prep_trotter)
native_both=record_cost(f'M{SELECTED_PREP_M} prep + M{SELECTED_EVOLUTION_M} evolution L=1',experiment_circuit(1,GAP))
''')
    md('''### 課題4：同じ条件で費用を比較する

1. 単独のM6準備60 CXを、そのままBPDE全体の費用と呼べるか。
2. M16が168 CXなのに、本冊の発展ブロックが216 CXになる理由は何か。
3. Lを増やすと、準備もL回になるか。
''')
    answer('''1. 呼べません。枝を作るcoherent A、逆準備、発展、位相補正が必要です。
既定M6/M16ではAとA†は各258 CX、変換込みの発展は216 CXで、L=1の全体は732 CXです。
2. 準備と発展のモード順を合わせる符号付き変換の費用を含めるためです。
3. 各ショットの準備と逆準備は1回ずつ、発展だけがL回です。直列計数はCX(A)+CX(A†)+L×CX(Uapp)です。既定6/16では516+216L CXです。
全体最適化や実機ルーティングにより実際の値は変わるので、これは費用の見積もりです。''')
    md(r'''## 5．四つの条件で、L=1の全0曲線を比較する

近似の原因を切り分けるため、次の4条件を同じH・時間・読出しで比べます。

| CASE | 準備する状態 | 一歩の発展 |
|---|---|---|
| `exact` | FCIの厳密S₀/T | exp(−iHDT) |
| `prep_only` | 指定Mの両状態 | exp(−iHDT) |
| `evolution_only` | FCIの厳密S₀/T | 指定Mの保存ユニタリ |
| `both` | 指定Mの両状態 | 指定Mの保存ユニタリ |

$$a_k=\langle\psi_k|U^L|\psi_k\rangle,\quad
P(0^9)=\tfrac14|a_0+e^{i\varepsilon t}a_T|^2.$$
固定したLでは$a_0,a_T$が定数なので、候補εに対する曲線は
$$P(\varepsilon)=\frac{|a_0|^2+|a_T|^2}{4}
+\frac12\mathrm{Re}[e^{i\varepsilon t}a_0^*a_T]$$
という一本の余弦です。近似状態が多くの固有成分を含んでも、この**固定Lでのε掃引**が多周波数になるわけではありません。
影響は高さ・位相に表れ、Lを変えたときの$a_k(L)$には複数成分の時間依存が現れ得ます。

次セルの左図は周期を見る広い範囲、中央は真値付近、右は厳密信号との差です。
FCIギャップを中心に拡大しているのは回路の診断のためで、未知の差を推定した結果ではありません。
''')
    code('''wide=np.linspace(GAP-np.pi/DT,GAP+np.pi/DT,401)
zoom=np.linspace(GAP-.05,GAP+.05,301)
fig,axes=plt.subplots(1,3,figsize=(14,3.5))
for case in CASES:
    axes[0].plot(wide,probability(case,1,wide),label=case)
    axes[1].plot(zoom,probability(case,1,zoom),label=case)
    if case!='exact':axes[2].plot(zoom,probability(case,1,zoom)-probability('exact',1,zoom),label=case)
for ax in axes:
    ax.axvline(GAP,color='black',ls='--');ax.set_xlabel('Trial gap / Ha');ax.legend(fontsize=7)
axes[0].set(ylabel='P(all zero)',title='Full fringe, L=1')
axes[1].set(ylabel='P(all zero)',title='Near FCI gap, L=1')
axes[1].ticklabel_format(axis='y',style='plain',useOffset=False)
axes[2].set(ylabel='P(case) - P(exact)',title='Signal difference')
plt.tight_layout();plt.show();plt.close(fig)
''')
    md(r'''### 5.1．fidelityと、短時間の高さ・位相を数値で結ぶ

各枝の自己相関の絶対値が1に近ければ、全0信号の最大値も1に近くなります。
しかし自己相関$|\langle\psi|U|\psi\rangle|^2$と、参照固有状態とのfidelityは別の重なりです。
厳密Hでの短時間なら$|\langle\psi|e^{-iHt}|\psi\rangle|^2=1-\mathrm{Var}_{\psi}(H)t^2+O(t^4)$なので、
固有状態との重なりが十分でなくても短時間ではほぼ元に戻り得ます。

$z=a_0^*a_T$とすると、ピーク位置は$\varepsilon_{\rm peak}=-\arg(z)/t$（周期$2\pi/t$）で、
最大値は$(|a_0|+|a_T|)^2/4$です。
次セルはL=1の両者を表示します。ピークは選んだ周期の**一回の干渉位相**であり、Hamiltonianの固有値やBayesian推定値ではありません。
''')
    code('''L1_METRICS={}
print('Case             Pmax          peak/Ha       peak-FCI/Ha    |a0|^2      |aT|^2')
for case,(states,u) in CASES.items():
    a0,at=[np.vdot(psi,u@psi) for psi in states]
    z=a0.conjugate()*at
    peak=-float(np.angle(z))/DT
    pmax=float((abs(a0)+abs(at))**2/4)
    L1_METRICS[case]=dict(peak=peak,pmax=pmax,bias=peak-GAP,
                        ground_survival=float(abs(a0)**2),triplet_survival=float(abs(at)**2))
    np.testing.assert_allclose(probability(case,1,peak),pmax,atol=1e-10)
    print(f'{case:15s} {pmax:.10f} {peak:13.9f} {peak-GAP:+14.9f} {abs(a0)**2:11.8f} {abs(at)**2:11.8f}')
np.testing.assert_allclose(L1_METRICS['both']['peak'],SELECTED_RESULTS['fringe']['apparent_gap'],atol=1e-10)
np.testing.assert_allclose(L1_METRICS['both']['pmax'],SELECTED_RESULTS['fringe']['Pmax'],atol=1e-10)
''')
    cells.append(copy.deepcopy(checks))
    md('''### 課題5：高い信号が意味すること

1. 既定M6/M16でS₀のstate fidelity≈0.9386なのに、L=1のPmax≈0.99997となるのは矛盾か。
2. `both`のL=1のピークは、FCIギャップと一致しているか。
3. operator fidelityが1に近いことと、ギャップの相対誤差が小さいことは同じか。
''')
    answer('''1. 矛盾しません。前者は参照固有状態との重なり、後者は短時間発展後に近似準備状態へ戻る干渉信号です。
短時間では異なる固有成分もまだ大きく位相が離れていないことがあります。
2. 既定M6/M16ではFCI≈0.01590496 Haに対し、bothのピークは約0.00922888 Ha、ずれは約−0.00667608 Haです。
高いPmaxだけから「正しいエネルギー差が得られた」とは判断できません。
3. 同じではありません。小さなギャップでは、小さな位相の偏りも大きな相対誤差になり得ます。
演算子・準備・信号・推定をそれぞれ確認します。''')
    md(r'''## 6．測定モードを選ぶ：CASEとMODEは別の設定

CASEは「どの状態と発展を使うか」、MODEは「その回路をどう測るか」です。

| 設定例 | 得る値 | 含まれる誤差 |
|---|---|---|
| CASE=exact, MODE=exact | 厳密準備・厳密発展の確率 | この模型内の参照 |
| CASE=both, MODE=exact | 選択したMの回路の厳密確率 | 回路近似、shot noiseなし |
| CASE=both, MODE=shots | 二項模擬測定の全0頻度 | 回路近似とshot noise（実機ノイズなし） |
| CASE=both, MODE=ibm | 実機の全0頻度 | 回路近似、shots、実機での影響 |

MODE=shotsでは$n_{0^9}\sim\mathrm{Binomial}(R,P)$、表示値は$n_{0^9}/R$です。
各探索点にSHOTS回を使うので、9点なら合計9×SHOTS回です。SEEDは模擬測定の再現性を決め、実機の乱数は制御しません。
MODE=exactではSHOTSを変えても結果は変わりません。SHOTSを増やしても回路近似の偏りは消えません。

**操作手順：** 次の設定セルのMODE・CASE・SHOTS・SEEDを変更したら、その設定セルから下を再実行します。
これによりmeterと乱数を作り直します。測定セルだけを繰り返す場合は、同じmeterの乱数系列の続きを使います。
初めはMODE=exact・CASE=both、その後MODE=shotsへ進めてください。
''')
    cells.append(copy.deepcopy(hardwaretext))
    md('''### 6.1．実機モードの実行順

1. MODE=ibm、CASE=both、BACKEND_NAMEと必要ならINSTANCEを設定する。
2. RUN_HARDWARE=Falseで設定セルを実行し、L=1の配置後回路の費用を確認する。この状態では測定をスキップする。
3. 送信する場合はRUN_HARDWARE=Trueへ変更し、設定セルと次の測定セルを再実行する。
4. ジョブIDとcountsを保存し、結果を回収する。同じ条件の再実行は保存結果を使う。

実機では選択したMの実回路を使うためCASE=bothに限定します。
MODE=ibmへ変更しても、CASEをexactのまま残さないようにします。
MAX_HARDWARE_Lは発展反復の上限であり、QPU使用量や信号を保証する値ではありません。
実機接続制約を含む2量子ビット操作数は、全結合のCX数と区別してください。
''')
    cells.append(copy.deepcopy(config))
    md('''### 6.2．同じL=1曲線を9点で測る

次セルはL=1固定、広い1周期の9探索点です。線は選んだ回路の理想確率、点はMODEで選んだ測定値です。
9点から細かいピーク位置を精密推定するのではなく、信号が残るか、測定頻度がどう変わるかを診断します。
FCIギャップを中心にするのはこの検証用の既知の模型だからです。
実機モードでは9回路を1ジョブで送り、待ち時間中もジョブIDを保持します。
''')
    code('''if MEASUREMENT_READY:
    L=1
    trial=np.linspace(GAP-np.pi/DT,GAP+np.pi/DT,9)
    measured=meter.measure(CASE,L,trial,'exercise4_L1')
    fig,ax=plt.subplots(figsize=(8,3))
    dense_trial=np.linspace(trial[0],trial[-1],401)
    ax.plot(dense_trial,probability(CASE,L,dense_trial),label='Ideal circuit probability')
    ax.plot(trial,measured,'o',label=MODE)
    ax.set(xlabel='Trial gap / Ha',ylabel='P(all zero)');ax.legend();plt.show();plt.close(fig)
    print('Trial gap / Ha    observed probability')
    for eps,p in zip(trial,measured):print(f'{eps:14.6f} {p:21.8f}')
    print('Source:',meter.records[-1]['source'])
    meter.save('outputs/exercise4_measurements.json')
else:
    print('Preview only; no measurement was submitted.')
''')
    md('''### 課題6：モード切替と測定回数

1. MODE=exactのままSHOTSを128→8192へ変えると、何が変わるか。
2. MODE=shotsで同じ変更をすると、何が変わるか。
3. MODE=ibm、RUN_HARDWARE=Falseでは送信されるか。
4. 同じ保存先で同条件のIBM測定を再実行すると、新しいショットが追加されるか。
''')
    answer('''1. SHOTSを使わないので確率は変わりません。CASE=bothなら回路近似は残っています。
2. 各探索点の頻度の標準偏差√[P(1−P)/SHOTS]が減ります。9点の総shotsは1152から73728へ増えます。
   実機ノイズを加えたシミュレーションにはならず、回路近似の偏りも残ります。
3. 設定セルは接続・コンパイルのプレビューを行いますが、測定セルは送信せずスキップします。
4. 同じ設定・保存先なら既存ジョブ／countsを再利用します。独立した再測定は別の保存先を指定し、元の記録を保持します。''')
    md(r'''## 7．演習5への接続：L=1の診断から反復推定へ

本冊の4条件比較とfidelity評価は、測定へ渡す回路を理解するための診断です。
次冊では**同じH、同じMの選択回路、全0確率、MODE/CASEの意味**を使い、測定後の知識からLと探索点を更新します。

```text
本冊: L=1固定 → 全0曲線・fidelity・費用を診断
次冊: 初期事前 → Lと探索点を選ぶ → 測定 → 事後 → 次の事前 → 繰り返す
```

本冊の9点を次冊の初期事前へ自動的に変換するわけではありません。次冊は独立した初期事前から開始します。
まず次冊のCASE=exact・MODE=exactで更新手順を見て、CASE=bothへ切り替えて本冊と同じ近似回路へ戻ります。
比較では初期事前・SHOTS・SEEDなどを揃えます。各Notebookは独立して設定セルを実行してください。

演習5のGaussian更新は、確率曲線をGaussianにフィットし、事前Gaussianと掛け合わせる著者コードの近似手順です。
演習1・2の二項尤度によるグリッド更新とは異なります。分布が狭くなることと、FCIとの差が小さくなることを別々に見ます。
実機ではL増加で発展回路が長くなります。次のLの費用をpreviewし、上限までの結果でも考察します。

### 課題7：次冊で確認すべきこと

1. 本冊でL=1のピークが高ければ、次冊の推定も正しいギャップへ収束すると断定できるか。
2. 次冊のCASE=exactとCASE=bothの比較で、何を揃えるべきか。
3. L=1の全0確率と、次冊の事後分布は同じ量か。
''')
    # Visible answer at the handoff: avoid relying on a collapsed element at the end.
    cells[-1].source+='''

### 解答例（課題7：演習5への接続）

1. **断定できません。** 短時間の高さは正しい固有位相の保証ではありません。長時間での近似位相、混入成分、推定モデルも確認します。
2. 初期事前、測定方式、shots、seed、探索点数などを揃え、回路近似の有無を変えます。
   反復中のLや探索点は過去の測定に依存するため、同じ初期条件でも途中から異なり得ます。
3. **異なる量です。** 全0確率は回路から予測する測定事象の確率、事後分布は測定後の未知ギャップについての知識を表します。
'''
    return cells
