"""Exercise 5: expose every step of the Gaussian-fit adaptive measurement loop."""
import copy
import nbformat as nb
from ibm_setup import ibm_setup_cells


def build_h4_bayes_cells(install, init, modeltext, config, hardwaretext, checks):
    cells=[]
    def md(text):cells.append(nb.v4.new_markdown_cell(text.strip()))
    def code(text):cells.append(nb.v4.new_code_cell(text.strip()))
    def answer(text):cells[-1].source+='\n\n<details><summary>解答例</summary>\n\n'+text.strip()+'\n\n</details>'

    md(r'''# 演習5：H₄で測定計画とGaussian更新を反復する

演習4と同じH・全0測定を使い、準備Mと発展Mを独立に選びます（既定6・16）。
本冊では、**現在の推定から次の測定を選び、測定後に推定を更新する**過程を追います。
最初に1ラウンドを段階ごとに実行し、その後、同じ処理をループで繰り返します。
ラウンド数を2に固定しません。まず5ラウンド、次に回数や初期分布・shotsを変えて比較します。

```text
事前 Gaussian(μ,v)
    ↓ μ,vから時間LΔtと探索点εを選ぶ
BPDEを各εで測定 → 全0確率または頻度
    ↓ 確率曲線をGaussianにフィット → 中心m、幅²w
事前Gaussian × フィットGaussian → 事後 Gaussian(μ',v')
    ↓ 次の事前(μ,v)へそのまま渡す
次ラウンド
```

| 段階 | コード | 各ラウンドの図 |
|---|---|---|
| 事前→測定計画 | `choose_measurement_settings` | A：事前・探索窓・探索点 |
| 測定→フィット | `meter.measure`, `fit_gaussian_signal` | B：観測とフィット、C：残差 |
| フィット→事後 | `update_gaussian_prior` | D：事前・正規化したフィット・事後 |
| 次ラウンド | `iterative_bayes` | 表の事後が次の行の事前と一致 |

著者コードのGaussianフィットとGaussian積の更新を参照します。
**演習1・2の二項尤度をそのまま使うグリッド更新とは異なる近似手順**です。
全ての問いに解答を付けています。Notebook単独で動き、最適化は再実行しません。
''')
    cells.append(copy.deepcopy(install))
    cells.extend(ibm_setup_cells())
    cells.append(copy.deepcopy(init))
    md('### 回路層数Mを選ぶ\n\n準備はS₀・T₀の両枝に共通の `PREP_M`（4・6・8）、時間発展は独立な `EVOLUTION_M`（8・16・32）を選びます。\n同じ模型で保存済みの9通りを利用します。任意の整数の回路を新規最適化するセルではありません。\n既定値は6と16です。変更したら、この設定セル以降を上から再実行してください。\nCASE=`exact`はMによらない対照、他のCASEは指定した部品を使います。Mと時間反復Lは別です。\n実機モードも同じ選択回路を使います。Mの異なるジョブを再利用せず、保存履歴にもMを記録します。\n演習4と5を比較するときは両冊で同じMを指定します。')
    code("PREP_M = 6         # 4, 6, 8（S0/T0共通）\nEVOLUTION_M = 16   # 8, 16, 32\nselect_circuits(PREP_M, EVOLUTION_M)")
    molecule_intro=copy.deepcopy(modeltext)
    molecule_intro.source=molecule_intro.source.split('| MODE |')[0].rstrip()
    molecule_intro.source=molecule_intro.source.replace('## 共通条件と測定方式','## 共通の分子・回路条件')
    cells.append(molecule_intro)
    md(r'''## 0．何を測り、何を推定するか

未知変数は正のT₀−S₀エネルギー差$\Delta E$（Hartree）です。
補助q₈と系q₀…q₇を全て測り、**全9ビットが0だった事象**を数えます。
時間$t=L\Delta t$、補助の補正位相$\beta=\varepsilon t$を使うと、
$$a_k(L)=\langle\psi_k|U^L|\psi_k\rangle,\quad
P(0^9|\varepsilon,L)=|a_0+e^{i\varepsilon t}a_T|^2/4.$$
厳密状態と厳密発展なら$[1+\cos((\varepsilon-\Delta E)t)]/2$です。
近似回路では高さ・位相が変わり得ます。L=1で高い信号があったことだけでは、正しいギャップは保証されません。

```text
a=q8: |0> ──┌─────┐────────P(εt)──┌─────┐──M──┐
            │  A  │               │ A†  │     ├─ 全0回数 / 全shots
s=q0..7: ───└─────┘──[U]^L─────────└─────┘──M──┘
```

本冊のFCIギャップ≈0.01590496 Haは、図・誤差の診断のためだけに使います。
測定計画の関数には正解を渡しません。演習4の9点を初期事前へ自動変換することもありません。
''')
    md(r'''## 1．事前分布と測定モードを設定する

事前は$\mathcal N(\Delta E;\mu,v)$、標準偏差は$\sqrt v$です。
初期値μ=0、v=4は著者コードに合わせた数値で、正解を中心に置いていません。
Gaussianは実数全体に定義されるため、正のギャップを対象としても負の候補を含みます。
本手順は各更新後に0で切り取る処理をしていません。

| MODE | フィットへ渡す値 | SHOTS |
|---|---|---|
| `exact` | 選んだ回路の厳密確率P | 不使用。擬似countsを作らない |
| `shots` | Binomial(R,P)で得る全0頻度n/R | 各探索点のR |
| `ibm` | 実機の全0頻度n/R | 各探索点の要求shots |

CASE=`exact`は厳密準備・厳密発展、CASE=`both`は演習4と同じ選択Mの回路です。
MODE=`exact`とCASE=`exact`は別の設定です。回路近似だけを調べるにはMODE=exactのままCASEを変えます。
MODE=shotsには実機のゲートノイズを加えていません。
MODE・CASE・SHOTS・SEEDを変えたら、**設定セルから下を再実行**してください。
''')
    md('''### 実機の認証とプレビュー

冒頭の接続準備で指定したCRN・バックエンドと認証を引き継ぎます。
MODE=ibm、CASE=bothに変更します。BACKEND_NAMEとINSTANCEは必要な場合だけ変更します。
まずRUN_HARDWARE=Falseで接続し、第2節で初回Lの費用を確認してから送信を選びます。
ジョブの再利用・保存・停止後の復旧は第9節で説明します。

[IBM Sampler入出力](https://eu-de.quantum.cloud.ibm.com/docs/en/guides/sampler-input-output)
''')
    local_config=copy.deepcopy(config)
    local_config.source=local_config.source.replace("CASE = 'both'","CASE = 'exact'")
    cells.append(local_config)
    md('''### 1.1．反復の上限と初期値

ROUNDSは最大ラウンド数です。L上限、フィット失敗、指定した目標SDによって早く止まる場合があります。
POINTSは1ラウンドの探索点数です。21点なら、shotsモードの1ラウンドは21×SHOTS回です。
実機では前ラウンドの結果を回収してから次の条件を決めます。
''')
    code('''ROUNDS = 5
INITIAL_MU = 0.0
INITIAL_VARIANCE = 4.0  # v。標準偏差はsqrt(v)
POINTS = 21
MAX_L = 20000          # この値を超える次ラウンドは実行しない
TARGET_SD = None       # 例: 0.001。GaussianのSDによる停止条件
''')
    md('''### 課題1：測定と回路の設定を分ける

1. MODE=exact、CASE=bothなら、何の誤差が残るか。
2. MODE=ibmへ変更するとき、CASE=exactを残してよいか。
3. POINTS=21、SHOTS=1024の1ラウンドは合計何shotsか。
''')
    answer('''1. shot noiseはありませんが、選択した準備・発展の回路近似は残ります。
2. 実機では選択したMの実回路を使うため、CASE=bothへ変えます。不適合なCASEは接続前に止めます。
3. 21×1024=21504 shotsです。ラウンドが増えれば合計shotsと発展の総回数も増えます。''')
    md(r'''## 2．事前(μ,v)から測定時間と探索点を選ぶ

著者コードの**数値規約**を使います。
$$t_{\rm req}=1.8/v,\quad L=\max(1,\lfloor t_{\rm req}/\Delta t\rfloor),\quad t=L\Delta t,$$
$$\varepsilon_j=\mu-v+\frac{2vj}{N-1},\quad j=0,\ldots,N-1,\qquad\beta_j=\varepsilon_j t.$$
NはPOINTSです。時間はΔtの整数倍へ切り下げます（最低1回）。
**探索幅に使うのはvで、標準偏差√vではありません。** Gaussianの分散と、探索の数値規約を区別します。
これは著者コードの時間・窓の選択で、一般のBayesian推定の普遍的な法則ではありません。

分散が減るほど時間を延ばし、窓を狭める設計です。長時間では余弦の感度が上がる一方、
周期$2\pi/t$ごとの別候補が生じます。事前で絞りながら時間を延ばす意味を意識します。
次セルは最初のラウンドの計画だけを作ります。接続・送信はしません。
実機ではRUN_HARDWARE=Falseのままこの計画を確認し、`hardware_preview(backend,L=plan1['L'])`で
**初回に実際に使うL**の費用を確認してから送信を選んでください。設定セルのL=1プレビューと、初回のL=4は区別します。
''')
    code('''plan1=choose_measurement_settings(INITIAL_MU,INITIAL_VARIANCE,POINTS)
eps1=np.asarray(plan1['trial_gaps'])
print('Initial prior mean, variance, SD:',INITIAL_MU,INITIAL_VARIANCE,np.sqrt(INITIAL_VARIANCE))
print('Requested time:',plan1['requested_time'],'L:',plan1['L'],'actual time:',plan1['time'])
print('Scan endpoints:',eps1[0],eps1[-1],'spacing:',eps1[1]-eps1[0])
print('Trial gaps / Ha:',eps1)
print('Ancilla phases / rad:',eps1*plan1['time'])
''')
    md('''### 課題2：初回の数値を読む

μ=0、v=4、Δt=0.1、POINTS=21で、要求時間、実際の時間、L、探索区間、点の間隔はいくつか。
''')
    answer('''要求時間は1.8/4=0.45、L=floor(0.45/0.1)=4、実時間は0.4です。
探索区間は[−4,4] Ha、21点なので間隔は8/20=0.4 Haです。
事前の標準偏差は2 Haですが、コードの探索幅にはv=4の数値を使っています。
初回の探索点にFCIの正解を固定していません。''')
    md(r'''## 3．選んだ点で測定する：全0countsと確率を揃える

次セルは計画した各εに対して、同じLのBPDEを測定します。
MODE=exactでは観測値yⱼ=Pⱼ、MODE=shots/ibmではyⱼ=nⱼ/Rⱼです。
$$n_j\sim\mathrm{Binomial}(R_j,P_j),\qquad \mathrm{SD}(n_j/R_j)=\sqrt{P_j(1-P_j)/R_j}.$$
全0にならなかった結果も分母に含めます。これは成功したショットへの条件付けではありません。

実機ではRUN_HARDWARE=Falseならこの段階をスキップします。
ラウンド1のジョブ名・条件は後のループと同じにするため、保存先が同じなら既存ジョブを再利用します。
''')
    code('''FIRST_READY=MEASUREMENT_READY and plan1['L']<=MAX_L and (
    MODE!='ibm' or plan1['L']<=MAX_HARDWARE_L)
FIRST_FIT_OK=False
if FIRST_READY:
    observed_1=meter.measure(CASE,plan1['L'],plan1['trial_gaps'],'round_1')
    record_1=meter.records[-1]
    print('Measurement source:',record_1['source'])
    print('Trial gap / Ha     observed y     all-zero count / total')
    for j,(eps,y) in enumerate(zip(eps1,observed_1)):
        if MODE=='exact':count_text='not sampled'
        elif MODE=='shots':count_text=f"{record_1['counts_all_zero'][j]} / {SHOTS}"
        else:
            counts=record_1['counts'][j];count_text=f"{counts.get('0'*9,0)} / {sum(counts.values())}"
        print(f'{eps:14.6f} {y:14.8f}   {count_text}')
    meter.save('outputs/exercise5_walkthrough_measurements.json')
else:
    print('First-round walkthrough skipped: preview-only mode or L limit.')
''')
    md('''### 課題3：shotsを増やす意味

1. MODE=exactのままSHOTSを増やせばフィット入力が変わるか。
2. MODE=shotsで増やすと何が減るか。回路近似による信号のずれも消えるか。
''')
    answer('''1. 変わりません。厳密確率を直接使い、SHOTSも擬似countsも使いません。
2. 測定頻度の統計変動は一般に1/√SHOTSで減ります。回路の系統的なずれは消えません。
MODE=shotsは理想回路の有限shotsだけを再現するモードです。''')
    md(r'''## 4．測定した確率曲線をGaussianへフィットする

観測値yⱼへ、次の曲線を非線形最小二乗で当てはめます。
$$g(\varepsilon)=A\exp[-(\varepsilon-m)^2/(2w)],\qquad w=s^2.$$
`curve_fit`のパラメータは(A,m,s)、初期値は(1,μ,√v)です。Aは高さ、mは中心、sは幅です。
この節のA,m,sはフィットのパラメータで、回路A、QPEの補助数m、Trotter分割数sとは別の量です。
コードの`amplitude`、`fit_mu`、`fit_sd`に対応します。
元の信号は余弦なので、Gaussianは近似です。実機では信号消失やノイズで形が変わることもあります。

**wはフィットした曲線の幅の二乗です。** `curve_fit`が返すパラメータの共分散や、中心mの推定誤差ではありません。
このwを次の更新に使うのが本手順です。フィットのAがわずかに1を超えても、それは物理確率を生成し直す値ではなく近似曲線のパラメータです。

次セルで入力・フィット値・残差RMSEを確認します。
RMSE=√mean[(yⱼ−g(εⱼ))²]は曲線との不一致の診断で、厳密な統計的合否判定ではありません。
''')
    code('''if FIRST_READY:
    try:
        fit1=fit_gaussian_signal(plan1['trial_gaps'],observed_1,INITIAL_MU,INITIAL_VARIANCE)
        FIRST_FIT_OK=True
        print('Fit amplitude A:',fit1['amplitude'])
        print('Fit center m / Ha:',fit1['fit_mu'],'width s / Ha:',fit1['fit_sd'],'w:',fit1['fit_variance'])
        print('Residual RMSE:',fit1['residual_rmse'],'warnings:',fit1['warnings'])
    except (ValueError,RuntimeError,FloatingPointError) as exc:
        print('Fit failed; measured data retained, no update:',exc)
''')
    md('''### 課題4：Gaussian幅とフィットの良否

1. wは「フィット中心の標準誤差の二乗」か。
2. shot noiseがないとRMSEは必ず0になるか。
3. フィットが数値的に終了すれば、信号の形を十分に表しているか。
''')
    answer('''1. 違います。wは確率曲線をGaussianで近似したときの幅²です。パラメータ共分散を更新へ使う処理ではありません。
2. なりません。余弦をGaussianで近似することによる残差が、理想データにも残ります。
3. 数値的終了だけでは保証できません。図Bの形と図Cの残差も確認します。幅が非常に大きい平坦なフィットなどは特に注意して読みます。''')
    md(r'''## 5．事前とフィットGaussianを掛け、事後へ更新する

フィット曲線をGaussian因子として使い、事前と掛け合わせて正規化します。
$$q_{\rm post}(\Delta E)\propto\mathcal N(\Delta E;\mu,v)
\exp[-(\Delta E-m)^2/(2w)].$$
Aは正規化で消えます。平方完成すると
$$v'=(1/v+1/w)^{-1}=\frac{vw}{v+w},\qquad
\mu'=v'(\mu/v+m/w)=\frac{vm+w\mu}{v+w}.$$
次セルの`update_gaussian_prior`がこの式を計算します。

このGaussian因子は、全観測countsの二項尤度を直接計算したものではありません。
演習1・2の厳密な二項モデルの更新と、今回の**Gaussian近似による更新**を区別してください。
正のv,wを使う限りv′<vとなるため、幅が縮むこと自体は正解へ近づいた証拠になりません。
''')
    code('''if FIRST_FIT_OK:
    mu1,v1=update_gaussian_prior(INITIAL_MU,INITIAL_VARIANCE,fit1['fit_mu'],fit1['fit_variance'])
    print('Prior mean, variance:',INITIAL_MU,INITIAL_VARIANCE)
    print('Posterior mean, variance, SD:',mu1,v1,np.sqrt(v1))
    # 同じ式を直接書いて関数の出力と照合。
    w=fit1['fit_variance'];m=fit1['fit_mu'];v=INITIAL_VARIANCE;mu=INITIAL_MU
    np.testing.assert_allclose([mu1,v1],[(v*m+w*mu)/(v+w),v*w/(v+w)],atol=1e-12)
''')
    md(r'''## 6．事後を次の事前へ、そのまま渡す

更新後は$\mu_{r+1}=\mu'_r$、$v_{r+1}=v'_r$です。
次の時間・探索点を決める際にFCIの正解で中心を補正したり、初期分布へ戻したりしません。
次セルは更新後の分布からラウンド2の計画だけを作ります。
''')
    code('''if FIRST_FIT_OK:
    plan2=choose_measurement_settings(mu1,v1,POINTS)
    print('Next prior:',mu1,v1)
    print('Next L,time:',plan2['L'],plan2['time'])
    print('Next scan window:',plan2['trial_gaps'][0],plan2['trial_gaps'][-1])
''')
    md('''### 課題5：引継ぎと停止条件

1. 初回のフィット中心mを、そのまま次の事前平均へ渡すのか。
2. vが小さくなると、次の時間と探索窓はどう変わるか。
3. 小さな目標SDへ達したら、実機や回路近似を含む全誤差が小さいと断定できるか。
''')
    answer('''1. 渡すのはm単独ではなく、事前とフィットを掛けた後のμ′です。分散もv′を渡します。
2. 数値規約t_req=1.8/vにより時間が長くなり、探索窓μ±vが狭くなります。整数化でLが同じままのこともあります。
3. 断定できません。目標SDはこのGaussian更新の内部の幅です。回路近似、測定モデル、実機の影響を含む誤差保証ではありません。''')
    md(r'''## 7．同じ処理を反復し、各ラウンドの図と照合する

ここからは`iterative_bayes`が第2〜6節をループで実行します。
初期分布からやり直し、meterも同じ設定・SEEDで作り直すので、模擬測定の初回は先ほどの段階実行と一致します。
実機の初回は同条件の保存ジョブを再利用します。新しいshotsを追加する処理ではありません。

各ラウンドの図は次の対応です。

- **A**：事前の密度、μ±vの探索窓、実際に選んだ点、Lと実時間。
- **B**：観測値とGaussianフィット。破線は選んだ回路の理想信号で、実機データとは区別します。
- **C**：観測−フィットの残差。曲線の形の不一致やフィット異常を確認します。
- **D**：事前、面積1に正規化したフィットGaussian、事後。全て同じギャップ軸です。

FCI線は診断用です。図の線が測定計画や更新に入力されることはありません。
実機は次のLがMAX_HARDWARE_Lを超える前に停止し、それまでの履歴を表示・保存します。
''')
    code('''if MEASUREMENT_READY:
    meter=Measurement(MODE,SHOTS,SEED,service,backend,RUN_HARDWARE,MAX_HARDWARE_L)
    history=iterative_bayes(meter,case=CASE,rounds=ROUNDS,initial_mu=INITIAL_MU,
        initial_variance=INITIAL_VARIANCE,points=POINTS,max_L=MAX_L,target_sd=TARGET_SD)
else:
    history=[]
    print('Preview only; no measurement was submitted.')
''')
    md(r'''### 7.1．数値履歴と、幅・誤差・時間の推移

次セルの表は各図と同じ数値です。前行の事後(μ′,v′)と、次行の事前(μ,v)の一致をassertで確認します。
mean±SDは今回のGaussianの平均と幅です。95%信用区間や全誤差保証としては表示していません。
絶対誤差|mean−FCI|とSDを同じ図に描き、「幅が減る」と「真値に近づく」を分けます。
Lと、選択回路の概算CX費用 CX(A)+CX(A†)+L×CX(Uapp)も確認します（既定6/16は516+216L）。実機配置・全体最適化前のブロック直列見積もりです。
''')
    code('''if history:
    print('r  prior mu    prior v     L     t       fit m        fit w       post mu      post v      RMSE')
    for row in history:
        print(f"{row['round']:1d} {row['prior_mu']:10.6f} {row['prior_variance']:9.5g} {row['L']:5d} "
              f"{row['time']:7.2f} {row['fit_mu']:11.7f} {row['fit_variance']:11.5g} "
              f"{row['posterior_mu']:11.7f} {row['posterior_variance']:11.5g} {row['residual_rmse']:.4g}")
    for a,b in zip(history,history[1:]):
        assert a['posterior_mu']==b['prior_mu'] and a['posterior_variance']==b['prior_variance']
    rounds=np.arange(len(history)+1)
    means=np.array([INITIAL_MU]+[row['posterior_mu'] for row in history])
    sd=np.sqrt([INITIAL_VARIANCE]+[row['posterior_variance'] for row in history])
    fig,axes=plt.subplots(1,3,figsize=(13,3))
    axes[0].errorbar(rounds,means,yerr=sd,fmt='o-');axes[0].axhline(GAP,color='black',ls='--')
    axes[0].set(xlabel='Round',ylabel='Mean +/- Gaussian SD / Ha')
    axes[1].semilogy(rounds,sd,'o-',label='Gaussian SD')
    axes[1].semilogy(rounds,np.maximum(abs(means-GAP),1e-15),'s-',label='Absolute error vs FCI')
    axes[1].set(xlabel='Round',ylabel='Ha');axes[1].legend(fontsize=8)
    axes[2].plot(rounds[1:],[row['L'] for row in history],'o-')
    axes[2].set(xlabel='Round',ylabel='Evolution repetitions L')
    plt.tight_layout();plt.show();plt.close(fig)
    print('Posterior -> next prior checked. Source:',MODE,'CASE:',CASE)
    print(f'M{SELECTED_PREP_M}/M{SELECTED_EVOLUTION_M} serial CX estimates:',
          [SERIAL_CX['preparation']+SERIAL_CX['evolution']*row['L'] for row in history])
else:
    print('No updates. Check preview-only mode, fit failure or iteration limits above.')
''')
    md('''### 課題6：ラウンドの図と表を読む

1. 初期事前を固定し、ROUNDSを2→5→8へ増やすと何が変わるか。
2. MODE=shotsでSHOTSやSEEDを変えると、Lの列も変わり得るか。
3. 図Dの事後が狭くても、図Cの残差が大きい場合はどう解釈すべきか。
''')
    answer('''1. 更新を追加します。vが減れば要求Lが大きくなり、L上限に達して指定回数の前に止まる場合があります。
先行ラウンドを同じ設定で再実行した部分は、同じSEEDなら同じ結果です。
2. 測定値が変わるとフィットと事後が変わり、次のLや探索点も変わり得ます。適応的な計画だからです。
3. Gaussianの幅縮小だけを成功とみなせません。曲線モデルやフィットがデータを十分に表しているか、観測の出所・信号・初期条件を確認します。''')
    md(r'''## 8．shot noiseを外して、回路近似によるずれを比較する

次の対照は**常にMODE=exactのシミュレーション**です。主実験がibmでも、このセルから実機へ送信しません。
演習4と同じ4条件を同じ初期事前・更新規約で反復し、準備と発展の影響を切り分けます。
主実験のmeter・history・保存ファイルは書き換えません。

CASE=exactでも、有限回数・Gaussian近似・フィット不調による誤差は残り得ます。
CASE=bothでは、元の固有状態でない準備や別の近似ユニタリの位相も影響します。
近似が相殺する場合もあり、「両方近似」が常に「片側だけ近似」より悪いとは限りません。

COMPARE_ROUNDS=8までを試しますが、次のLが上限を超えたらそこまでです。
各CASEで独立に適応するため、途中からL・探索点・実行ラウンド数が異なり得ます。
これは同じ初期条件からの比較で、同じ総発展時間や測定予算を保証する比較ではありません。
固定Lでの純粋な信号比較は演習4を参照してください。
''')
    code('''RUN_COMPARISON=True
COMPARE_ROUNDS=8
comparison={};comparison_meters={}
if RUN_COMPARISON:
    for case in ['exact','prep_only','evolution_only','both']:
        control=Measurement('exact',shots=None,seed=SEED)
        rows=iterative_bayes(control,case=case,rounds=COMPARE_ROUNDS,initial_mu=INITIAL_MU,
            initial_variance=INITIAL_VARIANCE,points=POINTS,max_L=MAX_L,
            show_plots=False,output_prefix=None,verbose=False)
        comparison[case]=rows;comparison_meters[case]=control
    print('No-shot controls: case / completed rounds / mean / SD / signed error / last residual RMSE')
    fig,axes=plt.subplots(1,3,figsize=(14,3.5))
    for i,(case,rows) in enumerate(comparison.items()):
        if not rows:continue
        final=rows[-1];sd=np.sqrt(final['posterior_variance'])
        print(case,len(rows),final['posterior_mu'],sd,final['posterior_mu']-GAP,final['residual_rmse'])
        axes[0].errorbar(i,final['posterior_mu'],yerr=sd,fmt='o')
        rr=[row['round'] for row in rows]
        axes[1].semilogy(rr,[max(abs(row['posterior_mu']-GAP),1e-15) for row in rows],'o-',label=case)
        axes[2].semilogy(rr,[max(row['residual_rmse'],1e-15) for row in rows],'o-',label=case)
    axes[0].set_xticks(range(4),list(comparison),rotation=25);axes[0].axhline(GAP,color='black',ls='--')
    axes[0].set(ylabel='Final mean +/- Gaussian SD / Ha',title='Last completed round per case')
    axes[1].set(xlabel='Round',ylabel='Absolute error vs FCI / Ha');axes[1].legend(fontsize=7)
    axes[2].set(xlabel='Round',ylabel='Gaussian-fit residual RMSE');axes[2].legend(fontsize=7)
    plt.tight_layout();plt.show();plt.close(fig)
''')
    md('''### 課題7：偏りとGaussian更新の限界

1. MODE=exactで残るずれをshot noiseと説明してよいか。
2. 一つのCASEの残差が急に増えたラウンドで、平均だけを見て収束を主張してよいか。
3. CASE=bothの平均が偶然FCIに近ければ、準備と発展がそれぞれ正確だったといえるか。
''')
    answer('''1. shot noiseは加えていません。回路近似、有限反復、Gaussian近似、フィットの不調を検討します。
特にevolution_onlyの比較は、厳密準備を保ったまま近似発展の位相の影響を調べる対照です。
2. 主張できません。観測曲線とフィットを見て、残差の増加が形の不一致や平坦な局所解などによらないか確認します。
3. いえません。二つの影響が相殺する場合があります。prep_only/evolution_only、演習4の各fidelityと信号も合わせて読みます。''')
    md(r'''### 8.1．Mごとの反復推定を比較する

準備M=4・6・8 × 発展M=8・16・32を、**CASE=both・MODE=exact**で比較します。
主実験のMODEがshots/ibmでも、この比較は常にshot noiseなしで、実機へ送信しません。
同じINITIAL_MU・INITIAL_VARIANCE・POINTS・MAX_Lと更新回数から開始します。
各Mの測定結果から独立に次のL・探索点を選ぶため、途中から測定計画や完了回数が異なり得ます。
M_COMPARE_ROUNDSを変更すると比較の更新上限を変えられます。

表の「平均−FCI」は、この有限反復の結果の符号付き誤差です。
ここで比較する偏りには回路近似だけでなくGaussian近似・有限反復・フィットの影響も含まれます。
厳密準備・厳密発展の対照も同じ初期条件で計算し、図の破線で示します。
SDはGaussian更新の幅で、全誤差を覆う保証ではありません。L=1のピークずれとも区別します。
比較後は元のMへ戻し、主実験のmeter・history・保存ファイルを保持します。
''')
    code(r'''from IPython.display import Markdown, display
M_COMPARE_ROUNDS=8

def compare_saved_m_bayes():
    selected=(SELECTED_PREP_M,SELECTED_EVOLUTION_M)
    results={}
    kwargs=dict(rounds=M_COMPARE_ROUNDS,initial_mu=INITIAL_MU,initial_variance=INITIAL_VARIANCE,
                points=POINTS,max_L=MAX_L,show_plots=False,output_prefix=None,verbose=False)
    reference=iterative_bayes(Measurement('exact',shots=None),case='exact',**kwargs)
    try:
        for pm in AVAILABLE_PREP_M:
            for em in AVAILABLE_EVOLUTION_M:
                select_circuits(pm,em)
                results[(pm,em)]=iterative_bayes(Measurement('exact',shots=None),case='both',**kwargs)
    finally:
        select_circuits(*selected)
    return results,reference

M_BAYES_COMPARISON,M_BAYES_REFERENCE=compare_saved_m_bayes()
header='| 準備M | 発展M | 完了round | 最終L | 最終平均 / Ha | Gaussian SD / Ha | 平均−FCI / Ha | 最終fit RMSE |\n|---|---|---|---|---|---|---|---|'
lines=[]
for (pm,em),rows in M_BAYES_COMPARISON.items():
    if not rows:
        lines.append(f'| {pm} | {em} | 0 | — | — | — | — | — |')
        continue
    r=rows[-1]
    lines.append(f"| {pm} | {em} | {len(rows)} | {r['L']} | {r['posterior_mu']:.8f} | {np.sqrt(r['posterior_variance']):.3g} | {r['posterior_mu']-GAP:+.8f} | {r['residual_rmse']:.3g} |")
display(Markdown(header+'\n'+'\n'.join(lines)))
fig,axes=plt.subplots(1,3,figsize=(13,3.5),sharey=True)
for ax,pm in zip(axes,AVAILABLE_PREP_M):
    for em in AVAILABLE_EVOLUTION_M:
        rows=M_BAYES_COMPARISON[(pm,em)]
        ax.plot([r['round'] for r in rows],[r['posterior_mu']-GAP for r in rows],
                'o-',label=f'evolution M={em}')
    ax.plot([r['round'] for r in M_BAYES_REFERENCE],
            [r['posterior_mu']-GAP for r in M_BAYES_REFERENCE],'k--',label='exact prep/evolution')
    ax.axhline(0,color='gray',lw=.7)
    ax.set(title=f'Preparation M={pm}',xlabel='Round',ylabel='Posterior mean - FCI / Ha')
    ax.legend(fontsize=7)
plt.tight_layout();plt.show();plt.close(fig)
''')
    md('''### 課題7補足：Mを増やした推定を評価する

最終SDが最小のMを、最も偏りが小さいMとみなしてよいか。
また、同じroundの各点は同じL・同じ測定予算での比較か。
''')
    answer('''どちらもいえません。SDは近似Gaussian更新の幅であり、FCIからのずれは表と図の符号付き誤差で別に比較します。
厳密回路の破線にも有限反復やGaussianフィット由来のずれが残り得ます。
各Mが測定結果に応じてLと探索点を選ぶため、同じroundでもL・回路費用は揃いません。
上限によって完了round数が異なる場合は最終値の単純な順位付けも避け、共通のroundと各L・残差を確認します。
固定LでMの信号を比較する対照は演習4の表と図です。''')
    md(r'''## 9．実機で反復するときの操作と、保存結果の読み方

1ラウンドごとにPOINTS本を1ジョブで送り、countsを取得・フィットしてから次の条件を決めます。
実機の待ち時間やノイズがあるため、シミュレーションと同じL列になるとは限りません。

- RUN_HARDWARE=Falseはプレビューのみです。段階実行も反復測定もスキップします。
- 送信直後に`outputs/ibm_jobs/`へジョブIDを保存します。同じ条件なら再投入せず回収します。
- 手順を最初からやり直すときも、保存先と条件が同じなら初回から再利用し、計画を再計算します。
- 新しい独立測定が必要なら、元の記録を保持したまま`Measurement(directory=...)`へ別の保存先を渡します。
- 主実験の履歴は`outputs/exercise5_updates.json`、測定値は`outputs/exercise5_measurements.json`です。
- MAX_HARDWARE_Lを超えるラウンドは送信前に止まり、そこまでの事後を残します。次のLは`hardware_preview`で費用を確認します。

フィット失敗時は観測を保存し、不適切な事後を追加せずに終了します。
残差が大きい場合は自動的に正解へ補正しません。SDだけでなく曲線とデータの出所を確認してください。
Colabでは切断に備え、ジョブ記録をDrive等へ保存します。記録を失うと自動再利用ができません。
本教材の検証ではIBM実機へ送信していません。

### 課題8：反復実験の結果を説明する

事前・測定条件・counts/確率・フィット・事後・次の計画を、保存履歴と図で対応させて説明してください。
同時に、表示SDを何の保証として扱ってはいけないか述べてください。
''')
    cells[-1].source+='''

### 解答例（課題8：全体の振り返り）

事前のμ,vからt_req・L・探索点を選び、全0確率または全shotsに対する頻度を得ます。
その曲線のGaussianフィットから中心mと幅²wを求め、事前との積でμ′,v′を更新し、次の事前へ渡します。
図Aが計画、B/Cが測定とフィット、Dが積更新です。表では前行の事後と次行の事前が一致します。
表示SDはこのGaussian近似更新の幅であり、二項尤度の厳密な事後幅や、回路近似・実機ノイズを含む全誤差の保証ではありません。
'''
    cells.append(copy.deepcopy(checks))
    md('''### 出所と再現範囲

[Kanno et al., PNAS (2025)](https://doi.org/10.1073/pnas.2425026122)と、
[固定コミットの著者コード](https://github.com/sk888ks/TQPDE_open/blob/9cd88464c5c50dda783150661813609092d7f2b6/src/QPDE_Bayes_session_MPSBWMPO_RD.py)のGaussianフィット・更新を参照します。
本冊はH₄と選択した保存回路を用いた説明用の計算で、論文の分子・実機実験全体の再現ではありません。
''')
    return cells
