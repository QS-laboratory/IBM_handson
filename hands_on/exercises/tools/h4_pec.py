"""Simulation-only H4 PEC helpers embedded in 04; no IBM dependencies."""
import hashlib
import io
import json
import os
from pathlib import Path
from urllib.request import urlopen
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from scipy.interpolate import PchipInterpolator
from qiskit import qpy
from qiskit.quantum_info import Statevector, Operator
from IPython.display import Markdown, display

DT = 0.1
DATA_COMMIT = '4ab6ef3d06f4256b5c818feb8c3fa2af4381dd43'
DATA_URL = f'https://raw.githubusercontent.com/QS-laboratory/IBM_handson/{DATA_COMMIT}/hands_on/exercises/data/h4_pec_v1/'
_MODEL_CACHE, _STATE_CACHE, _SPECTRUM_CACHE = {}, {}, {}


def data_bytes(relative):
    if relative not in H4_INDEX['files']:
        raise ValueError('Unknown snapshot file: ' + relative)
    roots = []
    if os.environ.get('H4_DATA_DIR'):
        roots.append(Path(os.environ['H4_DATA_DIR']).expanduser())
    for parent in [Path.cwd(), *Path.cwd().parents]:
        roots.extend([parent/'data/h4_pec_v1', parent/'hands_on/exercises/data/h4_pec_v1'])
    cache = Path('outputs/h4_pec_v1')
    roots.append(cache)
    for root in roots:
        path = root/relative
        if path.is_file():
            raw = path.read_bytes()
            break
    else:
        with urlopen(DATA_URL + relative, timeout=90) as response:
            raw = response.read()
        if hashlib.sha256(raw).hexdigest() != H4_INDEX['files'][relative]:
            raise ValueError('Downloaded snapshot checksum mismatch: ' + relative)
        path = cache/relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(raw)
    if hashlib.sha256(raw).hexdigest() != H4_INDEX['files'][relative]:
        raise ValueError('Snapshot checksum mismatch: ' + relative)
    return raw


def saved_model(geometry):
    if geometry not in _MODEL_CACHE:
        with np.load(io.BytesIO(data_bytes(geometry+'/model.npz')), allow_pickle=False) as data:
            _MODEL_CACHE[geometry] = {k:data[k].copy() for k in data.files}
    return _MODEL_CACHE[geometry]


def check_choice(value, choices, name):
    if value == 'exact':
        return
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value not in choices:
        raise ValueError(f'{name}: exact または {choices} を指定してください')


def saved_states(geometry, prep_m):
    check_choice(prep_m, [4,6,8], '準備M')
    key = geometry, prep_m
    if key not in _STATE_CACHE:
        model = saved_model(geometry)
        if prep_m == 'exact':
            states = [model['S0'].copy(), model['T'].copy()]
        else:
            states = []
            for state in ['S0','T']:
                qc = qpy.load(io.BytesIO(data_bytes(f'{geometry}/prep_{prep_m}_{state}.qpy')))[0]
                states.append(Statevector.from_instruction(qc).data)
        for psi in states:
            np.testing.assert_allclose(np.vdot(psi,psi), 1, atol=1e-10)
        _STATE_CACHE[key] = states
    return _STATE_CACHE[key]


def evolution_spectrum(geometry, evolution_m):
    check_choice(evolution_m, [8,16,32], '発展M')
    key = geometry, evolution_m
    if key not in _SPECTRUM_CACHE:
        if evolution_m == 'exact':
            energy, vectors = eigh(saved_model(geometry)['hamiltonian'])
            phases = -energy*DT
        else:
            qc = qpy.load(io.BytesIO(data_bytes(f'{geometry}/evolution_{evolution_m}.qpy')))[0]
            u = Operator(qc).data
            # All saved one-step phases lie in (-pi/2, pi/2): sin(theta) is one-to-one.
            assert np.min(np.linalg.eigvalsh((u+u.conj().T)/2)) > 0
            _, vectors = eigh((u-u.conj().T)/(2j))
            phases = np.angle(np.diag(vectors.conj().T@u@vectors))
            np.testing.assert_allclose((vectors*np.exp(1j*phases))@vectors.conj().T, u, atol=1e-9)
        _SPECTRUM_CACHE[key] = phases, vectors
    return _SPECTRUM_CACHE[key]


def expectation_pec(prep_m=8):
    check_choice(prep_m, [4,6,8], '準備M')
    rows = []
    for point in sorted(H4_INDEX['geometries'], key=lambda row:row['R']):
        model = saved_model(point['id'])
        psi = saved_states(point['id'], prep_m)
        electronic = [float(np.vdot(s, model['hamiltonian']@s).real) for s in psi]
        enuc = float(model['enuc'])
        rows.append(dict(geometry=point['id'], R=point['R'], prep_m=prep_m,
            S0_expect=electronic[0]+enuc, T_expect=electronic[1]+enuc, gap_expect=electronic[1]-electronic[0],
            S0_FCI=float(model['energy_S0'])+enuc, T_FCI=float(model['energy_T'])+enuc,
            gap_FCI=float(model['energy_T']-model['energy_S0'])))
    return rows


def continuous_phase(phases, weights, lmax):
    # Exact H uses its physical eigenphases. Saved U uses the principal-log interpolation.
    derivative = np.abs(phases)@weights
    for subdivisions in [1,2,4,8,16,32,64,128,256,512,1024]:
        angles = []
        minimum = np.ones(2)
        previous = np.ones(2,complex)
        accumulated = np.zeros(2)
        coarse_previous = np.ones(2,complex)
        coarse_phase = np.zeros(2)
        captured = []
        for start in range(0,lmax*subdivisions+1,2048):
            ix = np.arange(start,min(start+2048,lmax*subdivisions+1))
            a = np.exp(1j*(ix/subdivisions)[:,None]*phases)@weights
            increments = np.angle(a*np.conj(np.vstack([previous,a[:-1]])))
            unwrapped = accumulated + np.cumsum(increments,axis=0)
            accumulated = unwrapped[-1];previous = a[-1]
            minimum = np.minimum(minimum,abs(a).min(axis=0))
            chosen = ix%subdivisions == 0
            angles.extend(unwrapped[chosen]);captured.extend(a[chosen])
            if subdivisions>1:
                coarse = a[ix%2==0]
                if len(coarse):
                    coarse_phase += np.angle(coarse*np.conj(np.vstack([coarse_previous,coarse[:-1]]))).sum(axis=0)
                    coarse_previous = coarse[-1]
        bound = float(np.max(derivative/(subdivisions*np.maximum(minimum,1e-300))))
        agreement = float(np.max(abs(accumulated-coarse_phase))) if subdivisions>1 else 0.0
        if bound<0.2 and agreement<1e-8:
            return np.array(captured), np.array(angles), dict(subdivisions=subdivisions,
                minimum_amplitude=minimum.tolist(),derivative_ratio=bound,phase_check_error=agreement)
    raise ValueError('位相の枝を確認できません。L_MAXを小さくしてください。')


def phase_pec_scan(prep_m=8, evolution_m=32, lmax=64):
    check_choice(prep_m, [4,6,8], '準備M');check_choice(evolution_m, [8,16,32], '発展M')
    if isinstance(lmax,(bool,np.bool_)) or not isinstance(lmax,(int,np.integer)) or not 1<=lmax<=1024:
        raise ValueError('L_MAXは1〜1024の整数です')
    beta = np.arange(9)*2*np.pi/9
    design = np.c_[np.ones(9), np.cos(beta), np.sin(beta)]
    inverse = np.linalg.pinv(design)
    rows = []
    expectation = expectation_pec(prep_m)
    for point in expectation:
        geometry = point['geometry']
        phases,vectors = evolution_spectrum(geometry,evolution_m)
        states = saved_states(geometry,prep_m)
        weights = np.stack([abs(vectors.conj().T@s)**2 for s in states],axis=1)
        np.testing.assert_allclose(weights.sum(axis=0),1,atol=1e-9)
        amplitudes, phase, certificate = continuous_phase(phases,weights,lmax)
        for L in range(1,lmax+1):
            a0,at = amplitudes[L]
            bpe = (1+np.real(np.exp(1j*beta)*a0))/2
            bpde = abs(a0+np.exp(1j*beta)*at)**2/4
            coefficients = np.stack([inverse@bpe,inverse@bpde])
            if np.min(np.hypot(coefficients[:,1],coefficients[:,2]))<1e-10:
                raise ValueError(f'信号が平坦で位相を読めません: {geometry}, L={L}')
            principal = np.arctan2(coefficients[:,2],coefficients[:,1])
            anchors = np.array([-phase[L,0],phase[L,0]-phase[L,1]])
            lifted = principal+2*np.pi*np.rint((anchors-principal)/(2*np.pi))
            np.testing.assert_allclose(lifted,anchors,atol=1e-8)
            fit_error = float(np.max(abs(design@coefficients.T-np.c_[bpe,bpde])))
            energy0,gap = lifted/(L*DT)
            enuc = float(saved_model(geometry)['enuc'])
            rows.append(dict(point, evolution_m=evolution_m,L=L,t=L*DT,
                S0_phase=float(energy0+enuc),T_phase=float(energy0+gap+enuc),gap_phase=float(gap),
                max_fit_residual=fit_error,phase_certificate=certificate))
        print(f'完了 {geometry}: R={point["R"]:.5f} Å, L=1..{lmax}',flush=True)
    return rows


def save_rows(rows, filename):
    path = Path(filename);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(rows,ensure_ascii=False,indent=2)+'\n')


def plot_pec(rows, phase=False, show_l=None):
    points = rows if not phase else [r for r in rows if r['L']==1]
    radii = np.array([r['R'] for r in points]);fine = np.linspace(min(radii),max(radii),400)
    fig,axes = plt.subplots(2,3,figsize=(14,7))
    for j,(key,label) in enumerate([('S0','S0 total energy'),('T','T0 total energy'),('gap','T0 - S0')]):
        reference = np.array([r[key+'_FCI'] for r in points])
        baseline = np.array([r[key+'_expect'] for r in points])
        axes[0,j].plot(fine,PchipInterpolator(radii,reference)(fine),'k-',label='FCI')
        axes[0,j].plot(radii,baseline,'s--',ms=3,label='<H>')
        axes[1,j].plot(radii,1000*(baseline-reference),'s--',ms=3,label='<H> - FCI')
        if phase:
            for L in show_l:
                selected = [r for r in rows if r['L']==L]
                estimated = np.array([r[key+'_phase'] for r in selected])
                line, = axes[0,j].plot(fine,PchipInterpolator(radii,estimated)(fine),label=f'L={L}')
                axes[0,j].plot(radii,estimated,'o',ms=2,color=line.get_color())
                axes[1,j].plot(radii,1000*(estimated-reference),'o-',ms=2,label=f'L={L}')
        axes[0,j].set(title=label,ylabel='Energy / Ha')
        axes[1,j].axhline(0,color='black',lw=.6);axes[1,j].set(ylabel='Signed error / mHa')
        for ax in axes[:,j]:ax.set_xlabel('Adjacent H-H distance / Angstrom');ax.legend(fontsize=8)
    plt.tight_layout();plt.show();plt.close(fig)


def print_pec(rows, kind='expect'):
    text=['| R / Å | S₀ / Ha | T₀ / Ha | T₀−S₀ / Ha |','|---:|---:|---:|---:|']
    text += [f"| {r['R']:.8f} | {r['S0_'+kind]:.10f} | {r['T_'+kind]:.10f} | {r['gap_'+kind]:.10f} |" for r in rows]
    display(Markdown('\n'.join(text)))


def plot_l_errors(rows):
    ls = sorted({r['L'] for r in rows})
    fig,axes = plt.subplots(1,3,figsize=(13,3.3))
    summary=[]
    for key,ax in zip(['S0','T','gap'],axes):
        errors = [1000*np.sqrt(np.mean([(r[key+'_phase']-r[key+'_FCI'])**2 for r in rows if r['L']==L])) for L in ls]
        baseline = 1000*np.sqrt(np.mean([(r[key+'_expect']-r[key+'_FCI'])**2 for r in rows if r['L']==1]))
        ax.plot(ls,errors,label='BPE / BPDE');ax.axhline(baseline,ls='--',label='<H>')
        ax.set(xlabel='L',ylabel='PEC RMSE / mHa',title=key);ax.legend()
        summary.append((key,baseline,errors[-1]))
    plt.tight_layout();plt.show();plt.close(fig)
    print('Quantity   <H> RMSE / mHa   final-L RMSE / mHa')
    for key,a,b in summary:print(f'{key:8s} {a:17.9g} {b:21.9g}')
