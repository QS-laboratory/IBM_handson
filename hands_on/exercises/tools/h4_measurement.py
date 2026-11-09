"""Embedded runtime for exercises 4/5. ASSETS is supplied by the builder."""
import base64
import hashlib
import io
import json
import os
from pathlib import Path
import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm, eigh
from scipy.optimize import curve_fit, OptimizeWarning
from qiskit import QuantumCircuit, qpy, transpile
from qiskit.quantum_info import Operator, Statevector, SparsePauliOp
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.synthesis import SuzukiTrotter

DT = .1

def load_array_asset(name):
    return np.load(io.BytesIO(base64.b64decode(ASSETS[name])), allow_pickle=False)

def load_circuit_asset(name):
    key = globals().get('ACTIVE_ASSET_KEYS', {}).get(name, name)
    return qpy.load(io.BytesIO(base64.b64decode(ASSETS[key])))[0]

with load_array_asset('model') as saved:
    MODEL = {k: saved[k].copy() for k in saved.files}
H = MODEL['hamiltonian']
GAP = float(MODEL['energy_T'] - MODEL['energy_S0'])
EXACT_STATES = [MODEL['S0'], MODEL['T']]
U_EXACT = expm(-1j * H * DT)
A_EXACT = prepare_pair(*EXACT_STATES)
SELECTED_GEOMETRY = 'R_2p00000000'
AVAILABLE_PREP_M = (4,6,8)
AVAILABLE_EVOLUTION_M = (8,16,32)


def select_circuits(prep_m=6, evolution_m=16):
    """Select compatible saved circuits. Rerun following cells after changing M."""
    global SELECTED_PREP_M, SELECTED_EVOLUTION_M, ACTIVE_ASSET_KEYS
    global PREP_INFO, EVOLUTION_INFO, SELECTED_RESULTS, CIRCUIT_DIGEST, SERIAL_CX
    global A_APPROX, EVOLUTION_APPROX, U_APPROX, APPROX_STATES, CASES, SPECTRA
    for value, choices, name in [(prep_m, AVAILABLE_PREP_M, 'PREP_M'),
                                  (evolution_m, AVAILABLE_EVOLUTION_M, 'EVOLUTION_M')]:
        if value != 'exact' and (not isinstance(value, (int, np.integer)) or isinstance(value, (bool, np.bool_))
                or value not in choices):
            raise ValueError(f'{name} must be exact or one of {choices}; received {value!r}')
    prep = CIRCUIT_CATALOG['preparation'][prep_m]
    evo = CIRCUIT_CATALOG['evolution'][evolution_m]
    pair = CIRCUIT_CATALOG['pairs'][f'{prep_m}/{evolution_m}']
    keys = dict(prep['keys'], evolution=evo['key'])
    circuits = {k:qpy.load(io.BytesIO(base64.b64decode(ASSETS[v])))[0] for k,v in keys.items()}
    states = [Statevector.from_instruction(circuits[k]).data for k in ['S0','T']]
    uapp = Operator(circuits['evolution']).data
    cases = {'exact':(EXACT_STATES,U_EXACT), 'prep_only':(states,U_EXACT),
             'evolution_only':(EXACT_STATES,uapp), 'both':(states,uapp)}
    # Hermitian phase decomposition accelerates U^L classically; circuits repeat L times.
    spectra = {}
    for label, (_, u) in cases.items():
        assert np.linalg.norm(u.conj().T@u-np.eye(256)) < 1e-9
        assert np.min(np.linalg.eigvalsh((u+u.conj().T)/2)) > 0
        _, vectors = eigh((u-u.conj().T)/(2j))
        phases = np.angle(np.diag(vectors.conj().T@u@vectors))
        np.testing.assert_allclose((vectors*np.exp(1j*phases))@vectors.conj().T,u,atol=1e-10)
        spectra[label] = (phases, vectors)
    SELECTED_PREP_M, SELECTED_EVOLUTION_M = prep_m, evolution_m
    ACTIVE_ASSET_KEYS = keys
    PREP_INFO, EVOLUTION_INFO = prep, evo
    SELECTED_RESULTS = dict(fidelities=prep['fidelities'], fringe=pair['fringe'])
    PROVENANCE['results'] = SELECTED_RESULTS  # legacy teaching diagnostics
    CIRCUIT_DIGEST = pair['digest']
    SERIAL_CX = dict(preparation=prep['resources']['A']['cx']+prep['resources']['inverse']['cx'],
                     evolution=evo['resources']['cx'])
    A_APPROX, EVOLUTION_APPROX = circuits['A'], circuits['evolution']
    APPROX_STATES, U_APPROX, CASES, SPECTRA = states, uapp, cases, spectra
    cost = ('native evolution CX not synthesized (dense exact reference)' if SERIAL_CX['evolution'] is None
            else f"serial CX = {SERIAL_CX['preparation']} + {SERIAL_CX['evolution']} L")
    print(f'{SELECTED_GEOMETRY}: PREP_M={prep_m} (S0/T), EVOLUTION_M={evolution_m}; {cost}')


select_circuits()

def probability(case, L, trial_gaps):
    if not isinstance(L, (int, np.integer)) or isinstance(L, (bool, np.bool_)) or L < 1:
        raise ValueError('L must be a positive integer')
    states, _ = CASES[case]
    phases, vectors = SPECTRA[case]
    amplitudes = [np.sum(abs(vectors.conj().T @ psi)**2 * np.exp(1j*phases*L)) for psi in states]
    return np.clip(abs(amplitudes[0] + np.exp(1j*np.asarray(trial_gaps)*L*DT)*amplitudes[1])**2 / 4, 0, 1)

def experiment_circuit(L, trial_gap):
    """Selected C-readout circuit; no controlled evolution or dense U^L shortcut."""
    if not isinstance(L, (int, np.integer)) or isinstance(L, (bool, np.bool_)) or L < 1:
        raise ValueError('L must be a positive integer')
    qc = QuantumCircuit(9)
    qc.compose(A_APPROX, inplace=True)
    for _ in range(int(L)):
        qc.compose(EVOLUTION_APPROX, range(8), inplace=True)
    qc.p(float(trial_gap)*L*DT, 8)
    qc.compose(A_APPROX.inverse(), inplace=True)
    return qc

def connect_ibm(backend_name, instance=''):
    prepared = globals().get('IBM_SERVICE')
    if prepared is not None and instance == globals().get('IBM_INSTANCE'):
        return prepared, prepared.backend(backend_name)
    from getpass import getpass
    from qiskit_ibm_runtime import QiskitRuntimeService
    token = os.environ.get('IBM_QUANTUM_API_KEY')
    if not token:
        try:
            from google.colab import userdata
            token = userdata.get('IBM_QUANTUM_API_KEY')
        except Exception:
            token = None
    if not token:
        token = getpass('IBM Quantum API key (非表示): ')
    kwargs = dict(channel='ibm_quantum_platform', token=token)
    if instance:
        kwargs['instance'] = instance
    service = QiskitRuntimeService(**kwargs)
    return service, service.backend(backend_name)

def hardware_preview(backend, L=1):
    if SELECTED_EVOLUTION_M == 'exact':
        raise ValueError('Dense exact evolution is a simulation reference; select EVOLUTION_M=8/16/32 for IBM')
    qc = experiment_circuit(L, GAP)
    qc.measure_all()
    isa = transpile(qc, backend=backend, optimization_level=1, seed_transpiler=7)
    report = {'geometry': SELECTED_GEOMETRY, 'prep_m': SELECTED_PREP_M, 'evolution_m': SELECTED_EVOLUTION_M, 'backend': backend.name, 'L': L, 'depth': isa.depth(),
              'two_qubit_gates': sum(i.operation.num_qubits == 2 for i in isa.data),
              'operations': dict(isa.count_ops())}
    print(report)
    return report

class Measurement:
    """One interface; hardware jobs and counts are persisted before the next iteration."""
    def __init__(self, mode='exact', shots=1024, seed=42, service=None, backend=None,
                 run_hardware=False, max_hardware_L=8, directory='outputs/ibm_jobs'):
        if mode not in ['exact', 'shots', 'ibm']:
            raise ValueError('MODE must be exact, shots, or ibm')
        if mode != 'exact' and (not isinstance(shots, int) or isinstance(shots, bool) or shots < 1):
            raise ValueError('SHOTS must be a positive integer')
        self.mode, self.shots = mode, shots
        self.rng = np.random.default_rng(seed)
        self.service, self.backend = service, backend
        self.run_hardware, self.max_hardware_L = run_hardware, max_hardware_L
        self.directory = Path(directory)
        self.records = []
        self.selection = dict(geometry=SELECTED_GEOMETRY, prep_m=SELECTED_PREP_M, evolution_m=SELECTED_EVOLUTION_M,
                              circuit_digest=CIRCUIT_DIGEST)

    def measure(self, case, L, trial_gaps, label):
        if self.selection['circuit_digest'] != CIRCUIT_DIGEST:
            raise ValueError('Circuit selection changed; recreate Measurement by rerunning the mode setting cell')
        if case not in CASES:
            raise ValueError(f'Unknown CASE: {case}')
        if not isinstance(L, (int, np.integer)) or isinstance(L, (bool, np.bool_)) or L < 1:
            raise ValueError('L must be a positive integer')
        trial_gaps = np.asarray(trial_gaps, dtype=float)
        if trial_gaps.ndim != 1 or trial_gaps.size == 0 or not np.all(np.isfinite(trial_gaps)):
            raise ValueError('trial_gaps must be a nonempty finite one-dimensional array')
        if self.mode != 'ibm':
            p = probability(case, L, trial_gaps)
            counts = None if self.mode == 'exact' else self.rng.binomial(self.shots, p)
            values = p if counts is None else counts/self.shots
            self.records.append(dict(source='exact_probabilities' if counts is None else 'simulated_counts',
                label=label, case=case, L=int(L), **self.selection, trial_gaps=trial_gaps.tolist(),
                shots=None if counts is None else self.shots,
                counts_all_zero=None if counts is None else counts.tolist(), values=values.tolist()))
            return values
        if SELECTED_EVOLUTION_M == 'exact':
            raise ValueError('Dense exact evolution is a simulation reference; select EVOLUTION_M=8/16/32 for IBM')
        if case != 'both':
            raise ValueError('IBM mode uses the selected both circuit; set CASE=both')
        if self.backend is None or self.service is None:
            raise ValueError('Connect IBM and inspect hardware_preview first')
        if L > self.max_hardware_L:
            raise ValueError(f'L={L} exceeds MAX_HARDWARE_L={self.max_hardware_L}; inspect cost before increasing')
        if not self.run_hardware:
            raise ValueError('After hardware_preview, explicitly set RUN_HARDWARE=True to submit/retrieve')
        self.directory.mkdir(parents=True, exist_ok=True)
        spec = dict(assets=ASSET_DIGEST, **self.selection, backend=self.backend.name, L=int(L),
                    trial_gaps=trial_gaps.tolist(), shots=self.shots, label=label)
        key = hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()
        path = self.directory / f'{key}.json'
        if path.exists():
            record = json.loads(path.read_text())
            if 'counts' in record:
                self.records.append(record)
                return np.array(record['values'])
            if 'job_id' not in record:
                raise RuntimeError(f'Unresolved submission marker: {path}. Check IBM jobs before removing it.')
            job = self.service.job(record['job_id'])
        else:
            from qiskit_ibm_runtime import SamplerV2
            circuits = []
            for gap in trial_gaps:
                qc = experiment_circuit(L, gap); qc.measure_all(); circuits.append(qc)
            isa = transpile(circuits, backend=self.backend, optimization_level=1, seed_transpiler=7)
            print('Submitting', len(isa), 'circuits;', self.shots, 'shots each; L=', L)
            record = dict(source='ibm_submission', specification=spec)
            # An exclusive marker prevents a rerun from resubmitting after a submission-time interruption.
            with path.open('x') as stream:
                json.dump(record, stream)
            job = SamplerV2(mode=self.backend).run(isa, shots=self.shots)
            record['job_id'] = job.job_id()
            path.write_text(json.dumps(record, indent=2))
            print('Saved job ID:', record['job_id'], 'in', path)
        raw = [pub.data.meas.get_counts() for pub in job.result()]
        if len(raw) != len(trial_gaps) or any(sum(c.values()) <= 0 for c in raw):
            raise ValueError('Invalid or incomplete IBM counts')
        values = [c.get('0'*9, 0)/sum(c.values()) for c in raw]
        record.update(source='ibm_raw_counts', counts=raw, values=values)
        path.write_text(json.dumps(record, indent=2))
        self.records.append(record)
        return np.array(values)

    def save(self, path):
        path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.records, indent=2))


def gaussian(x, amplitude, mean, sd):
    return amplitude*np.exp(-(x-mean)**2/(2*sd**2))

def choose_measurement_settings(mu, variance, points=21):
    """Keep upstream's numerical convention: time=1.8/v, window=mu +/- v."""
    if not np.isfinite(mu) or not np.isfinite(variance) or variance <= 0:
        raise ValueError('Finite mean and positive variance required')
    if not isinstance(points, int) or points < 5:
        raise ValueError('points must be an integer >=5')
    requested_time = 1.8/variance
    L = max(1, int(requested_time/DT))
    return dict(requested_time=float(requested_time), L=L, time=L*DT,
                trial_gaps=np.linspace(mu-variance, mu+variance, points).tolist())


def fit_gaussian_signal(trial_gaps, observed, prior_mu, prior_variance):
    """Unweighted signal fit; fitted width is not the parameter covariance."""
    eps, y = np.asarray(trial_gaps, float), np.asarray(observed, float)
    if eps.ndim != 1 or y.shape != eps.shape or not np.all(np.isfinite(y)):
        raise ValueError('Finite observed values matching one-dimensional scan points required')
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always', OptimizeWarning)
        params, _ = curve_fit(gaussian, eps, y,
            p0=[1., prior_mu, np.sqrt(prior_variance)], maxfev=10000)
    amplitude, mean, sd = params
    if not np.all(np.isfinite(params)) or amplitude <= 0 or sd == 0:
        raise ValueError(f'Invalid Gaussian fit: {params}')
    fitted_variance = float(sd**2)
    if not np.isfinite(fitted_variance) or fitted_variance <= 0:
        raise ValueError('Invalid fitted variance')
    fitted = gaussian(eps, *params)
    return dict(amplitude=float(amplitude), fit_mu=float(mean), fit_sd=float(abs(sd)),
                fit_variance=fitted_variance, fitted_probability=fitted.tolist(),
                residual_rmse=float(np.sqrt(np.mean((y-fitted)**2))),
                warnings=[str(w.message) for w in caught])


def update_gaussian_prior(prior_mu, prior_variance, fit_mu, fit_variance):
    """Gaussian product; amplitude cancels under normalization."""
    if (not np.all(np.isfinite([prior_mu, prior_variance, fit_mu, fit_variance]))
            or prior_variance <= 0 or fit_variance <= 0):
        raise ValueError('Finite means and positive variances required')
    mu = (prior_variance*fit_mu + fit_variance*prior_mu)/(prior_variance+fit_variance)
    variance = prior_variance*fit_variance/(prior_variance+fit_variance)
    return float(mu), float(variance)


def plot_bayes_round(row, case, mode):
    """A: chosen scan, B: signal fit, C: residual, D: Gaussian product."""
    mu, v = row['prior_mu'], row['prior_variance']
    post_mu, post_v = row['posterior_mu'], row['posterior_variance']
    eps, observed = np.asarray(row['trial_gaps']), np.asarray(row['observed'])
    density = lambda x, mean, var: gaussian(x, 1/np.sqrt(2*np.pi*var), mean, np.sqrt(var))
    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    low, high = min(mu-4*np.sqrt(v), eps[0], GAP), max(mu+4*np.sqrt(v), eps[-1], GAP)
    x = np.linspace(low, high, 2001)
    axes[0,0].plot(x, density(x, mu, v), label='Prior')
    axes[0,0].axvspan(eps[0], eps[-1], alpha=.12, color='C1', label='Scan window')
    axes[0,0].plot(eps, np.zeros_like(eps), '|', color='C1', ms=9, label='Scan points')
    axes[0,0].set(title=f'A. Prior -> scan: L={row["L"]}, t={row["time"]:.4g}', ylabel='Prior density')
    dense = np.linspace(eps[0], eps[-1], 501)
    axes[0,1].plot(dense, probability(case, row['L'], dense), ls='--', label=f'Ideal circuit ({case})')
    axes[0,1].plot(eps, observed, 'o', ms=4, label=f'Observed ({mode})')
    axes[0,1].plot(dense, gaussian(dense,row['amplitude'],row['fit_mu'],row['fit_sd']), label='Gaussian fit')
    axes[0,1].set(title='B. Measurements -> Gaussian fit', ylabel='P(all zero)')
    residual = observed-np.asarray(row['fitted_probability'])
    axes[1,0].axhline(0, color='gray', lw=1)
    axes[1,0].plot(eps, residual, 'o-')
    axes[1,0].set(title=f'C. Fit residual: RMSE={row["residual_rmse"]:.3g}', ylabel='Observed - fit')
    low = min(mu-4*np.sqrt(v), post_mu-4*np.sqrt(post_v), GAP)
    high = max(mu+4*np.sqrt(v), post_mu+4*np.sqrt(post_v), GAP)
    x = np.linspace(low, high, 2001)
    for mean, var, label in [(mu,v,'Prior'),(row['fit_mu'],row['fit_variance'],'Fitted Gaussian (normalized)'),
                             (post_mu,post_v,'Posterior')]:
        axes[1,1].plot(x, density(x,mean,var), label=label)
    axes[1,1].set(title='D. Prior x fitted Gaussian -> posterior', ylabel='Gaussian density')
    for ax in axes.flat:
        ax.axvline(GAP, color='black', ls=':', label='FCI gap (diagnostic)')
        ax.set_xlabel('Trial gap / Ha'); ax.legend(fontsize=7)
    fig.suptitle(f'Round {row["round"]}: MODE={mode}, CASE={case}, R={row["geometry"]}, M={row["prep_m"]}/{row["evolution_m"]}')
    plt.tight_layout();plt.show();plt.close(fig)


def iterative_bayes(measurement, case='exact', rounds=5, initial_mu=0., initial_variance=4.,
                    points=21, max_L=20000, target_sd=None, *, show_plots=True,
                    output_prefix='outputs/exercise5', verbose=True):
    """Upstream Gaussian fit/product iteration, with optional teaching diagnostics.

    output_prefix=None disables output files, for independent noiseless controls.
    No changes to the upstream time/window, initialization or product formulas.
    """
    if (not isinstance(rounds, int) or rounds < 1 or not isinstance(points, int) or points < 5
            or not np.isfinite(initial_mu) or not np.isfinite(initial_variance) or initial_variance <= 0):
        raise ValueError('Use rounds>=1, points>=5, initial_variance>0')
    mu, variance = float(initial_mu), float(initial_variance)
    history = []
    def say(*args):
        if verbose: print(*args)
    def save_history():
        if output_prefix is not None:
            path = Path(str(output_prefix)+'_updates.json');path.parent.mkdir(parents=True,exist_ok=True)
            path.write_text(json.dumps(history,indent=2))
    def save_measurements():
        if output_prefix is not None: measurement.save(str(output_prefix)+'_measurements.json')
    save_history()
    for iteration in range(1, rounds+1):
        plan = choose_measurement_settings(mu, variance, points)
        L = plan['L']
        if L > max_L:
            say(f'Stop before round {iteration}: requested L={L} > MAX_L={max_L}')
            break
        if measurement.mode == 'ibm' and L > measurement.max_hardware_L:
            say(f'Stop before round {iteration}: requested L={L} > MAX_HARDWARE_L={measurement.max_hardware_L}; '
                'previous results are retained. Inspect hardware_preview before increasing the limit.')
            break
        observed = measurement.measure(case, L, plan['trial_gaps'], f'round_{iteration}')
        save_measurements()
        try:
            fitted = fit_gaussian_signal(plan['trial_gaps'], observed, mu, variance)
        except (ValueError, RuntimeError, FloatingPointError) as exc:
            say('Fit failed; keep measured data, stop updates:', exc)
            break
        post_mu, post_variance = update_gaussian_prior(mu, variance, fitted['fit_mu'], fitted['fit_variance'])
        row = dict(round=iteration, prior_mu=mu, prior_variance=variance, **plan,
                   observed=observed.tolist(), **fitted,
                   posterior_mu=post_mu, posterior_variance=post_variance,
                   mode=measurement.mode, case=case, **measurement.selection, source=measurement.records[-1]['source'],
                   shots=None if measurement.mode=='exact' else measurement.shots)
        history.append(row);save_history()
        say(f'round {iteration}: prior=({mu:.7g},{variance:.7g}), L={L}, t={plan["time"]:.4g}, '
            f'fit mean={fitted["fit_mu"]:.7g}, fit SD={fitted["fit_sd"]:.5g}, '
            f'posterior mean={post_mu:.7f}, SD={np.sqrt(post_variance):.5g}, '
            f'RMSE={fitted["residual_rmse"]:.4g}')
        for warning in fitted['warnings']:say('Gaussian fit warning:',warning)
        if show_plots:plot_bayes_round(row,case,measurement.mode)
        mu, variance = post_mu, post_variance
        if target_sd is not None and np.sqrt(variance) <= target_sd:
            say('Requested SD reached');break
    return history
