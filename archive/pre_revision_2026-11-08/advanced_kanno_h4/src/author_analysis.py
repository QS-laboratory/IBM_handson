"""Read upstream brick-wall gates and reproduce its QPDE Gaussian-fit updates.
The Julia optimizer is used unchanged. Only execution/plotting are isolated from
legacy IBM authentication. Arrays follow Qiskit's native little endian; ancilla q0.
"""
from pathlib import Path
import ast,json,hashlib
import warnings
import numpy as np
from scipy.linalg import expm
from scipy.optimize import curve_fit, OptimizeWarning
from qiskit import QuantumCircuit,transpile,qpy
from qiskit.quantum_info import Statevector,Operator
ROOT=Path(__file__).resolve().parents[1]


def gaussian_from_author():
    source=ROOT/'upstream/src/QPDE_Bayes_session_MPSBWMPO_RD.py'
    tree=ast.parse(source.read_text())
    node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='gaussian')
    env={'np':np}
    exec(compile(ast.Module(body=[node],type_ignores=[]),str(source),'exec'),env)
    return env['gaussian']


gaussian=gaussian_from_author()


def load_meta(name='canonical',geometry=0):
    return json.loads((ROOT/f'data/geometry_{geometry}/{name}.json').read_text())


def load_arrays(name='canonical',geometry=0):
    with np.load(ROOT/f'data/geometry_{geometry}/{name}.npz',allow_pickle=False) as data:
        return {k:data[k].copy() for k in data.files}


def brickwall(name,kind,geometry=0,depth=8):
    if kind not in ('prep','evol'):raise ValueError(kind)
    qc=QuantumCircuit(9 if kind=='prep' else 8)
    root=ROOT/f'results/geometry_{geometry}'
    for layer in range(depth):
        count=4 if kind=='prep' or layer%2==0 else 3
        for unit in range(count):
            suffix='MPSBW' if kind=='prep' else 'MPO'
            path=root/f'{name}_{suffix}_depth{layer+1}_{unit+1}.npy'
            matrix=np.load(path,allow_pickle=False)
            if np.linalg.norm(matrix.conj().T@matrix-np.eye(4))>1e-8:
                raise ValueError(f'Gate is not unitary: {path}')
            start=2*unit+layer%2
            qc.unitary(matrix,[start,start+1])
    return qc


def qpde_circuit(name,epsilon,steps,geometry=0,measure=False):
    prep=brickwall(name,'prep',geometry)
    evol=brickwall(name,'evol',geometry)
    qc=QuantumCircuit(9)
    qc.compose(prep,inplace=True)
    for _ in range(int(steps)):qc.compose(evol,qubits=list(range(1,9)),inplace=True)
    qc.p(.1*steps*epsilon,0)
    qc.compose(prep.inverse(),inplace=True)
    if measure:qc.measure_all()
    return qc


def dense_second_trotter_adjoint(arrays, slices=100, dt=.1):
    """Independent dense reference with the same Pauli order, without MPO/SVD truncation."""
    h=arrays['hamiltonian'];dim=len(h);n=int(np.log2(dim))
    matrix=np.eye(dim,dtype=complex);indices=np.arange(dim)
    terms=list(zip(arrays['pauli_labels'],arrays['pauli_coefficients']))
    for label,coefficient in terms+list(reversed(terms)):
        mask=0;phase=np.ones(dim,complex)
        for q,axis in enumerate(reversed(label)):
            bit=(indices>>q)&1
            if axis in 'XY':mask|=1<<q
            if axis=='Y':phase*=1j*(1-2*bit)
            elif axis=='Z':phase*=1-2*bit
        applied=np.empty_like(matrix);applied[indices^mask]=phase[:,None]*matrix
        angle=float(coefficient.real)*dt/slices/2
        matrix=np.cos(angle)*matrix+1j*np.sin(angle)*applied
    return np.linalg.matrix_power(matrix,slices)


def evaluate(name,geometry=0):
    arrays=load_arrays(name,geometry);h=arrays['hamiltonian']
    folder=ROOT/f'results/geometry_{geometry}'
    prep=brickwall(name,'prep',geometry)
    evolved=brickwall(name,'evol',geometry)
    vector=Statevector.from_instruction(prep).data
    target=np.load(folder/f'{name}_target_mps_vector.npy',allow_pickle=False)
    u=Operator(evolved).data
    ref=np.load(folder/f'{name}_reference_mpo_dense.npy',allow_pickle=False)
    exact=expm(-1j*h*.1)
    overlap=abs(np.vdot(target,vector))
    history=np.load(folder/f'{name}_prep_history.npy',allow_pickle=False)
    delta_history=np.load(folder/f'{name}_evol_history.npy',allow_pickle=False)
    # Author MPO stores the adjoint reference in its contraction orientation.
    ref_adjoint_error=float(np.linalg.norm(ref-exact.conj().T)/16)
    dense_trotter=dense_second_trotter_adjoint(arrays)
    dense_trotter_error=float(np.linalg.norm(dense_trotter-exact.conj().T)/16)
    mpo_truncation_error=float(np.linalg.norm(ref-dense_trotter)/16)
    score=float(np.trace(ref@u).real)
    delta=float(np.sqrt(max(0,2-max(0,score)**(1/8))))
    ground=target[::2]*np.sqrt(2);excited=target[1::2]*np.sqrt(2)
    eg=float(np.vdot(ground,h@ground).real);ee=float(np.vdot(excited,h@excited).real)
    number=np.array([i.bit_count() for i in range(256)])
    sp=np.zeros((256,256))
    sz=np.zeros(256)
    for i in range(256):
        for orbital in range(4):
            up=(i>>(2*orbital))&1; down=(i>>(2*orbital+1))&1
            sz[i]+=.5*(up-down)
            if down and not up:sp[i^(3<<(2*orbital)),i]+=1
    s2=np.diag(sz**2)+.5*(sp@sp.T+sp.T@sp)
    metrics={'geometry_id':geometry,'orbitals':name,'depth_prep':8,'depth_evol':8,'optimization_sweeps':len(history)-1,
             'dmrg_sweeps':200,'dt':.1,'trotter_slices':100,'seed':1,
             'preparation_overlap':float(overlap),'preparation_fidelity_squared':float(overlap**2),
             'author_preparation_final_cost':float(history[-1]),'author_evolution_final_delta':float(delta_history[-1]),
             'evolution_delta_recomputed':delta,'evolution_normalized_frobenius_error':float(np.linalg.norm(u-exact)/16),
             'reference_mpo_adjoint_error_vs_matrix_exponential':ref_adjoint_error,
             'dense_second_trotter_error_vs_exponential':dense_trotter_error,
             'reference_mpo_vs_dense_trotter_error':mpo_truncation_error,
             'reference_mpo_unitarity_defect':float(np.linalg.norm(ref.conj().T@ref-np.eye(256))/16),
             'prepared_ground_dmrg_energy':eg,'prepared_excited_dmrg_energy':ee,'dmrg_gap':ee-eg,
             'exact_gap':load_meta(name,geometry)['gap_exact_hartree'],
             'dmrg_ground_residual':float(np.linalg.norm(h@ground-eg*ground)),
             'dmrg_excited_residual':float(np.linalg.norm(h@excited-ee*excited)),
             'dmrg_particle_numbers':[float(np.vdot(v,number*v).real) for v in (ground,excited)],
             'dmrg_s_squared':[float(np.vdot(v,s2@v).real) for v in (ground,excited)],
             'mapped_qubit_exact_gap':float(np.diff(arrays['full_eigenvalues'][:2])[0])}
    # Prove that NumPy gate export -> Qiskit uses the actual author convention.
    if abs(overlap-history[-1])>1e-6:raise AssertionError(f'Preparation export/order mismatch: {overlap} vs {history[-1]}')
    if abs(delta-delta_history[-1])>1e-6:raise AssertionError(f'Evolution export/order mismatch: {delta} vs {delta_history[-1]}')
    assert dense_trotter_error<1e-6, dense_trotter_error
    # The unchanged author MPO pipeline additionally truncates; quantify rather than call it exact.
    assert ref_adjoint_error<1e-3, ref_adjoint_error
    assert abs(ee-eg-metrics['exact_gap'])<1e-7
    for kind,qc in [('prep',prep),('evol',evolved)]:
        compiled=transpile(qc,basis_gates=['rz','sx','x','cx'],optimization_level=3,seed_transpiler=1)
        metrics[f'{kind}_two_qubit_blocks']=len(qc.data)
        metrics[f'{kind}_cx_count_all_to_all']=int(compiled.count_ops().get('cx',0))
        metrics[f'{kind}_depth_all_to_all']=int(compiled.depth())
        with (folder/f'{name}_{kind}.qpy').open('wb') as stream:qpy.dump(qc,stream)
    (folder/f'{name}_metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
    return metrics


def run_qpde(name,geometry=0,shots=10000,max_iterations=30,seed=1,exact_probabilities=False,reference_control=False):
    """Same t=1.8/variance, mu±variance, 21 points, fit and Gaussian update as upstream.
    All-zero counts use their exact binomial marginal, matching the shot distribution
    of an ideal statevector simulator; seed streams differ from Aer.
    """
    if reference_control:
        vector=np.load(ROOT/f'results/geometry_{geometry}/{name}_target_mps_vector.npy',allow_pickle=False)
        u=expm(-1j*load_arrays(name,geometry)['hamiltonian']*.1)
    else:
        vector=Statevector.from_instruction(brickwall(name,'prep',geometry)).data
        u=Operator(brickwall(name,'evol',geometry)).data
    v0,v1=vector[::2],vector[1::2]
    rng=np.random.default_rng(seed);mu=0.;variance=4.;rows=[];status='max_iterations'
    for iteration in range(1,max_iterations+1):
        requested_time=1.8/variance
        steps=int(requested_time/.1)
        time=.1*steps
        eps=np.linspace(mu-variance,mu+variance,21)
        power=np.linalg.matrix_power(u,steps)
        a,b=np.vdot(v0,power@v0),np.vdot(v1,power@v1)
        probability=np.clip(abs(a+np.exp(1j*eps*time)*b)**2,0,1)
        counts=None if exact_probabilities else rng.binomial(shots,probability)
        observed=probability if exact_probabilities else counts/shots
        warning_messages=[]
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always',OptimizeWarning)
                params,covariance=curve_fit(gaussian,eps,observed,p0=[1.,mu,np.sqrt(variance)],maxfev=10000)
                warning_messages=[str(w.message) for w in caught]
        except (RuntimeError,ValueError,FloatingPointError) as exc:
            status='gaussian_fit_failed'
            rows.append({'iteration':iteration,'prior_mean':mu,'prior_variance':variance,'time':time,'steps':steps,
                         'epsilon':eps.tolist(),'observed_probability':observed.tolist(),'error':str(exc)})
            break
        amplitude,mean,std=params
        likelihood_variance=std**2
        posterior_mu=(variance*mean+likelihood_variance*mu)/(variance+likelihood_variance)
        posterior_variance=variance*likelihood_variance/(variance+likelihood_variance)
        rows.append({'iteration':iteration,'prior_mean':float(mu),'prior_variance':float(variance),'requested_time':requested_time,
                     'time':time,'steps':steps,'epsilon':eps.tolist(),'exact_probability':probability.tolist(),
                     'counts_all_zero':None if counts is None else counts.tolist(),'observed_probability':observed.tolist(),
                     'gaussian_fit_amplitude':float(amplitude),'likelihood_mean':float(mean),'likelihood_variance':float(likelihood_variance),
                     'posterior_mean':float(posterior_mu),'posterior_variance':float(posterior_variance),
                     'warnings':warning_messages})
        mu,variance=posterior_mu,posterior_variance
        if not np.isfinite(mu+variance) or variance<=0:status='invalid_gaussian_fit';break
        if variance<.005:status='variance_threshold_reached';break
    result={'source':'noiseless_exact_probabilities' if exact_probabilities else 'noiseless_simulated_shot_counts',
            'orbitals':name,'geometry_id':geometry,'shots':shots,'seed':seed,'status':status,'iterations':rows,'reference_control':reference_control,
            'exact_gap':load_meta(name,geometry)['gap_exact_hartree'],
            'method':'Author Gaussian fit + Gaussian product update, not the main workshop grid-binomial posterior.',
            'interpretation':'posterior variance is the author Gaussian algorithm parameter, NOT a full-error confidence guarantee',
            'execution_adapter':'All-zero amplitude evaluation and binomial marginal; no IBM/Fire Opal/authentication calls.'}
    suffix=('reference_control_' if reference_control else '')+('exact' if exact_probabilities else 'shots')
    (ROOT/f'results/geometry_{geometry}/{name}_qpde_{suffix}.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
