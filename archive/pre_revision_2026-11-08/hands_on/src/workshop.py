"""Participant-facing molecular workshop tools. No IBM jobs at import time."""
from pathlib import Path
from functools import lru_cache
import json
import hashlib
from datetime import datetime, timezone
import numpy as np
from scipy.special import logsumexp
from scipy.linalg import expm
from qiskit import QuantumCircuit,qpy,transpile
from qiskit.quantum_info import Statevector,Operator

ROOT=Path(__file__).resolve().parents[1]
HARTREE_TO_EV=27.211386245988


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')
    return path


def load_model(name='ethylene'):
    if name not in ('ethylene','butadiene'):raise ValueError(name)
    with np.load(ROOT/f'data/reference/{name}.npz',allow_pickle=False) as source:
        arrays={key:source[key].copy() for key in source.files}
    return read_json(ROOT/f'data/reference/{name}.json'),arrays


def load_circuit(name,key):
    if name not in ('ethylene','butadiene'):raise ValueError(name)
    if key not in ('prep_reference','prep_compact','evol_reference','evol_compact'):raise ValueError(key)
    path=ROOT/f'data/circuits/{name}_{key}.qpy'
    with path.open('rb') as stream:return qpy.load(stream)[0]


def circuit_summary(name):
    return read_json(ROOT/f'data/circuits/{name}.json')


def qpde_circuit(name,epsilon,steps,variant='compact',measure=True):
    if variant not in ('reference','compact'):raise ValueError(variant)
    if not isinstance(steps,(int,np.integer)) or steps<1:raise ValueError('steps must be a positive integer')
    prep=load_circuit(name,f'prep_{variant}')
    evol=load_circuit(name,f'evol_{variant}')
    dt=circuit_summary(name)['dt_atomic_units']
    n=evol.num_qubits
    qc=QuantumCircuit(n+1)
    qc.compose(prep,inplace=True)
    for _ in range(steps):qc.compose(evol,qubits=list(range(n)),inplace=True)
    qc.p(float(epsilon*steps*dt),n)
    qc.compose(prep.inverse(),inplace=True)
    qc.name=f'{name}_{variant}_m{steps}_eps{epsilon:.4f}'
    qc.metadata={'molecule':name,'variant':variant,'steps':int(steps),'epsilon_hartree':float(epsilon),
                 'time_atomic_units':float(dt*steps)}
    if measure:qc.measure_all()
    return qc


@lru_cache(maxsize=8)
def _prepared_state(name,variant):
    return Statevector.from_instruction(load_circuit(name,f'prep_{variant}')).data


@lru_cache(maxsize=8)
def _step_operator(name,variant):
    return Operator(load_circuit(name,f'evol_{variant}')).data


def zero_probability(name,epsilon,steps,variant='compact'):
    # Exactly the same prep and repeated step circuits as qpde_circuit; evaluate only its all-zero amplitude.
    prep=_prepared_state(name,variant)
    u=_step_operator(name,variant)
    u=np.linalg.matrix_power(u,int(steps))
    d=len(u);v0,v1=prep[:d],prep[d:]
    time=circuit_summary(name)['dt_atomic_units']*steps
    amplitude=np.vdot(v0,u@v0)+np.exp(1j*epsilon*time)*np.vdot(v1,u@v1)
    return float(np.clip(abs(amplitude)**2,0,1))


def make_plan(name='ethylene',variant='compact',team='team01',shots=1024,step_counts=None):
    """Fixed broad scan; no exact gap is used to choose the scan center."""
    if name not in ('ethylene','butadiene'):raise ValueError(name)
    if variant not in ('compact','reference'):raise ValueError(variant)
    if shots<1:raise ValueError('shots must be positive')
    if step_counts is None:step_counts=[2,4,8] if name=='ethylene' else [8,16,32]
    summary=circuit_summary(name)
    rows=[]
    for iteration,steps in enumerate(step_counts,1):
        for epsilon in np.linspace(0,.5,9):
            rows.append({'iteration':iteration,'steps':int(steps),'epsilon_hartree':float(epsilon),
                         'time_atomic_units':float(steps*summary['dt_atomic_units']),'shots':int(shots)})
    fingerprints={key:hashlib.sha256((ROOT/f'data/circuits/{name}_{key}.qpy').read_bytes()).hexdigest()
                  for key in (f'prep_{variant}',f'evol_{variant}')}
    return {'schema_version':1,'molecule':name,'variant':variant,'team':team,'rows':rows,
            'system_qubits':4 if name=='ethylene' else 8,'total_qubits':5 if name=='ethylene' else 9,
            'circuit_sha256':fingerprints,'source':'plan_only','scan_range_hartree':[0,.5],
            'maximum_time_atomic_units':max(r['time_atomic_units'] for r in rows)}


def simulate_results(plan,seed=20261109,noise_rate=.15):
    if not 0<=noise_rate<=1:raise ValueError('noise_rate must be in [0,1]')
    rng=np.random.default_rng(seed)
    result=json.loads(json.dumps(plan));d=2**plan['total_qubits']
    # Cache circuit-to-array conversion once per run; avoid repeating costly decompositions.
    v=Statevector.from_instruction(load_circuit(plan['molecule'],f"prep_{plan['variant']}")).data
    u=Operator(load_circuit(plan['molecule'],f"evol_{plan['variant']}")).data
    s=len(u);v0,v1=v[:s],v[s:];amplitudes={}
    for steps in {row['steps'] for row in plan['rows']}:
        power=np.linalg.matrix_power(u,steps)
        amplitudes[steps]=(np.vdot(v0,power@v0),np.vdot(v1,power@v1))
    for row in result['rows']:
        a,b=amplitudes[row['steps']]
        p=abs(a+np.exp(1j*row['epsilon_hartree']*row['time_atomic_units'])*b)**2
        p=float(np.clip((1-noise_rate)*p+noise_rate/d,0,1))
        row['successes']=int(rng.binomial(row['shots'],p))
        row['simulated_probability']=p
    result.update(source='simulated_counts',seed=seed,noise_model={'name':'global_depolarization_after_circuit',
                  'rate':noise_rate,'not_ibm_measurements':True})
    return result


def validate_counts(data):
    if data.get('schema_version')!=1:raise ValueError('Unsupported data schema')
    if data.get('source') not in ('simulated_counts','ibm_raw_counts'):
        raise ValueError('Counts data required; a plan or mitigated quasiprobabilities cannot be used as binomial counts')
    if not data.get('rows'):raise ValueError('No rows')
    for row in data['rows']:
        shots,k=row['shots'],row['successes']
        if not isinstance(shots,int) or not isinstance(k,int) or not 0<=k<=shots or shots<1:
            raise ValueError('Invalid integer counts')
        if not (np.isfinite(row['epsilon_hartree']) and row['time_atomic_units']>0):raise ValueError('Invalid setting')
    return data


def infer_gap(data,noise_aware=True):
    """Sequential grid Bayes, marginalizing shared nuisance parameters consistently.
    p = floor + visibility*(1+cos((gap-epsilon)*t))/2.
    Nuisance parameters are shared across iterations, a teaching assumption.
    """
    validate_counts(data)
    grid=np.linspace(0,.5,1001)
    if noise_aware:
        visibilities=np.linspace(.2,1.,17)
        floors=np.array([0.,.005,.02,.05])
    else:visibilities=np.array([1.]);floors=np.array([0.])
    valid=visibilities[:,None]+floors[None,:]<=1.+1e-12
    joint=np.zeros((len(grid),len(visibilities),len(floors)))
    joint[:,~valid]=-np.inf
    joint-=logsumexp(joint)
    history=[]
    iterations=sorted({r['iteration'] for r in data['rows']})
    for iteration in iterations:
        for row in [r for r in data['rows'] if r['iteration']==iteration]:
            fringe=(1+np.cos((grid-row['epsilon_hartree'])*row['time_atomic_units']))/2
            p=floors[None,None,:]+visibilities[None,:,None]*fringe[:,None,None]
            p=np.clip(p,1e-12,1-1e-12)
            k,n=row['successes'],row['shots']
            joint+=k*np.log(p)+(n-k)*np.log1p(-p)
        joint-=logsumexp(joint)
        marginal=np.exp(logsumexp(joint,axis=(1,2)))
        cdf=np.cumsum(marginal)
        low,high=np.interp([.025,.975],cdf,grid)
        mean=float(np.sum(grid*marginal))
        history.append({'iteration':iteration,'grid':grid.copy(),'posterior':marginal,
                        'mean_hartree':mean,'ci95_hartree':[float(low),float(high)],
                        'edge_mass':float(marginal[:20].sum()+marginal[-20:].sum())})
    return history


def atomic_write_record(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+'.tmp')
    write_json(temporary,data);temporary.replace(path)
    return path


def submit_ibm(plan,service,backend_name,record_path):
    """Explicit participant action. Never called by default notebook execution.
    Save an execution record immediately after each job submission. One job/iteration.
    """
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
    from qiskit_ibm_runtime import SamplerV2
    backend=service.backend(backend_name)
    if backend.num_qubits<plan['total_qubits']:raise ValueError('Backend too small')
    record=json.loads(json.dumps(plan))
    record.update(source='ibm_submission',backend=backend_name,submitted_at=datetime.now(timezone.utc).isoformat(),jobs=[])
    if Path(record_path).exists():raise FileExistsError('Record exists; recover the existing job or choose a new path explicitly')
    manager=generate_preset_pass_manager(backend=backend,optimization_level=3,seed_transpiler=11)
    sampler=SamplerV2(mode=backend,options={'default_shots':int(plan['rows'][0]['shots'])})
    for iteration in sorted({r['iteration'] for r in plan['rows']}):
        indices=[i for i,r in enumerate(plan['rows']) if r['iteration']==iteration]
        circuits=[qpde_circuit(plan['molecule'],plan['rows'][i]['epsilon_hartree'],plan['rows'][i]['steps'],plan['variant']) for i in indices]
        isa=manager.run(circuits)
        job=sampler.run(isa,shots=int(plan['rows'][indices[0]]['shots']))
        entry={'iteration':iteration,'job_id':job.job_id(),'row_indices':indices,
               'isa_depths':[q.depth() for q in isa],
               'isa_two_qubit_counts':[sum(1 for op in q.data if len(op.qubits)==2) for q in isa]}
        record['jobs'].append(entry)
        atomic_write_record(record_path,record)
    return record


def fetch_ibm(record,service):
    """Fetch completed raw Sampler jobs without resubmitting; fail on partial/pending data."""
    if record.get('source')!='ibm_submission':raise ValueError('Submission record required')
    result=json.loads(json.dumps(record))
    for entry in record['jobs']:
        job=service.job(entry['job_id'])
        status=str(job.status()).upper()
        if 'DONE' not in status:raise RuntimeError(f"Job {entry['job_id']} status={status}; retry later without resubmitting")
        pubs=job.result()
        if len(pubs)!=len(entry['row_indices']):raise ValueError('Job result size mismatch')
        for index,pub in zip(entry['row_indices'],pubs):
            counts=pub.data.meas.get_counts()
            total=sum(counts.values())
            zeros='0'*record['total_qubits']
            result['rows'][index]['successes']=int(counts.get(zeros,0))
            result['rows'][index]['shots']=int(total)
            result['rows'][index]['counts']=counts
    retrieved=sorted({i for job in record['jobs'] for i in job['row_indices']})
    result['rows']=[result['rows'][i] for i in retrieved]
    result['requested_row_count']=len(record['rows'])
    result['retrieved_row_count']=len(result['rows'])
    result['partial_plan']=len(result['rows'])!=len(record['rows'])
    result['source']='ibm_raw_counts'
    return validate_counts(result)


def plot_molecule(name,ax=None):
    import matplotlib.pyplot as plt
    meta,_=load_model(name)
    if ax is None:_,ax=plt.subplots(figsize=(7,3))
    atoms=meta['atoms_angstrom']
    for a,b,order in meta['bonds']:
        pa=np.array(atoms[a][1][:2]);pb=np.array(atoms[b][1][:2])
        normal=np.array([-(pb-pa)[1],(pb-pa)[0]])
        normal/=np.linalg.norm(normal)
        for offset in ([-.045,.045] if order==2 else [0]):
            ax.plot([pa[0]+offset*normal[0],pb[0]+offset*normal[0]],
                    [pa[1]+offset*normal[1],pb[1]+offset*normal[1]],color='gray',lw=1.6,zorder=0)
    for atom,position in atoms:
        ax.scatter(*position[:2],s=180 if atom=='C' else 80,c='#334155' if atom=='C' else '#dbeafe',edgecolors='white',zorder=2)
        ax.text(position[0],position[1]+.15,atom,ha='center',fontsize=10)
    positions=np.array([entry[1][:2] for entry in atoms])
    ax.set_xlim(positions[:,0].min()-.35,positions[:,0].max()+.35)
    ax.set_ylim(positions[:,1].min()-.35,positions[:,1].max()+.65)
    ax.set_aspect('equal');ax.axis('off');ax.set_title(name+' (fixed teaching geometry)',pad=12)
    return ax


def plot_posterior(history,reference=None):
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7,3.5))
    for row in history:
        grid=row['grid'];density=row['posterior']/(grid[1]-grid[0])
        ax.plot(grid*HARTREE_TO_EV,density/HARTREE_TO_EV,label=f"After round {row['iteration']}")
    if reference is not None:ax.axvline(reference*HARTREE_TO_EV,color='black',ls=':',label='Classical model reference')
    ax.set(xlabel='Singlet-triplet gap (eV)',ylabel='Posterior density (1/eV)')
    ax.legend();fig.tight_layout()
    return fig


def inference_diagnostics(data, history):
    """Exploratory signal and model checks; not a guarantee of hardware accuracy."""
    from scipy.stats import chi2
    validate_counts(data)
    deviances=[]
    def deviance(k,n,p):
        k=np.asarray(k,dtype=float);n=np.asarray(n,dtype=float);p=np.clip(p,1e-12,1-1e-12)
        expected=n*p
        first=np.zeros_like(k);last=np.zeros_like(k)
        nonzero=k>0;notfull=n-k>0
        first[nonzero]=k[nonzero]*np.log(k[nonzero]/expected[nonzero])
        last[notfull]=(n-k)[notfull]*np.log((n-k)[notfull]/(n-expected)[notfull])
        return max(0.,float(2*np.sum(first+last)))
    for iteration in sorted({r['iteration'] for r in data['rows']}):
        rows=[r for r in data['rows'] if r['iteration']==iteration]
        k=np.array([r['successes'] for r in rows]);n=np.array([r['shots'] for r in rows])
        d=deviance(k,n,np.full(len(k),k.sum()/n.sum()))
        deviances.append({'iteration':iteration,'flat_signal_deviance':d,
                         'approx_flat_signal_p':float(chi2.sf(d,max(1,len(k)-1)))})
    # Minimize binomial deviance over the same nuisance grid at posterior mean gap.
    gap=history[-1]['mean_hartree'];rows=data['rows']
    k=np.array([r['successes'] for r in rows]);n=np.array([r['shots'] for r in rows])
    fringe=np.array([(1+np.cos((gap-r['epsilon_hartree'])*r['time_atomic_units']))/2 for r in rows])
    best=(np.inf,None,None)
    for visibility in np.linspace(.2,1.,17):
        for floor in [0.,.005,.02,.05]:
            if visibility+floor>1+1e-12:continue
            d=deviance(k,n,floor+visibility*fringe)
            if d<best[0]:best=(d,float(visibility),floor)
    tail=float(chi2.sf(best[0],max(1,len(rows)-3)))
    return {'round_signal_checks':deviances,'detectable_scan_variation':any(r['approx_flat_signal_p']<.01 for r in deviances),
            'deviance_at_posterior_mean_with_best_nuisance':float(best[0]),
            'fitted_visibility':best[1],'fitted_floor':best[2],
            'approx_model_tail_probability':tail,'model_mismatch_warning':tail<.01,
            'interpretation':'Exploratory chi-square checks; low counts, boundaries and finite grids limit calibration. Scan variation alone does not establish a correct gap.'}
