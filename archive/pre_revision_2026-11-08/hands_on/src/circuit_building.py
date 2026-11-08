"""Offline circuit compilation for the small molecular teaching models."""
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import expm
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import StatePreparation
from qiskit.quantum_info import Statevector, Operator


def single_ry(state, angle, q):
    out=state.copy().reshape(-1,2,1<<q)
    a=out[:,0,:].copy(); b=out[:,1,:].copy()
    c,s=np.cos(angle/2),np.sin(angle/2)
    out[:,0,:]=c*a-s*b; out[:,1,:]=s*a+c*b
    return out.reshape(state.shape)


def ry_derivative_post(state,q):
    out=state.copy().reshape(-1,2,1<<q)
    a=out[:,0,:].copy(); b=out[:,1,:].copy()
    out[:,0,:]=-.5*b; out[:,1,:]=.5*a
    return out.reshape(state.shape)


def prep_operations(n,layers):
    operations=[]; index=0
    for layer in range(layers):
        for q in range(n):
            operations.append(('ry',q,index)); index+=1
        # Alternating nearest-neighbor sweeps allow both directions of information flow.
        pairs=list(range(n-1)) if layer%2==0 else list(reversed(range(n-1)))
        for q in pairs:
            operations.append(('cx',q,q+1) if layer%2==0 else ('cx',q+1,q))
    return operations


def cx_permutation(n,control,target):
    indices=np.arange(1<<n)
    return indices ^ (((indices>>control)&1)<<target)


def preparation_circuit(parameters,n,layers):
    qc=QuantumCircuit(n)
    for kind,a,b in prep_operations(n,layers):
        if kind=='ry': qc.ry(float(parameters[b]),a)
        else: qc.cx(a,b)
    return qc


def fit_preparation(target,layers,maxiter=1800,seed=19):
    n=int(np.log2(len(target))); ops=prep_operations(n,layers)
    perms={(a,b):cx_permutation(n,a,b) for kind,a,b in ops if kind=='cx'}
    target=np.asarray(target,dtype=float)
    def loss(parameters):
        state=np.zeros(1<<n);state[0]=1.
        history=[]
        for kind,a,b in ops:
            state=single_ry(state,parameters[b],a) if kind=='ry' else state[perms[a,b]]
            history.append(state)
        adj=state-target
        value=float(np.dot(adj,adj))
        gradient=np.zeros_like(parameters)
        for (kind,a,b),post in zip(reversed(ops),reversed(history)):
            if kind=='ry':
                gradient[b]=2*np.dot(adj,ry_derivative_post(post,a))
                adj=single_ry(adj,-parameters[b],a)
            else: adj=adj[perms[a,b]]
        return value,gradient
    rng=np.random.default_rng(seed)
    best=None
    for attempt in range(2):
        fit=minimize(loss,rng.normal(0,.8,n*layers),jac=True,method='L-BFGS-B',
                     options={'maxiter':maxiter,'ftol':1e-13,'gtol':1e-8,'maxls':30})
        if best is None or fit.fun<best.fun: best=fit
        if best.fun<1e-10: break
    qc=preparation_circuit(best.x,n,layers)
    fidelity=float(abs(np.vdot(target,Statevector.from_instruction(qc).data))**2)
    return qc,best.x,{'layers':layers,'fidelity':fidelity,'optimizer_success':bool(best.success),
                     'iterations':int(best.nit),'squared_amplitude_error':float(best.fun)}


def append_pauli(qc,label,angle):
    """Append exp(-i angle P) using basis changes, CNOTs and RZ(2 angle)."""
    support=[]
    for q,axis in enumerate(reversed(label)):
        if axis=='I': continue
        support.append(q)
        if axis=='X':qc.h(q)
        elif axis=='Y':qc.sdg(q);qc.h(q)
    if not support:
        qc.global_phase-=angle;return
    for a,b in zip(support[:-1],support[1:]):qc.cx(a,b)
    qc.rz(2*float(angle),support[-1])
    for a,b in reversed(list(zip(support[:-1],support[1:]))):qc.cx(a,b)
    for q in support:
        axis=label[-1-q]
        if axis=='X':qc.h(q)
        elif axis=='Y':qc.h(q);qc.s(q)


def ordered_terms(labels,coeffs):
    order=sorted(range(len(labels)),key=lambda i:(not set(labels[i])<=set('IZ'),labels[i]))
    return [(str(labels[i]),float(coeffs[i]),int(not set(labels[i])<=set('IZ'))) for i in order]


def evolution_circuit(labels,coeffs,times):
    n=len(labels[0]);qc=QuantumCircuit(n)
    terms=ordered_terms(labels,coeffs)
    for pair in np.asarray(times).reshape(-1,2):
        for label,c,group in terms:append_pauli(qc,label,c*pair[group])
    return qc


def pauli_action_metadata(label):
    n=len(label);indices=np.arange(1<<n);mask=0;phase=np.ones(1<<n,dtype=complex)
    for q,axis in enumerate(reversed(label)):
        bit=(indices>>q)&1
        if axis in 'XY':mask|=1<<q
        if axis=='Y':phase*=1j*(1-2*bit)
        elif axis=='Z':phase*=1-2*bit
    return indices^mask,phase


def apply_pauli(value,meta):
    permutation,phase=meta
    out=np.empty_like(value)
    out[permutation]=phase[:,None]*value if value.ndim==2 else phase*value
    return out


def fit_evolution(h,labels,coeffs,dt,layers=1,seed=23):
    """Fit grouped Trotter times against dense reference actions on fixed probes.
    Dense reference/probes are teaching tools, not scalable tensor contractions.
    """
    dim=len(h);rng=np.random.default_rng(seed)
    # Exact full operator objective for 4 qubits, 12 fixed random probes for 8 qubits.
    if dim<=16:probes=np.eye(dim,dtype=complex)
    else:
        probes=rng.normal(size=(dim,12))+1j*rng.normal(size=(dim,12))
        probes/=np.linalg.norm(probes,axis=0)
    target=expm(-1j*h*dt)@probes
    terms=ordered_terms(labels,coeffs)
    actions=[pauli_action_metadata(label) for label,_,_ in terms]
    scale=probes.shape[1]
    def loss(parameters):
        value=probes.copy();history=[]
        for layer in range(layers):
            for (_,c,group),action in zip(terms,actions):
                angle=c*parameters[2*layer+group]
                value=np.cos(angle)*value-1j*np.sin(angle)*apply_pauli(value,action)
                history.append(value)
        adj=value-target
        objective=float(np.vdot(adj,adj).real/scale)
        gradient=np.zeros_like(parameters);index=len(history)-1
        for layer in reversed(range(layers)):
            for (_,c,group),action in reversed(list(zip(terms,actions))):
                post=history[index];index-=1
                gradient[2*layer+group]+=2*np.vdot(adj,-1j*c*apply_pauli(post,action)).real/scale
                angle=c*parameters[2*layer+group]
                adj=np.cos(angle)*adj+1j*np.sin(angle)*apply_pauli(adj,action)
        return objective,gradient
    initial=np.full(2*layers,dt/layers)
    fit=minimize(loss,initial,jac=True,method='L-BFGS-B',bounds=[(0,2*dt)]*len(initial),
                 options={'maxiter':150,'ftol':1e-13,'gtol':1e-9})
    params=fit.x if fit.fun<=loss(initial)[0] else initial
    return evolution_circuit(labels,coeffs,params),params,{'layers':layers,'probe_columns':scale,
         'training_error':float(loss(params)[0]),'initial_training_error':float(loss(initial)[0]),
         'optimizer_success':bool(fit.success)}


def compile_basic(qc):
    return transpile(qc,basis_gates=['rz','sx','x','cx'],optimization_level=2,seed_transpiler=11)


def metrics(qc):
    return {'qubits':qc.num_qubits,'depth':qc.depth(),'cx_count':int(qc.count_ops().get('cx',0)),
            'operations':dict(qc.count_ops())}
