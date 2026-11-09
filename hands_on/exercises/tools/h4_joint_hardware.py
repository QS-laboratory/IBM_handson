"""Exercise 08 helpers, embedded verbatim in the notebook for standalone Colab use.

The pinned H4 bootstrap supplies np, QuantumCircuit, Statevector, qpy, transpile,
A_APPROX, U_APPROX, EVOLUTION_APPROX and the selected model. No account access
or job submission occurs on loading this file.
"""
import hashlib
import json
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
from qiskit import QuantumCircuit, qpy, transpile
from qiskit.circuit import Parameter, Delay
from qiskit.circuit.library import XGate
from qiskit.quantum_info import Statevector, Operator, Pauli
from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import (
    ALAPScheduleAnalysis, ConstrainedReschedule, PadDelay, PadDynamicalDecoupling,
)

CHAIN = (8, 0, 1, 2, 3, 4, 5, 6, 7)
ARMS = ("paired", "joint", "suppressed")


def digest_json(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def joint_mps(target, chi=16):
    """TT-SVD of the nine-qubit joint state in ancilla-first chain order."""
    target = np.asarray(target, complex)
    if target.shape != (512,) or not np.isclose(np.linalg.norm(target), 1):
        raise ValueError("A normalized nine-qubit joint target is required")
    if not isinstance(chi, int) or isinstance(chi, bool) or chi < 1:
        raise ValueError("chi must be a positive integer")
    axes = [8-q for q in CHAIN]
    tensor = target.reshape([2]*9).transpose(axes)
    cores, rank = [], 1
    for _ in range(8):
        u, s, vh = np.linalg.svd(tensor.reshape(rank*2, -1), full_matrices=False)
        keep = min(chi, len(s))
        cores.append(u[:, :keep].reshape(rank, 2, keep))
        tensor = s[:keep, None]*vh[:keep]
        rank = keep
    cores.append(tensor.reshape(rank, 2, 1))
    rebuilt = cores[0][0]
    for core in cores[1:]:
        rebuilt = np.tensordot(rebuilt, core, axes=(-1, 0))
    state = rebuilt[..., 0].transpose(np.argsort(axes)).reshape(-1)
    retained = float(np.linalg.norm(state)**2)
    state /= np.sqrt(retained)
    return cores, state, retained


def two_site_view(state, pair):
    # Qiskit's two-qubit matrix basis is |q_second q_first>.
    a, b = pair
    axes = [8-b, 8-a]+[i for i in range(9) if i not in (8-b, 8-a)]
    return state.reshape([2]*9).transpose(axes).reshape(4, -1), axes


def apply_two(state, gate, pair):
    matrix, axes = two_site_view(state, pair)
    return (gate@matrix).reshape([2]*9).transpose(np.argsort(axes)).reshape(-1)


def initial_joint_state():
    # A product determinant times |+>_ancilla, prepared without controlled gates.
    state = np.zeros(512, complex)
    state[165] = state[165+256] = 1/np.sqrt(2)
    return state


def compress_joint(target, layers=8, sweeps=2000, seeds=(3, 17, 42), chi=16):
    if layers < 1 or sweeps < 1 or not seeds:
        raise ValueError("Positive layers/sweeps and at least one seed are required")
    cores, mps_target, retained = joint_mps(target, chi)
    pairs = [(CHAIN[j], CHAIN[j+1]) for layer in range(layers)
             for j in range(layer % 2, 8, 2)]
    best, runs = None, []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        gates = []
        for _ in pairs:
            z = np.eye(4)+.05*(rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4)))
            u, _, vh = np.linalg.svd(z)
            gates.append(u@vh)
        history = []
        for sweep in range(sweeps):
            if sweep % 2 == 0:
                suffix = [None]*(len(gates)+1)
                suffix[-1] = mps_target
                for j in reversed(range(len(gates))):
                    suffix[j] = apply_two(suffix[j+1], gates[j].conj().T, pairs[j])
                state = initial_joint_state()
                for j, pair in enumerate(pairs):
                    x, _ = two_site_view(state, pair)
                    y, _ = two_site_view(suffix[j+1], pair)
                    u, _, vh = np.linalg.svd(y@x.conj().T)
                    gates[j] = u@vh
                    state = apply_two(state, gates[j], pair)
            else:
                prefix = [initial_joint_state()]
                for gate, pair in zip(gates, pairs):
                    prefix.append(apply_two(prefix[-1], gate, pair))
                state = mps_target.copy()
                for j in reversed(range(len(gates))):
                    x, _ = two_site_view(prefix[j], pairs[j])
                    y, _ = two_site_view(state, pairs[j])
                    u, _, vh = np.linalg.svd(y@x.conj().T)
                    gates[j] = u@vh
                    state = apply_two(state, gates[j].conj().T, pairs[j])
            prepared = initial_joint_state()
            for gate, pair in zip(gates, pairs):
                prepared = apply_two(prepared, gate, pair)
            history.append(float(abs(np.vdot(mps_target, prepared))**2))
        fidelity = float(abs(np.vdot(target, prepared))**2)
        runs.append(dict(seed=seed, fidelity=fidelity, history=history))
        print(f"Joint compression: seed={seed}, fidelity={fidelity:.8f}", flush=True)
        if best is None or fidelity > best["fidelity"]:
            best = dict(fidelity=fidelity, gates=np.array(gates), seed=seed)
    circuit = QuantumCircuit(9, name="joint_MPS_compressed")
    for q in (0, 2, 5, 7):
        circuit.x(q)
    circuit.h(8)
    for gate, pair in zip(best["gates"], pairs):
        circuit.unitary(gate, list(pair))
    actual = Statevector.from_instruction(circuit).data
    np.testing.assert_allclose(abs(np.vdot(target, actual))**2, best["fidelity"], atol=1e-11)
    fingerprint = hashlib.sha256(best["gates"].tobytes()+np.array(pairs).tobytes()).hexdigest()
    report = dict(layers=layers, sweeps=sweeps, seeds=list(seeds), selected_seed=best["seed"],
                  chi=chi, bond_dimensions=[c.shape[2] for c in cores[:-1]],
                  mps_retained_norm_squared=retained, fidelity=best["fidelity"],
                  circuit_digest=fingerprint, two_qubit_blocks=len(pairs), runs=runs)
    return circuit, report


def logical_phase_circuit(prep, length):
    beta = Parameter("beta")
    circuit = QuantumCircuit(9)
    circuit.compose(prep, inplace=True)
    for _ in range(length):
        circuit.compose(EVOLUTION_APPROX, range(8), inplace=True)
    circuit.p(beta, 8)
    circuit.compose(prep.inverse(), inplace=True)
    circuit.measure_all()
    return circuit


def ideal_signal(prep, length, beta):
    phi = Statevector.from_instruction(prep).data
    evolution = np.linalg.matrix_power(U_APPROX, length)
    a, b = [np.vdot(branch, evolution@branch) for branch in [phi[:256], phi[256:]]]
    return np.clip(abs(a+np.exp(1j*np.asarray(beta))*b)**2, 0, 1)


def instruction_cost(circuit, backend):
    """Calibration-based heuristic, not an estimated circuit success probability."""
    score, missing = 0., 0
    for item in circuit.data:
        if item.operation.name in ("barrier", "delay"):
            continue
        qubits = tuple(circuit.find_bit(q).index for q in item.qubits)
        entries = backend.target.get(item.operation.name, {})
        prop = entries.get(qubits)
        error = getattr(prop, "error", None)
        if error is None:
            error = 0. if item.operation.name == "rz" else .01
            missing += 1
        score -= np.log1p(-min(max(float(error), 0), .999999))
    return float(score), missing


def calibrated_paths(backend, count=4, beam_width=256):
    edges = {}
    neighbors = {q: set() for q in range(backend.num_qubits)}
    for name in backend.target.operation_names:
        operation = backend.target.operation_from_name(name)
        if operation.num_qubits != 2:
            continue
        for pair, prop in backend.target[name].items():
            if pair is None:
                continue
            a, b = pair
            error = .01 if prop is None or prop.error is None else prop.error
            key = tuple(sorted(pair))
            edges[key] = min(edges.get(key, 1.), -np.log1p(-min(error, .999999)))
            neighbors[a].add(b)
            neighbors[b].add(a)
    readout = {}
    for q in neighbors:
        prop = backend.target.get("measure", {}).get((q,))
        error = .01 if prop is None or prop.error is None else prop.error
        readout[q] = -np.log1p(-min(error, .999999))
    beam = [(readout[q], (q,)) for q in neighbors if neighbors[q]]
    for _ in range(8):
        expanded = [(cost+edges[tuple(sorted((path[-1], q)))]+readout[q], path+(q,))
                    for cost, path in beam for q in sorted(neighbors[path[-1]]) if q not in path]
        beam = sorted(expanded)[:beam_width]
        if not beam:
            raise ValueError("No nine-qubit path found in this target")
    paths, seen = [], set()
    for _, path in beam:
        canonical = min(path, path[::-1])
        if canonical not in seen:
            paths.append(path)
            seen.add(canonical)
        if len(paths) >= count:
            break
    return paths


def chain_layout(path):
    layout = [None]*9
    for logical, physical in zip(CHAIN, path):
        layout[logical] = physical
    return layout


def active_qubits(circuit):
    return sorted({circuit.find_bit(q).index for item in circuit.data
                   if item.operation.name not in ("barrier", "delay") for q in item.qubits})


def twirl_native(circuit, seed):
    """Explicit exact Pauli twirls for native CX/CZ/ECR, including global phases."""
    rng = np.random.default_rng(seed)
    labels = [a+b for a in "IXYZ" for b in "IXYZ"]
    paulis = [Pauli(label).to_matrix() for label in labels]
    cache = {}
    out = circuit.copy_empty_like()
    out.global_phase = circuit.global_phase

    def append_pauli(label, qubits):
        for letter, q in zip(label[::-1], qubits):
            if letter != "I":
                getattr(out, letter.lower())(q)

    for item in circuit.data:
        op = item.operation
        if op.num_qubits != 2 or op.name == "barrier":
            out.append(op, item.qubits, item.clbits)
            continue
        if op.name not in ("cx", "cz", "ecr"):
            raise ValueError(f"Explicit Pauli twirling supports CX/CZ/ECR, received {op.name}")
        if op.name not in cache:
            u = Operator(op).data
            table = []
            for pre in paulis:
                conjugated = u@pre@u.conj().T
                overlaps = np.array([np.trace(p.conj().T@conjugated)/4 for p in paulis])
                index = int(np.argmax(abs(overlaps)))
                if not np.isclose(abs(overlaps[index]), 1., atol=1e-10):
                    raise ValueError("Native gate is not Clifford")
                table.append((labels[index], float(np.angle(overlaps[index]))))
            cache[op.name] = table
        index = int(rng.integers(16))
        post, phase = cache[op.name][index]
        append_pauli(labels[index], item.qubits)
        out.append(op, item.qubits, item.clbits)
        append_pauli(post, item.qubits)
        out.global_phase += phase
    return out


def schedule_native(circuit, backend, dd=False):
    active = active_qubits(circuit)
    target = backend.target
    unbarred = circuit.copy_empty_like()
    unbarred.global_phase = circuit.global_phase
    for item in circuit.data:
        if item.operation.name != "barrier":
            unbarred.append(item.operation, item.qubits, item.clbits)
    padding = (PadDynamicalDecoupling(target=target, dd_sequence=[XGate(), XGate()],
               qubits=active, pulse_alignment=target.pulse_alignment, skip_reset_qubits=False)
               if dd else PadDelay(target=target))
    scheduled = PassManager([ALAPScheduleAnalysis(target=target),
                            ConstrainedReschedule(target=target), padding]).run(unbarred)
    # Padding unused device wires would make a local density simulation enormous.
    # They never interact with the active circuit and are excluded in both modes.
    clean = scheduled.copy_empty_like()
    clean.global_phase = scheduled.global_phase
    for item in scheduled.data:
        if item.operation.name == "delay" and scheduled.find_bit(item.qubits[0]).index not in active:
            continue
        clean.append(item.operation, item.qubits, item.clbits)
    return clean


def circuit_signature(circuit):
    return digest_json(dict(qubits=circuit.num_qubits, clbits=circuit.num_clbits,
        phase=str(circuit.global_phase), instructions=[
        [i.operation.name, [str(p) for p in i.operation.params],
         getattr(i.operation, "unit", None),
         [circuit.find_bit(q).index for q in i.qubits],
         [circuit.find_bit(c).index for c in i.clbits]] for i in circuit.data]))


def compact_unitary_circuit(circuit):
    active = active_qubits(circuit)
    mapping = {q: j for j, q in enumerate(active)}
    compact = QuantumCircuit(len(active))
    compact.global_phase = circuit.global_phase
    for item in circuit.data:
        if item.operation.name not in ("measure", "barrier", "delay"):
            compact.append(item.operation, [mapping[circuit.find_bit(q).index] for q in item.qubits])
    return compact


def compile_comparison(backend, joint, compression, config):
    if compression["fidelity"] < config["minimum_joint_fidelity"]:
        raise ValueError("Joint compression fidelity is below the configured submission threshold")
    paired = transpile(logical_phase_circuit(A_APPROX, config["L"]), backend=backend,
                       optimization_level=1, seed_transpiler=7, num_processes=1)
    standard_layout = paired.layout.initial_index_layout(filter_ancillas=True)
    simple_joint = transpile(logical_phase_circuit(joint, config["L"]), backend=backend,
                            initial_layout=standard_layout, optimization_level=1,
                            seed_transpiler=7, approximation_degree=1., num_processes=1)
    candidates = []
    for path in calibrated_paths(backend, count=config["layout_candidates"]):
        compiled = transpile(logical_phase_circuit(joint, config["L"]), backend=backend,
                             initial_layout=chain_layout(path), optimization_level=3,
                             seed_transpiler=7, approximation_degree=1., num_processes=1)
        cost, missing = instruction_cost(compiled, backend)
        candidates.append((cost, missing, compiled))
    cost, missing, aware = min(candidates, key=lambda item: item[0])
    templates = dict(paired=paired, joint=simple_joint, suppressed=aware)
    ideals = {arm: ideal_signal(A_APPROX if arm == "paired" else joint,
                               config["L"], config["phases"]) for arm in ARMS}
    circuits, schedule, reports = [], [], []
    per_circuit = config["shots_per_phase"]//config["randomizations"]
    for arm in ARMS:
        ideal_errors = []
        for repeat in range(config["randomizations"]):
            for j, beta in enumerate(config["phases"]):
                template = templates[arm]
                bound = template.assign_parameters({next(iter(template.parameters)): beta})
                if arm == "suppressed":
                    # Reuse each random Pauli pattern across the entire phase scan.
                    bound = twirl_native(bound, config["seed"]+1000*repeat)
                    # Lower Pauli gates into the native basis without moving physical wires.
                    bound = transpile(bound, backend=backend, optimization_level=0,
                        initial_layout=list(range(backend.num_qubits)), routing_method="none",
                        seed_transpiler=7, num_processes=1)
                scheduled = schedule_native(bound, backend, dd=(arm == "suppressed"))
                if len(active_qubits(scheduled)) > config["max_active_qubits"]:
                    raise ValueError("Routing uses too many qubits for the configured local simulation limit")
                # Check every actual native circuit before it becomes eligible for submission.
                actual = abs(Statevector.from_instruction(compact_unitary_circuit(scheduled)).data[0])**2
                expected = float(ideals[arm][j])
                ideal_errors.append(float(abs(actual-expected)))
                if not np.isclose(actual, expected, atol=1e-7, rtol=0):
                    raise ValueError(f"Ideal circuit check failed: {arm}, beta={beta}, "
                                     f"actual={actual}, expected={expected}")
                circuits.append(scheduled)
                schedule.append(dict(arm=arm, repeat=repeat, phase_index=j, beta=beta))
        representative = circuits[-len(config["phases"])*config["randomizations"]]
        score, missing = instruction_cost(representative, backend)
        report = dict(arm=arm, depth=representative.depth(),
            two_qubit_gates=sum(i.operation.num_qubits == 2 for i in representative.data),
            x_gates=representative.count_ops().get("x", 0),
            delay_instructions=representative.count_ops().get("delay", 0),
            active_physical_qubits=active_qubits(representative),
            initial_layout=templates[arm].layout.initial_index_layout(filter_ancillas=True),
            calibration_cost=score, missing_calibration_entries=missing)
        report["ideal_probability_max_error"] = max(ideal_errors)
        reports.append(report)
        print(json.dumps(report), flush=True)
    order = np.random.default_rng(config["seed"]).permutation(len(circuits))
    circuits = [circuits[int(j)] for j in order]
    schedule = [schedule[int(j)] for j in order]
    properties = backend.properties()
    dense_phase = np.linspace(-np.pi, np.pi, 401)
    ideal_curves = {arm: ideal_signal(A_APPROX if arm == "paired" else joint,
                    config["L"], dense_phase).tolist() for arm in ARMS}
    used = sorted({q for report in reports for q in report["active_physical_qubits"]})
    qubit_calibration = [dict(qubit=q, t1=backend.target.qubit_properties[q].t1,
                             t2=backend.target.qubit_properties[q].t2) for q in used]
    spec = dict(configuration=config, backend=backend.name, num_qubits=backend.num_qubits,
        calibration_last_update=properties.last_update_date.isoformat(),
        compression_digest=compression["circuit_digest"], compression_fidelity=compression["fidelity"],
        reports=reports, schedule=schedule, shots_per_circuit=per_circuit,
        ideal_phase_grid=dense_phase.tolist(), ideal_curves=ideal_curves,
        qubit_calibration=qubit_calibration,
        isa_signatures=[circuit_signature(c) for c in circuits])
    key = digest_json(spec)
    print(f"Preview: {len(circuits)} circuits × {per_circuit} shots = {len(circuits)*per_circuit}")
    return dict(specification=spec, circuits=circuits, key=key)


def run_fake_comparison(preview, backend):
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import RelaxationNoisePass
    simulator = AerSimulator.from_backend(backend, method="density_matrix",
        max_parallel_threads=4, max_parallel_experiments=1, enable_truncation=True)
    if simulator.options.noise_model.is_ideal():
        raise ValueError("A nonideal fake backend noise model is required")
    props = backend.target.qubit_properties
    t1s = [p.t1 for p in props]
    t2s = [min(p.t2, 2*p.t1) for p in props]
    idle_noise = PassManager([RelaxationNoisePass(t1s=t1s, t2s=t2s,
                                                dt=backend.dt, op_types=[Delay])])
    raw = []
    for j, circuit in enumerate(preview["circuits"]):
        noisy = idle_noise.run(circuit)
        result = simulator.run(noisy, shots=preview["specification"]["shots_per_circuit"],
            seed_simulator=preview["specification"]["configuration"]["seed"]+j).result()
        if not result.success:
            raise RuntimeError(result.status)
        raw.append(dict(result.get_counts()))
        if (j+1) % 12 == 0:
            print(f"Local noisy simulation: {j+1}/{len(preview['circuits'])}", flush=True)
    return dict(source="fake_raw_counts", specification=preview["specification"], counts=raw,
        noise_description="Snapshot gate/readout noise plus T1/T2 relaxation on explicit delays; no added coherent drift")


def atomic_json(path, record):
    import os
    import tempfile
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False) as stream:
        temp = Path(stream.name)
        json.dump(record, stream, indent=2)
    try:
        os.replace(temp, path)
    finally:
        temp.unlink(missing_ok=True)


def acquire_real_comparison(preview, service, backend, instance, directory, config):
    if preview is None or service is None or backend is None:
        raise RuntimeError("Connect and inspect the real-backend preview before submission")
    spec = preview["specification"]
    if spec["configuration"] != config or spec["backend"] != backend.name:
        raise ValueError("Settings changed; regenerate the preview")
    spec = dict(spec, instance_sha256=hashlib.sha256(instance.encode()).hexdigest())
    key = digest_json(spec)
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory/(key+".json")
    if path.exists():
        record = json.loads(path.read_text())
        if record["specification"] != spec:
            raise ValueError("Stored specification differs")
        if "counts" in record:
            return record, path
        if "job_id" not in record:
            raise RuntimeError(f"Interrupted submission marker: {path}; inspect IBM jobs before recovery")
        job = service.job(record["job_id"])
    else:
        from qiskit_ibm_runtime import SamplerV2
        qpy_path = path.with_suffix(".qpy")
        with qpy_path.open("wb") as stream:
            qpy.dump(preview["circuits"], stream)
        record = dict(source="ibm_submission", specification=spec, isa_qpy=qpy_path.name,
                      created_at_utc=datetime.now(timezone.utc).isoformat())
        with path.open("x") as stream:
            json.dump(record, stream, indent=2)
        # Suppression is explicit in the submitted circuits. Avoid a second hidden layer.
        options = dict(dynamical_decoupling=dict(enable=False),
                       twirling=dict(enable_gates=False, enable_measure=False))
        job = SamplerV2(mode=backend, options=options).run(
            preview["circuits"], shots=spec["shots_per_circuit"])
        record["job_id"] = job.job_id()
        atomic_json(path, record)
        print("Saved job ID:", record["job_id"], flush=True)
    record.update(source="ibm_raw_counts", counts=[dict(pub.data.meas.get_counts()) for pub in job.result()],
                  retrieved_at_utc=datetime.now(timezone.utc).isoformat())
    atomic_json(path, record)
    return record, path


def aggregate_comparison(record):
    spec = record["specification"]
    config = spec["configuration"]
    schedule, raw = spec["schedule"], record["counts"]
    if len(raw) != len(schedule):
        raise ValueError("Incomplete counts")
    expected = {(arm, repeat, j) for arm in ARMS for repeat in range(config["randomizations"])
                for j in range(len(config["phases"]))}
    seen, totals = set(), {arm: [np.zeros(len(config["phases"]), int),
                                 np.zeros(len(config["phases"]), int)] for arm in ARMS}
    for row, counts in zip(schedule, raw):
        key = (row["arm"], row["repeat"], row["phase_index"])
        if key not in expected or key in seen or not np.isclose(row["beta"], config["phases"][key[2]], atol=1e-12):
            raise ValueError("Duplicate, missing, or mismatched PUB")
        seen.add(key)
        if not counts or any(not isinstance(s, str) or len(s) != 9 or set(s)-set("01")
            or not isinstance(v, int) or isinstance(v, bool) or v < 0 for s, v in counts.items()):
            raise ValueError("Invalid nine-bit counts")
        n = sum(counts.values())
        if n <= 0:
            raise ValueError("No actual shots")
        totals[row["arm"]][0][key[2]] += counts.get("0"*9, 0)
        totals[row["arm"]][1][key[2]] += n
    if seen != expected:
        raise ValueError("Incomplete schedule")
    # summarize_scan is the same conditional phase-dependence test used in 07.
    scans = [dict(arm=arm, **summarize_scan(config["phases"], *totals[arm], seed=42+j))
             for j, arm in enumerate(ARMS)]
    for scan in scans:
        scan["p_adjusted"] = min(1., len(scans)*scan["p_value"])
    return scans
