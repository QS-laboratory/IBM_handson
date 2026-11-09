"""Joint compression, explicit suppression, and persisted QPU workflow without QPU access."""
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CORE = ROOT / "tools/h4_joint_hardware.py"
BOOK = ROOT / "08_h4_joint_compression_suppression.ipynb"


def source(start):
    book = json.loads(BOOK.read_text())
    sources = ("".join(c["source"]) for c in book["cells"] if c["cell_type"] == "code")
    return next(s for s in sources if s.startswith(start))


@pytest.fixture(scope="module")
def namespace():
    ns = {}
    with pytest.MonkeyPatch.context() as patch:
        patch.setenv("H4_DATA_DIR", str(ROOT / "data/h4_pec_v1"))
        exec(source("H4_DATA_COMMIT ="), ns)
    exec(CORE.read_text(), ns)
    exec(source("def wilson_interval"), ns)
    return ns


def test_notebook_embeds_current_helpers():
    assert source('"""Exercise 08 helpers').strip() == CORE.read_text().strip()


def test_joint_mps_preserves_complex_state_and_qubit_order(namespace):
    rng = np.random.default_rng(12)
    target = rng.normal(size=512)+1j*rng.normal(size=512)
    target /= np.linalg.norm(target)
    cores, rebuilt, retained = namespace["joint_mps"](target, 16)
    np.testing.assert_allclose(rebuilt, target, atol=1e-12)
    assert abs(retained-1) < 1e-12
    assert max(c.shape[2] for c in cores) <= 16
    _, truncated, retained = namespace["joint_mps"](target, 2)
    assert 0 < retained < 1
    assert abs(np.linalg.norm(truncated)-1) < 1e-12


def test_two_qubit_application_matches_qiskit(namespace):
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Statevector, random_unitary
    state = namespace["initial_joint_state"]()
    gate = random_unitary(4, seed=9).data
    for pair in [(8, 0), (0, 1), (3, 4)]:
        circuit = QuantumCircuit(9)
        circuit.unitary(gate, list(pair))
        expected = Statevector(state).evolve(circuit).data
        np.testing.assert_allclose(namespace["apply_two"](state, gate, pair), expected, atol=1e-12)


@pytest.fixture(scope="module")
def trained(namespace):
    target = namespace["Statevector"].from_instruction(namespace["A_APPROX"]).data
    joint, report = namespace["compress_joint"](target, layers=8, sweeps=600, seeds=(3,))
    return joint, report


def test_svd_sweeps_improve_joint_fidelity_without_controlled_corrections(namespace, trained):
    joint, report = trained
    assert report["fidelity"] > .99
    assert np.min(np.diff(report["runs"][0]["history"])) > -1e-12
    assert report["two_qubit_blocks"] == 32
    assert max(i.operation.num_qubits for i in joint.data) == 2
    beta = np.array([-np.pi/2, 0., np.pi/2])
    logical = namespace["logical_phase_circuit"](joint, 1)
    for angle, expected in zip(beta, namespace["ideal_signal"](joint, 1, beta)):
        circuit = logical.assign_parameters({next(iter(logical.parameters)): angle})
        circuit.remove_final_measurements(inplace=True)
        actual = abs(namespace["Statevector"].from_instruction(circuit).data[0])**2
        np.testing.assert_allclose(actual, expected, atol=1e-11)


@pytest.mark.parametrize("gate", ["cx", "cz", "ecr"])
def test_explicit_pauli_twirls_preserve_full_unitary(namespace, gate):
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Operator
    circuit = QuantumCircuit(2)
    getattr(circuit, gate)(0, 1)
    for seed in range(16):
        twirled = namespace["twirl_native"](circuit, seed)
        np.testing.assert_allclose(Operator(twirled).data, Operator(circuit).data, atol=1e-11)


@pytest.fixture(scope="module")
def preview(namespace, trained):
    provider = pytest.importorskip("qiskit_ibm_runtime.fake_provider")
    backend = provider.FakeKawasaki()
    config = dict(L=1, phases=np.linspace(-np.pi, np.pi, 4, endpoint=False).tolist(),
                  seed=42, layout_candidates=2, randomizations=1, shots_per_phase=32,
                  minimum_joint_fidelity=.99, max_active_qubits=11)
    comparison = namespace["compile_comparison"](backend, *trained, config)
    return backend, comparison, config


def test_calibrated_compilation_and_dd_are_explicit_and_ideal_checked(preview):
    backend, comparison, config = preview
    reports = comparison["specification"]["reports"]
    assert len(comparison["circuits"]) == 12
    assert reports[0]["initial_layout"] == reports[1]["initial_layout"]
    assert reports[1]["two_qubit_gates"] < reports[0]["two_qubit_gates"]
    assert reports[2]["x_gates"] > reports[1]["x_gates"]
    assert all(r["ideal_probability_max_error"] < 1e-7 for r in reports)
    assert all(len(r["active_physical_qubits"]) == 9 for r in reports)
    # Unused device wires must not acquire delays or DD gates.
    for circuit in comparison["circuits"]:
        active = {circuit.find_bit(q).index for item in circuit.data
                  if item.operation.name != "delay" for q in item.qubits}
        assert all(circuit.find_bit(i.qubits[0]).index in active for i in circuit.data
                   if i.operation.name == "delay")


def test_fake_noise_includes_explicit_idle_relaxation(namespace, preview, monkeypatch):
    pytest.importorskip("qiskit_aer")
    import qiskit_aer.noise
    backend, comparison, config = preview
    real_pass = qiskit_aer.noise.RelaxationNoisePass
    captured = []

    def relaxation(**kwargs):
        captured.append(kwargs)
        return real_pass(**kwargs)

    monkeypatch.setattr(qiskit_aer.noise, "RelaxationNoisePass", relaxation)
    record = namespace["run_fake_comparison"](comparison, backend)
    assert captured and captured[0]["op_types"] == [namespace["Delay"]]
    assert len(captured[0]["t1s"]) == backend.num_qubits
    assert all(sum(c.values()) == 32 for c in record["counts"])
    scans = namespace["aggregate_comparison"](record)
    assert [s["arm"] for s in scans] == ["paired", "joint", "suppressed"]
    assert all(sum(s["actual_shots"]) == 128 for s in scans)
    record["specification"] = dict(record["specification"])
    record["specification"]["schedule"] = record["specification"]["schedule"][:-1]
    with pytest.raises(ValueError, match="Incomplete"):
        namespace["aggregate_comparison"](record)


def test_real_workflow_default_never_authenticates_or_submits(monkeypatch):
    import builtins
    monkeypatch.setattr(builtins, "input", lambda *a: pytest.fail("Unexpected prompt"))
    ns = {}
    exec(source("USE_IBM ="), ns)
    assert ns["IBM_SERVICE"] is None
    exec(source("RUN_HARDWARE ="), ns)
    assert ns["real_record"] is None


def test_real_job_persists_and_reuses_exact_precompiled_circuits(namespace, preview, tmp_path, monkeypatch):
    import qiskit_ibm_runtime
    backend, comparison, config = preview
    raw = [{"000000000": 1, "000000001": 31} for _ in comparison["circuits"]]
    calls = []
    job = SimpleNamespace(job_id=lambda: "mock-joint-job", result=lambda: [
        SimpleNamespace(data=SimpleNamespace(meas=SimpleNamespace(get_counts=lambda c=c: c)))
        for c in raw])

    class Sampler:
        def __init__(self, mode, options):
            assert mode is backend
            assert options["dynamical_decoupling"]["enable"] is False
            assert options["twirling"]["enable_gates"] is False

        def run(self, circuits, shots):
            assert circuits is comparison["circuits"] and shots == 32
            calls.append("submit")
            return job

    monkeypatch.setattr(qiskit_ibm_runtime, "SamplerV2", Sampler)
    service = SimpleNamespace(job=lambda job_id: job)
    acquire = namespace["acquire_real_comparison"]
    args = (comparison, service, backend, "test-instance", tmp_path, config)
    record, path = acquire(*args)
    assert json.loads(path.read_text())["counts"] == raw
    with path.with_suffix(".qpy").open("rb") as stream:
        assert len(namespace["qpy"].load(stream)) == 12
    acquire(*args)
    assert calls == ["submit"]
    pending = dict(record)
    pending.pop("counts")
    path.write_text(json.dumps(pending))
    acquire(*args)
    assert calls == ["submit"]
    pending.pop("job_id")
    path.write_text(json.dumps(pending))
    with pytest.raises(RuntimeError, match="Interrupted submission"):
        acquire(*args)
    with pytest.raises(ValueError, match="Settings changed"):
        acquire(comparison, service, backend, "test-instance", tmp_path, dict(config, L=4))


def test_compression_quality_guard_runs_before_compilation(namespace, trained):
    joint, report = trained
    with pytest.raises(ValueError, match="fidelity"):
        namespace["compile_comparison"](None, joint, dict(report, fidelity=.5),
                                        {"minimum_joint_fidelity": .99})
