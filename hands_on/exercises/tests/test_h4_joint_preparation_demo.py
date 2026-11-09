"""Compare 06 paired preparation and direct joint preparation without QPU access."""
import ast
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
BOOK = ROOT / "09_h4_joint_preparation_demo.ipynb"


def source(start, book=BOOK):
    data = json.loads(book.read_text())
    sources = ("".join(c["source"]) for c in data["cells"] if c["cell_type"] == "code")
    return next(s for s in sources if s.startswith(start))


def function_tree(text, name):
    return next(n for n in ast.parse(text).body if isinstance(n, ast.FunctionDef) and n.name == name)


def test_reuses_06_bootstrap_and_08_compression_without_suppression():
    assert source("H4_DATA_COMMIT =") == source("H4_DATA_COMMIT =", ROOT/"06_h4_hardware_demo.ipynb")
    helpers = source("CHAIN =")
    original = (ROOT/"tools/h4_joint_hardware.py").read_text()
    for name in ["joint_mps", "two_site_view", "apply_two", "initial_joint_state", "compress_joint"]:
        assert ast.dump(function_tree(helpers, name)) == ast.dump(function_tree(original, name))
    for unwanted in ["twirl_native", "schedule_native", "calibrated_paths"]:
        assert f"def {unwanted}" not in helpers


@pytest.fixture(scope="module")
def namespace(tmp_path_factory):
    ns = {}
    work = tmp_path_factory.mktemp("comparison09")
    with pytest.MonkeyPatch.context() as patch:
        patch.setenv("H4_DATA_DIR", str(ROOT/"data/h4_pec_v1"))
        patch.chdir(work)
        for start in ["H4_DATA_COMMIT =", "CHAIN =", "SHOTS =",
                      "def comparison_settings", "def compile_preparation_comparison"]:
            exec(source(start), ns)
    return ns


def test_l1_defaults_and_paired_preparation_is_unchanged(namespace):
    ns = namespace
    assert ns["L"] == 1 and ns["SHOTS"] == 1024
    assert ns["PHASE_POINTS"] == 12
    np.testing.assert_allclose(ns["BETA"], np.linspace(-np.pi, np.pi, 12, endpoint=False))
    # Full-period coverage separates a constant background from either fringe quadrature.
    design = np.column_stack([np.ones(12), np.cos(ns["BETA"]), np.sin(ns["BETA"])])
    np.testing.assert_allclose(design.T @ design, np.diag([12, 6, 6]), atol=1e-12)
    assert ns["comparison_settings"]()["scan"] == "full_period"
    np.testing.assert_allclose(ns["TRIAL_GAPS"]*ns["TIME"], ns["BETA"])
    assert ns["PREPARATIONS"]["paired"] == ns["A_APPROX"]
    assert ns["PREP_DIGESTS"]["paired"] == ns["BASELINE_DIGEST"]
    assert ns["PREP_DIGESTS"]["joint"] != ns["BASELINE_DIGEST"]
    assert ns["JOINT_REPORT"]["fidelity"] > .99
    assert ns["PREPARATIONS"]["joint"].count_ops()["cx"] < ns["A_PAIRED"].count_ops()["cx"]
    np.testing.assert_allclose(ns["preparation_probability"]("paired", 4, ns["TRIAL_GAPS"]),
                               ns["probability"]("both", 4, ns["TRIAL_GAPS"]), atol=1e-10)


@pytest.mark.parametrize("arm", ["paired", "joint"])
@pytest.mark.parametrize("length", [1, 4])
def test_probability_matches_actual_circuit_at_every_candidate(namespace, arm, length):
    ns = namespace
    gaps = ns["BETA"] / (length * ns["DT"])
    expected = ns["preparation_probability"](arm, length, gaps)
    for gap, probability in zip(gaps, expected):
        circuit = ns["preparation_circuit"](arm, length, float(gap))
        actual = abs(ns["Statevector"].from_instruction(circuit).data[0])**2
        np.testing.assert_allclose(actual, probability, atol=1e-10, rtol=0)


@pytest.fixture(scope="module")
def fake_preview(namespace):
    pytest.importorskip("qiskit_ibm_runtime.fake_provider")
    from qiskit_ibm_runtime.fake_provider import FakeKawasaki
    return namespace["compile_preparation_comparison"](FakeKawasaki())


def test_compiled_arms_share_initial_layout_and_match_ideal(namespace, fake_preview):
    from qiskit import QuantumCircuit
    ns = namespace
    spec = fake_preview["specification"]
    assert len(fake_preview["circuits"]) == 2 * namespace["PHASE_POINTS"]
    assert [s["arm"] for s in spec["schedule"][:4]] == ["paired", "joint", "joint", "paired"]
    assert all(r["initial_layout"] == spec["initial_layout"] for r in spec["resources"])
    for isa, entry in zip(fake_preview["circuits"], spec["schedule"]):
        active = sorted({isa.find_bit(q).index for item in isa.data
                         if item.operation.name not in ("barrier", "measure") for q in item.qubits})
        assert len(active) == 9
        mapping = {q: j for j, q in enumerate(active)}
        compact = QuantumCircuit(len(active))
        compact.global_phase = isa.global_phase
        for item in isa.data:
            if item.operation.name not in ("barrier", "measure"):
                compact.append(item.operation, [mapping[isa.find_bit(q).index] for q in item.qubits])
        actual = abs(ns["Statevector"].from_instruction(compact).data[0])**2
        expected = ns["preparation_probability"](entry["arm"], spec["L"],
                                                [spec["trial_gaps"][entry["point"]]])[0]
        np.testing.assert_allclose(actual, expected, atol=1e-10, rtol=0)


def test_fake_executes_both_arms_and_accounts_for_actual_shots(namespace, fake_preview, monkeypatch):
    pytest.importorskip("qiskit_aer")
    ns = namespace
    monkeypatch.setitem(ns, "SHOTS", 64)
    preview = dict(fake_preview, specification=dict(fake_preview["specification"], shots=64))
    monkeypatch.setitem(ns, "compile_preparation_comparison", lambda backend: preview)
    record = ns["fake_kawasaki_comparison"](64)
    assert len(record["counts"]) == 2 * ns["PHASE_POINTS"] and record["noise_model_is_ideal"] is False
    assert all(sum(counts.values()) == 64 for counts in record["counts"])
    signal = ns["analyze_comparison"](record)
    assert set(signal) == {"paired", "joint"}
    for arm in signal.values():
        assert arm["actual_shots"] == [64]*ns["PHASE_POINTS"]
        assert arm["p_value_adjusted"] == min(1., 2*arm["p_value"])


def test_fake_selection_blocks_ibm_connection_and_submission(monkeypatch):
    import builtins
    import subprocess
    monkeypatch.setattr(builtins, "input", lambda *a: pytest.fail("Unexpected prompt"))
    monkeypatch.setattr(subprocess, "check_call", lambda *a, **k: pytest.fail("Unexpected install"))
    ns = {"USE_FAKE_KAWASAKI": True}
    exec(source("USE_IBM =").replace("USE_IBM = False", "USE_IBM = True", 1), ns)
    assert ns["IBM_SERVICE"] is None
    exec(source("RUN_HARDWARE =").replace("RUN_HARDWARE = False", "RUN_HARDWARE = True", 1), ns)


def test_real_single_job_persists_id_and_counts_and_reuses_cache(namespace, fake_preview, tmp_path, monkeypatch):
    import qiskit_ibm_runtime
    ns = dict(namespace)
    raw = [{"000000000": 3, "000000001": 29} for _ in fake_preview["circuits"]]
    submitted, retrieved = [], []

    class Job:
        def job_id(self):
            return "mock-paired-joint"

        def result(self):
            # Job ID must already be saved before waiting for results.
            path = next((tmp_path/"ibm_jobs").glob("*.json"))
            assert json.loads(path.read_text())["job_id"] == "mock-paired-joint"
            return [SimpleNamespace(data=SimpleNamespace(meas=SimpleNamespace(get_counts=lambda c=c: c))) for c in raw]

    class Sampler:
        def __init__(self, mode, options):
            assert options["twirling"] == {"enable_gates": False, "enable_measure": False}
            assert options["dynamical_decoupling"] == {"enable": False}

        def run(self, circuits, shots):
            assert len(circuits) == len(fake_preview["circuits"]) and shots == 1024
            submitted.append(circuits)
            return Job()

    monkeypatch.setattr(qiskit_ibm_runtime, "SamplerV2", Sampler)
    service = SimpleNamespace(job=lambda job_id: (retrieved.append(job_id), Job())[1])
    backend = SimpleNamespace(name=fake_preview["specification"]["backend"])
    directory = tmp_path/"ibm_jobs"
    record, path = ns["acquire_ibm_comparison"](fake_preview, service, backend, directory=directory)
    assert record["counts"] == raw
    assert path.with_suffix(".qpy").exists()
    assert len(submitted) == 1
    ns["acquire_ibm_comparison"](fake_preview, service, backend, directory=directory)
    assert len(submitted) == 1 and not retrieved
    del record["counts"]
    ns["atomic_comparison_json"](path, record)
    ns["acquire_ibm_comparison"](fake_preview, service, backend, directory=directory)
    assert len(submitted) == 1 and retrieved == ["mock-paired-joint"]

    def plot_failure(*args):
        assert json.loads(path.read_text())["counts"] == raw
        raise RuntimeError("plot failed")

    ns.update(USE_FAKE_KAWASAKI=False, IBM_SERVICE=service, IBM_INSTANCE="", backend=backend,
              preview=fake_preview, plot_comparison=plot_failure)
    ns["acquire_ibm_comparison"] = lambda *a: (json.loads(path.read_text()), path)
    with pytest.raises(RuntimeError, match="plot failed"):
        exec(source("RUN_HARDWARE =").replace("RUN_HARDWARE = False", "RUN_HARDWARE = True", 1), ns)
    saved = json.loads(path.read_text())
    assert saved["counts"] == raw and set(saved["analysis"]) == {"paired", "joint"}


def test_changed_settings_and_interrupted_submission_do_not_submit(namespace, fake_preview, tmp_path):
    ns = dict(namespace)
    service = SimpleNamespace(job=lambda *a: pytest.fail("Unexpected job lookup"))
    backend = SimpleNamespace(name=fake_preview["specification"]["backend"])
    preview = dict(fake_preview, specification=dict(fake_preview["specification"], L=fake_preview["specification"]["L"]+1))
    with pytest.raises(ValueError, match="Settings changed"):
        ns["acquire_ibm_comparison"](preview, service, backend, directory=tmp_path)
    spec = dict(fake_preview["specification"], instance_sha256=ns["hashlib"].sha256(b"").hexdigest())
    key = ns["hashlib"].sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()
    (tmp_path/(key+".json")).write_text(json.dumps(dict(specification=spec)))
    with pytest.raises(RuntimeError, match="Interrupted submission"):
        ns["acquire_ibm_comparison"](fake_preview, service, backend, directory=tmp_path)
