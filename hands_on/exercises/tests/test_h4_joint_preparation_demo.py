"""Check 09's joint-only change to 06, with no QPU access."""
import ast
import json
from pathlib import Path

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


def test_reuses_06_bootstrap_and_fake_execution_recipe():
    old = ROOT / "06_h4_hardware_demo.ipynb"
    assert source("H4_DATA_COMMIT =") == source("H4_DATA_COMMIT =", old)
    for book in [BOOK, old]:
        node = function_tree(source("def fake_kawasaki_scan", book), "fake_kawasaki_scan")
        if book == BOOK:
            actual = ast.dump(node)
        else:
            assert ast.dump(node) == actual


def test_contains_joint_compression_only_without_suppression_helpers():
    helpers = source("CHAIN =")
    original = (ROOT / "tools/h4_joint_hardware.py").read_text()
    for name in ["joint_mps", "two_site_view", "apply_two", "initial_joint_state", "compress_joint"]:
        assert ast.dump(function_tree(helpers, name)) == ast.dump(function_tree(original, name))
    for unwanted in ["twirl_native", "schedule_native", "calibrated_paths", "compile_comparison"]:
        assert f"def {unwanted}" not in helpers


@pytest.fixture(scope="module")
def namespace(tmp_path_factory):
    ns = {}
    work = tmp_path_factory.mktemp("joint09")
    with pytest.MonkeyPatch.context() as patch:
        patch.setenv("H4_DATA_DIR", str(ROOT / "data/h4_pec_v1"))
        patch.chdir(work)
        exec(source("H4_DATA_COMMIT ="), ns)
        exec(source("CHAIN ="), ns)
        exec(source("SHOTS ="), ns)
        exec(source("def fake_kawasaki_scan"), ns)
    return ns


def test_joint_preparation_and_fingerprint_replace_paired_state(namespace):
    ns = namespace
    assert ns["SELECTED_PREP_M"] == "joint8"
    assert ns["SELECTED_EVOLUTION_M"] == 16
    assert ns["CIRCUIT_DIGEST"] != ns["BASELINE_DIGEST"]
    assert ns["JOINT_REPORT"]["fidelity"] > .99
    assert ns["SERIAL_CX"]["preparation"] < 516
    assert ns["L"] == 1 and len(ns["TRIAL_GAPS"]) == 12
    np.testing.assert_allclose(ns["TRIAL_GAPS"]*ns["TIME"], ns["BETA"])


@pytest.mark.parametrize("length", [1, 4])
def test_joint_probability_matches_actual_circuit_at_each_phase(namespace, length):
    ns = namespace
    gaps = ns["BETA"]/(length*ns["DT"])
    probabilities = ns["probability"]("both", length, gaps)
    for gap, expected in zip(gaps, probabilities):
        circuit = ns["experiment_circuit"](length, float(gap))
        actual = abs(ns["Statevector"].from_instruction(circuit).data[0])**2
        np.testing.assert_allclose(actual, expected, atol=1e-10, rtol=0)


def test_real_access_is_disabled_and_fake_selection_blocks_auth(monkeypatch):
    import builtins
    import subprocess
    monkeypatch.setattr(builtins, "input", lambda *a: pytest.fail("Unexpected prompt"))
    monkeypatch.setattr(subprocess, "check_call", lambda *a, **k: pytest.fail("Unexpected install"))
    ns = {"USE_FAKE_KAWASAKI": True}
    exec(source("USE_IBM =").replace("USE_IBM = False", "USE_IBM = True", 1), ns)
    assert ns["IBM_SERVICE"] is None
    exec(source("RUN_HARDWARE =").replace("RUN_HARDWARE = False", "RUN_HARDWARE = True", 1), ns)
    ns.update(USE_FAKE_KAWASAKI=False, backend=None, preview=None)
    with pytest.raises(RuntimeError, match="プレビュー"):
        exec(source("RUN_HARDWARE =").replace("RUN_HARDWARE = False", "RUN_HARDWARE = True", 1), ns)


def test_fake_executes_all_twelve_joint_circuits(namespace, monkeypatch):
    pytest.importorskip("qiskit_aer")
    pytest.importorskip("qiskit_ibm_runtime.fake_provider")
    from qiskit import QuantumCircuit
    captured = []
    transpile = namespace["transpile"]

    def capture(*args, **kwargs):
        circuits = transpile(*args, **kwargs)
        captured.extend(circuits)
        return circuits

    monkeypatch.setitem(namespace, "transpile", capture)
    record = namespace["fake_kawasaki_scan"](64, seed=42)
    assert record["prep_m"] == "joint8"
    assert record["circuit_digest"] == namespace["CIRCUIT_DIGEST"]
    assert len(record["counts"]) == 12
    assert record["actual_shots"] == [64]*12
    assert all(sum(counts.values()) == 64 for counts in record["counts"])
    assert record["noise_model_is_ideal"] is False
    for index in [0, 6, 11]:
        isa = captured[index]
        active = sorted({isa.find_bit(q).index for item in isa.data
                         if item.operation.name not in ("barrier", "measure") for q in item.qubits})
        mapping = {q: j for j, q in enumerate(active)}
        compact = QuantumCircuit(len(active))
        for item in isa.data:
            if item.operation.name not in ("barrier", "measure"):
                compact.append(item.operation, [mapping[isa.find_bit(q).index] for q in item.qubits])
        actual = abs(namespace["Statevector"].from_instruction(compact).data[0])**2
        expected = namespace["probability"]("both", 1, [record["trial_gaps"][index]])[0]
        np.testing.assert_allclose(actual, expected, atol=1e-10, rtol=0)


def test_real_counts_are_preserved_before_signal_plotting(namespace, tmp_path, monkeypatch):
    ns = dict(namespace)
    raw = [{"000000000": 3, "000000001": 29} for _ in range(12)]
    captured = []

    class Measurement:
        def __init__(self, mode, shots, seed, **kwargs):
            assert mode == "ibm" and kwargs["max_hardware_L"] == 1
            assert kwargs["directory"] == "outputs/exercise9/ibm_jobs"
            self.records = [{"counts": raw}]

        def measure(self, case, length, gaps, label):
            assert case == "both" and length == 1 and len(gaps) == 12 and label == "joint_signal"
            return np.full(12, 3/32)

        def save(self, path):
            path = Path(path)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(self.records))
            captured.append(path)

    def plot_failure(*args):
        raise RuntimeError("plot failed")

    monkeypatch.chdir(tmp_path)
    ns.update(USE_FAKE_KAWASAKI=False, backend=object(), IBM_SERVICE=object(),
              preview={"L": 1}, Measurement=Measurement, plot_signal_counts=plot_failure)
    with pytest.raises(RuntimeError, match="plot failed"):
        exec(source("RUN_HARDWARE =").replace("RUN_HARDWARE = False", "RUN_HARDWARE = True", 1), ns)
    assert captured
    assert json.loads(captured[0].read_text())[0]["counts"] == raw
    saved = json.loads((tmp_path/"outputs/exercise9/ibm_signal.json").read_text())
    assert saved["counts"] == raw and saved["signal"]["p_value"] == 1.
