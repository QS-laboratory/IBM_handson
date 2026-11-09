"""Exercise 06's optional local noise simulation, without IBM credentials."""
import ast
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
BOOK = ROOT / "06_h4_hardware_demo.ipynb"


def source_containing(text):
    book = json.loads(BOOK.read_text())
    sources = ("".join(c["source"]) for c in book["cells"] if c["cell_type"] == "code")
    return next(s for s in sources if text in s)


def test_fake_backend_is_optional_and_blocks_real_hardware(monkeypatch):
    import builtins
    import subprocess
    monkeypatch.setattr(builtins, "input", lambda *a: pytest.fail("Unexpected credential prompt"))
    monkeypatch.setattr(subprocess, "check_call", lambda *a, **k: pytest.fail("Unexpected install"))
    ns = {}
    exec(source_containing("USE_FAKE_KAWASAKI ="), ns)
    assert ns["USE_FAKE_KAWASAKI"] is False
    ns["USE_FAKE_KAWASAKI"] = True
    setup = source_containing("USE_IBM =")
    setup = setup.replace("USE_IBM = False", "USE_IBM = True", 1)
    exec(setup, ns)
    assert ns["IBM_SERVICE"] is None
    submit = source_containing("RUN_HARDWARE =")
    submit = submit.replace("RUN_HARDWARE = False", "RUN_HARDWARE = True", 1)
    # With fake mode selected, the real submission branch is skipped entirely.
    exec(submit, ns)


@pytest.mark.parametrize("shots", [0, -1, True, 1.5])
def test_invalid_shots_fail_before_simulator_import(shots):
    ns = {}
    exec(source_containing("def fake_kawasaki_scan"), ns)
    with pytest.raises(ValueError, match="正の整数"):
        ns["fake_kawasaki_scan"](shots)


@pytest.fixture(scope="module")
def noisy_scan():
    pytest.importorskip("qiskit_aer")
    runtime = pytest.importorskip("qiskit_ibm_runtime")
    ns = {}
    with pytest.MonkeyPatch.context() as patch:
        patch.setenv("H4_DATA_DIR", str(ROOT / "data/h4_pec_v1"))
        patch.setattr(runtime, "QiskitRuntimeService", lambda *a, **k: pytest.fail("IBM service accessed"))
        exec(source_containing("H4_DATA_COMMIT ="), ns)
        exec(source_containing("SHOTS = 1024"), ns)
        exec(source_containing("def fake_kawasaki_scan"), ns)
        transpile = ns["transpile"]
        compiled = []

        def capture(*args, **kwargs):
            isa = transpile(*args, **kwargs)
            compiled.extend(isa)
            return isa

        ns["transpile"] = capture
        record = ns["fake_kawasaki_scan"](64, seed=42)
    return ns, record, compiled


def test_local_counts_use_snapshot_noise_and_nine_bit_readout(noisy_scan):
    ns, record, compiled = noisy_scan
    assert record["source"] == "fake_kawasaki_noisy_counts"
    assert record["backend"] == "fake_kawasaki"
    assert record["noise_model_is_ideal"] is False
    assert record["method"] == "density_matrix"
    assert record["snapshot_last_update"]
    assert len(record["snapshot_sha256"]) == 2
    assert len(compiled) == len(record["counts"]) == 9
    assert record["actual_shots"] == [64]*9
    for counts, observed in zip(record["counts"], record["values"]):
        assert sum(counts.values()) == 64
        assert all(len(bits) == 9 and set(bits) <= set("01") for bits in counts)
        assert observed == counts.get("0"*9, 0)/64


def test_fake_backend_routing_preserves_ideal_probabilities(noisy_scan):
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Statevector
    ns, record, compiled = noisy_scan
    for index in [0, 4, 8]:
        isa = compiled[index]
        active = sorted({isa.find_bit(q).index for item in isa.data
                         if item.operation.name not in ["measure", "barrier"] for q in item.qubits})
        mapping = {q: i for i, q in enumerate(active)}
        compact = QuantumCircuit(len(active))
        for item in isa.data:
            if item.operation.name not in ["measure", "barrier"]:
                compact.append(item.operation, [mapping[isa.find_bit(q).index] for q in item.qubits])
        actual = abs(Statevector.from_instruction(compact).data[0])**2
        expected = ns["probability"]("both", 4, [record["trial_gaps"][index]])[0]
        np.testing.assert_allclose(actual, expected, atol=1e-10, rtol=0)


def test_raw_fake_counts_survive_plot_failure(noisy_scan, tmp_path, monkeypatch):
    ns, record, compiled = noisy_scan
    monkeypatch.chdir(tmp_path)
    run_ns = {**ns, "USE_FAKE_KAWASAKI": True, "SHOTS": 64,
              "fake_kawasaki_scan": lambda *a, **k: record}

    def plot_failure(*args):
        raise RuntimeError("plot failed")

    run_ns["plot_fake_kawasaki"] = plot_failure
    with pytest.raises(RuntimeError, match="plot failed"):
        exec(source_containing("fake_kawasaki_record = None"), run_ns)
    path = tmp_path / "outputs/exercise6/fake_kawasaki_measurements.json"
    assert json.loads(path.read_text())["counts"] == record["counts"]
