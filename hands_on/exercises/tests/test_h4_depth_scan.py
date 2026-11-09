"""Verify the depth experiment without credentials or hardware submissions."""
import ast
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
BOOK = ROOT / "07_h4_hardware_depth_scan.ipynb"


def source(start):
    book = json.loads(BOOK.read_text())
    sources = ("".join(c["source"]) for c in book["cells"] if c["cell_type"] == "code")
    return next(text for text in sources if text.startswith(start))


def definitions(text, namespace):
    tree = ast.parse(text)
    tree.body = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
    exec(compile(tree, str(BOOK), "exec"), namespace)


@pytest.fixture
def stats():
    ns = {"np": np}
    definitions(source("def wilson_interval"), ns)
    ns["CIRCUIT_DIGEST"] = "test-circuit"
    definitions(source("def atomic_json"), ns)
    definitions(source("def aggregate_hardware"), ns)
    return ns


def test_zero_counts_have_nonzero_uncertainty(stats):
    low, high = stats["wilson_interval"](np.array([0, 1024]), np.array([1024, 1024]))
    assert low[0] == 0 and 0.003 < high[0] < 0.004
    assert 0.996 < low[1] < 0.997 and high[1] == 1


@pytest.mark.parametrize("successes", [0, 2, 1024])
def test_flat_scan_does_not_create_a_peak(stats, successes):
    beta = np.linspace(-np.pi, np.pi, 12, endpoint=False)
    scan = stats["summarize_scan"](beta, np.full(12, successes), np.full(12, 1024))
    assert scan["p_value"] == 1
    assert scan["amplitude"] < 1e-12


def test_harmonic_scan_is_detected_with_unequal_shot_counts(stats):
    beta = np.linspace(-np.pi, np.pi, 12, endpoint=False)
    n = np.arange(2000, 3200, 100)
    k = np.rint(n*(0.25+0.2*np.cos(beta))).astype(int)
    scan = stats["summarize_scan"](beta, k, n)
    assert 0.199 < scan["amplitude"] < 0.201
    assert scan["p_value"] * 3 < 0.05


def hardware_record():
    phases = np.linspace(-np.pi, np.pi, 12, endpoint=False)
    config = dict(circuit_digest="test-circuit", phases=phases.tolist(),
                  L_values=[1, 2, 4], repeats=2)
    schedule, counts = [], []
    for repeat in range(2):
        for length in config["L_values"]:
            for j, beta in enumerate(phases):
                schedule.append(dict(repeat=repeat, L=length, phase_index=j, beta=float(beta)))
                counts.append({"000000000": 2, "000000001": 98+repeat})
    return dict(source="ibm_raw_counts", specification=dict(configuration=config, schedule=schedule),
                counts=counts)


def test_aggregation_uses_actual_shots_and_rejects_duplicate_pubs(stats):
    record = hardware_record()
    scans = stats["aggregate_hardware"](record)
    assert [s["L"] for s in scans] == [1, 2, 4]
    assert scans[0]["counts_all_zero"] == [4]*12
    assert scans[0]["actual_shots"] == [201]*12
    np.testing.assert_allclose(scans[0]["observed"], 4/201)
    record["specification"]["schedule"][-1] = record["specification"]["schedule"][0]
    with pytest.raises(ValueError, match="重複"):
        stats["aggregate_hardware"](record)


def test_default_hardware_cells_skip_prompts_and_submission(monkeypatch):
    import builtins
    import subprocess
    monkeypatch.setattr(builtins, "input", lambda *a: pytest.fail("Unexpected credential prompt"))
    monkeypatch.setattr(subprocess, "check_call", lambda *a, **k: pytest.fail("Unexpected install"))
    ns = {}
    exec(source("USE_IBM ="), ns)
    assert ns["IBM_SERVICE"] is None
    exec(source("RUN_HARDWARE ="), ns)
    assert ns["hardware_record"] is None
    ns.update(preview=None, backend=None, OUTPUT_DIR=Path("unused"))
    definitions(source("def atomic_json"), ns)
    with pytest.raises(RuntimeError, match="プレビュー"):
        exec(source("RUN_HARDWARE =").replace("RUN_HARDWARE = False", "RUN_HARDWARE = True"), ns)


@pytest.fixture(scope="module")
def compiled():
    from qiskit.providers.fake_provider import GenericBackendV2
    ns = {}
    with pytest.MonkeyPatch.context() as patch:
        patch.setenv("H4_DATA_DIR", str(ROOT / "data/h4_pec_v1"))
        exec(source("H4_DATA_COMMIT ="), ns)
    exec(source("SHOTS ="), ns)
    exec(source("USE_IBM ="), ns)
    exec(source("from qiskit.circuit import Parameter"), ns)
    exec(source("def atomic_json"), ns)
    # A connected line requires actual routing, rather than an all-to-all mock.
    edges = [(a, b) for a in range(8) for b in [a+1]]
    backend = GenericBackendV2(9, coupling_map=edges+[(b, a) for a, b in edges], seed=42)
    preview = ns["make_preview"](backend, "test-instance")
    return ns, backend, preview


def test_compiled_scan_matches_logical_circuits_on_one_initial_layout(compiled):
    ns, backend, preview = compiled
    reports = preview["specification"]["reports"]
    assert all(r["initial_layout"] == reports[0]["initial_layout"] for r in reports)
    assert len(preview["circuits"]) == 72
    schedule = preview["specification"]["schedule"]
    for length in [1, 2, 4]:
        for phase_index in [3, 6, 9]:
            index = next(i for i, row in enumerate(schedule)
                         if row["L"] == length and row["phase_index"] == phase_index)
            isa = preview["circuits"][index].remove_final_measurements(inplace=False)
            actual = abs(ns["Statevector"].from_instruction(isa).data[0])**2
            gap = schedule[index]["beta"]/(length*ns["DT"])
            expected = ns["probability"]("both", length, [gap])[0]
            np.testing.assert_allclose(actual, expected, atol=1e-10, rtol=0)


def test_submit_persist_reuse_resume_and_guard_changed_settings(compiled, tmp_path, monkeypatch):
    ns, backend, preview = compiled
    import sys
    calls = []
    raw = [{"000000000": 1, "000000001": 19+i} for i in range(72)]
    job = SimpleNamespace(job_id=lambda: "mock-depth-job", result=lambda: [
        SimpleNamespace(data=SimpleNamespace(meas=SimpleNamespace(get_counts=lambda c=c: c)))
        for c in raw])

    class Sampler:
        def __init__(self, mode):
            assert mode is backend

        def run(self, circuits, shots):
            assert circuits is preview["circuits"]
            assert shots == 1024
            calls.append("submit")
            return job

    monkeypatch.setitem(sys.modules, "qiskit_ibm_runtime", SimpleNamespace(SamplerV2=Sampler))
    service = SimpleNamespace(job=lambda job_id: job if job_id == "mock-depth-job" else None)
    args = (preview, service, backend, "test-instance", tmp_path)
    record, path = ns["acquire_scan"](*args)
    assert record["counts"] == raw and json.loads(path.read_text())["job_id"] == "mock-depth-job"
    assert path.with_suffix(".qpy").exists()
    with path.with_suffix(".qpy").open("rb") as stream:
        assert len(ns["qpy"].load(stream)) == 72
    ns["acquire_scan"](*args)
    assert calls == ["submit"]
    pending = dict(record)
    pending.pop("counts")
    path.write_text(json.dumps(pending))
    ns["acquire_scan"](*args)
    assert calls == ["submit"]  # Retrieve the original job, rather than resubmitting.
    pending.pop("job_id")
    path.write_text(json.dumps(pending))
    with pytest.raises(RuntimeError, match="送信中断"):
        ns["acquire_scan"](*args)
    monkeypatch.setitem(ns, "SHOTS", 2048)
    with pytest.raises(RuntimeError, match="設定または接続"):
        ns["acquire_scan"](*args)
