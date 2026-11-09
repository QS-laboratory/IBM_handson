"""Check setup without credentials, network calls, installations or QPU jobs."""
import ast
import builtins
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import nbformat
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from ibm_setup import SETUP_CODE


def test_simulator_skips_authentication_and_install(monkeypatch):
    monkeypatch.setattr(subprocess, 'check_call', lambda *a, **k: pytest.fail('Unexpected installation'))
    monkeypatch.setattr(builtins, 'input', lambda *a: pytest.fail('Unexpected prompt'))
    ns = {}
    exec(SETUP_CODE, ns)
    assert ns['IBM_SERVICE'] is None


def test_prepared_connection_reused_in_both_exercises(monkeypatch):
    monkeypatch.setattr(subprocess, 'check_call', lambda *a, **k: None)
    prompts = iter(['crn:test:workshop', 'ibm_test'])
    monkeypatch.setattr(builtins, 'input', lambda *a: next(prompts))
    monkeypatch.setattr('getpass.getpass', lambda *a: 'test-only-token')
    auth = []
    backend = SimpleNamespace(name='ibm_test')

    class Service:
        def __init__(self, **kwargs):
            auth.append(kwargs)
        def backends(self):
            return [backend]
        def backend(self, name):
            assert name == backend.name
            return backend

    monkeypatch.setitem(sys.modules, 'qiskit_ibm_runtime', SimpleNamespace(QiskitRuntimeService=Service))
    ns = {}
    exec(SETUP_CODE.replace('USE_IBM = False', 'USE_IBM = True'), ns)
    assert '_api_key' not in ns
    assert len(auth) == 1
    for path in ROOT.glob('0[45]*ipynb'):
        book = nbformat.read(path, 4)
        source = next(c.source for c in book.cells if c.cell_type == 'code' and c.source.startswith('H4_DATA_COMMIT ='))
        assert 'h4_measurement.py' in source
        runtime_source=(ROOT/'tools/h4_measurement.py').read_text()
        function = next(n for n in ast.parse(runtime_source).body if isinstance(n, ast.FunctionDef) and n.name == 'connect_ibm')
        exec(compile(ast.Module(body=[function], type_ignores=[]), '<connection>', 'exec'), ns)
        service, result = ns['connect_ibm'](ns['IBM_BACKEND_NAME'], ns['IBM_INSTANCE'])
        assert service is ns['IBM_SERVICE'] and result is backend
    assert len(auth) == 1  # No second authentication; no Sampler or job submission.
