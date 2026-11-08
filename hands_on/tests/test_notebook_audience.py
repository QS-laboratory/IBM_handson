"""Check separation in both distributed notebooks and regenerated copies."""
from pathlib import Path
import shutil
import subprocess
import sys
import nbformat

ROOT = Path(__file__).resolve().parents[2]
PARTICIPANT = [
    'hands_on/overview/00_overview.ipynb',
    'hands_on/00_preparation.ipynb',
    'hands_on/01_circuits_and_ibm.ipynb',
    'hands_on/02_inference.ipynb',
    'hands_on/03_mps_mpo_advanced.ipynb',
    'advanced_kanno_h4/H4_author_method.ipynb',
]
INSTRUCTOR_MARKERS = ['講師へのメモ', '講師用：', '講師向け', '講師・上級者', '配布前確認']


def check_participant(path):
    notebook = nbformat.read(path, as_version=4)
    nbformat.validate(notebook)
    assert notebook.metadata.get('audience') == 'participant'
    for cell in notebook.cells:
        if cell.cell_type == 'markdown':
            assert not any(marker in cell.source for marker in INSTRUCTOR_MARKERS)
        if cell.cell_type == 'code':
            assert 'REGENERATE_TENSORS' not in cell.source
            assert '配布時には講師が' not in cell.source


def test_distributed_audiences():
    for name in PARTICIPANT:
        check_participant(ROOT / name)
    teacher = nbformat.read(ROOT / 'advanced_kanno_h4/instructor/H4_regeneration.ipynb', 4)
    assert teacher.metadata.audience == 'instructor'
    assert any('REGENERATE_TENSORS=False' in c.source for c in teacher.cells)


def test_generators_keep_audiences_separate(tmp_path):
    for folder in ['hands_on/tools', 'advanced_kanno_h4/tools', 'versions/v1_2026-11-04']:
        shutil.copytree(ROOT / folder, tmp_path / folder)
    (tmp_path / 'hands_on/overview').mkdir()
    for script in ['hands_on/tools/build_notebooks.py', 'advanced_kanno_h4/tools/build_notebook.py']:
        subprocess.run([sys.executable, str(tmp_path / script)], check=True, capture_output=True)
    for name in PARTICIPANT:
        check_participant(tmp_path / name)
    teacher = nbformat.read(tmp_path / 'advanced_kanno_h4/instructor/H4_regeneration.ipynb', 4)
    assert teacher.metadata.audience == 'instructor'
