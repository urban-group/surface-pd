"""Execution and ownership tests for maintained tutorial notebooks."""

from pathlib import Path

import nbformat
from nbclient import NotebookClient

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PHASE_DIAGRAM_NOTEBOOK = (
    PROJECT_ROOT / "examples" / "phase-diagram-python-api.ipynb"
)


def test_phase_diagram_notebook_is_a_clean_maintained_tutorial():
    """The phase-diagram notebook should be portable and source-controlled."""
    notebook = nbformat.read(PHASE_DIAGRAM_NOTEBOOK, as_version=4)

    assert notebook.cells
    assert notebook.cells[0].cell_type == "markdown"
    assert "Phase-diagram Python API" in notebook.cells[0].source
    assert all(cell.get("outputs", []) == [] for cell in notebook.cells)
    assert all(cell.get("execution_count") is None for cell in notebook.cells)
    assert "/Users/" not in PHASE_DIAGRAM_NOTEBOOK.read_text()


def test_phase_diagram_notebook_executes_from_examples_directory():
    """The maintained tutorial should execute beside committed inputs."""
    notebook = nbformat.read(PHASE_DIAGRAM_NOTEBOOK, as_version=4)
    client = NotebookClient(
        notebook,
        timeout=120,
        kernel_name="python3",
        resources={
            "metadata": {"path": str(PROJECT_ROOT / "examples")}
        },
    )

    executed = client.execute()

    assert executed.cells[-1].cell_type == "code"
    assert not executed.cells[-1].get("outputs", [])


def test_jupyter_checkpoint_files_are_absent():
    """Ignored editor recovery files should not masquerade as examples."""
    assert not list((PROJECT_ROOT / "examples").rglob(".ipynb_checkpoints"))
