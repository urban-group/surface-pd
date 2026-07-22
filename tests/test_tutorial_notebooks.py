"""Execution and ownership tests for maintained tutorial notebooks."""

from pathlib import Path

import nbformat
from nbclient import NotebookClient

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PHASE_DIAGRAM_NOTEBOOK = (
    PROJECT_ROOT / "examples" / "phase-diagram-python-api.ipynb"
)
ENUMERATION_NOTEBOOK = (
    PROJECT_ROOT / "examples" / "enumeration-python-api.ipynb"
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
    outputs = executed.cells[-1].get("outputs", [])
    assert len(outputs) == 1
    assert outputs[0].output_type == "display_data"


def test_enumeration_notebook_is_a_clean_maintained_tutorial():
    """The enumeration notebook should be portable and source-controlled."""
    notebook = nbformat.read(ENUMERATION_NOTEBOOK, as_version=4)

    assert notebook.cells
    assert notebook.cells[0].cell_type == "markdown"
    assert "Surface-enumeration Python API" in notebook.cells[0].source
    assert all(cell.get("outputs", []) == [] for cell in notebook.cells)
    assert all(cell.get("execution_count") is None for cell in notebook.cells)

    source = ENUMERATION_NOTEBOOK.read_text()
    assert "/Users/" not in source
    assert "EnumerationSlab" in source
    assert "EnumerationSlab.from_file" in source
    assert ".analyze()" in source
    assert "SurfaceEnumerator" in source
    assert "enumlib" in source
    assert "symmetric" in source
    assert "max_cell_size" in source
    assert "DummySpecies" not in source
    assert "POSCAR_Li2O_110.vasp" in source
    assert "symmetric_slab" in source


def test_enumeration_notebook_executes_from_examples_directory():
    """The deterministic tutorial cells should execute beside their inputs."""
    notebook = nbformat.read(ENUMERATION_NOTEBOOK, as_version=4)
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
    outputs = executed.cells[-1].get("outputs", [])
    assert len(outputs) == 1
    assert outputs[0].output_type == "stream"
    assert "enumlib" in outputs[0].text


def test_jupyter_checkpoint_files_are_absent():
    """Ignored editor recovery files should not masquerade as examples."""
    assert not list((PROJECT_ROOT / "examples").rglob(".ipynb_checkpoints"))


def test_enumeration_config_paths_are_relative_to_examples_directory():
    """Maintained enumeration inputs should share one working directory."""
    input_directory = (
        PROJECT_ROOT / "examples" / "enumeration-examples" / "input"
    )

    for input_path in input_directory.glob("*.json"):
        contents = input_path.read_text()
        assert '"target_slab_path": "enumeration-examples/' in contents
        assert "./examples/" not in contents
