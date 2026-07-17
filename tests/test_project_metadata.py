"""Tests for project packaging and dependency metadata."""

import tomllib
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
STALE_DOC_COMMANDS = (
    "surface-enumeration.py",
    "surface-pd-plot.py",
    "discharge_pd_gene.py",
    "./scripts/surface-pd-plot.py",
)


def _dependency_names(dependencies):
    """Return normalized package names from dependency specifiers."""
    return {
        dependency.split(">=", maxsplit=1)[0]
        for dependency in dependencies
    }


def test_dev_extra_declares_configured_tools():
    """The dev extra should install the tools configured for local checks."""
    metadata = tomllib.loads((PROJECT_ROOT / "pyproject.toml").read_text())

    dev_dependencies = metadata["project"]["optional-dependencies"]["dev"]

    assert {
        "pytest",
        "pytest-cov",
        "ruff",
        "sphinx",
        "sphinx-rtd-theme",
    } <= _dependency_names(dev_dependencies)


def test_docs_requirements_match_docs_extra():
    """The RTD requirements file should stay aligned with the docs extra."""
    metadata = tomllib.loads((PROJECT_ROOT / "pyproject.toml").read_text())
    docs_extra = metadata["project"]["optional-dependencies"]["docs"]
    docs_requirements = [
        line.strip()
        for line in (PROJECT_ROOT / "docs" / "requirements.txt").read_text()
        .splitlines()
        if line.strip() and not line.startswith("#")
    ]

    assert docs_requirements == docs_extra


def test_readthedocs_installs_existing_docs_requirements():
    """Read the Docs should install the committed docs requirements file."""
    rtd_config = (PROJECT_ROOT / ".readthedocs.yaml").read_text()

    assert "requirements: docs/requirements.txt" in rtd_config
    assert (PROJECT_ROOT / "docs" / "requirements.txt").is_file()
    assert not (PROJECT_ROOT / "requirements.txt").exists()


def test_documented_installation_commands_match_metadata():
    """Documentation should use current Python and console-script metadata."""
    metadata = tomllib.loads((PROJECT_ROOT / "pyproject.toml").read_text())
    docs_text = "\n".join(
        path.read_text()
        for path in (PROJECT_ROOT / "docs" / "source").glob("*.rst")
    )

    assert metadata["project"]["requires-python"] == ">=3.11"
    assert "Python >= 3.11" in docs_text
    for script_name in metadata["project"]["scripts"]:
        assert script_name in docs_text
    for stale_command in STALE_DOC_COMMANDS:
        assert stale_command not in docs_text
    assert "pip install ***" not in docs_text


def test_package_exports_resolve_to_defined_names():
    """Package-level public exports should refer to importable symbols."""
    import surface_pd.analysis as analysis
    import surface_pd.error as error
    import surface_pd.plot as plot
    import surface_pd.util as util

    for module in (analysis, error, plot, util):
        assert all(hasattr(module, name) for name in module.__all__)
