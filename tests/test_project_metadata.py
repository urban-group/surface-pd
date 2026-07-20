"""Tests for project packaging and dependency metadata."""

import tomllib
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SUPPORTED_PYTHON_VERSIONS = ("3.11", "3.12", "3.13", "3.14")
REQUIRES_PYTHON = ">=3.11,<3.15"
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

    assert metadata["project"]["requires-python"] == REQUIRES_PYTHON
    assert "Python >= 3.11, < 3.15" in docs_text
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


def test_ci_workflow_runs_verification_commands():
    """The CI workflow should mirror the documented local checks."""
    workflow = PROJECT_ROOT / ".github" / "workflows" / "ci.yml"
    workflow_text = workflow.read_text()

    matrix_versions = ", ".join(f'"{version}"' for version in (
        SUPPORTED_PYTHON_VERSIONS
    ))
    assert f"python-version: [{matrix_versions}]" in workflow_text
    assert 'python -m pip install -e ".[dev]"' in workflow_text
    assert "ruff check ." in workflow_text
    assert "python -m pytest" in workflow_text
    assert "python -m sphinx -b html docs/source docs/_build/html" in (
        workflow_text
    )
    assert "python -m pip wheel --no-deps --wheel-dir dist ." in (
        workflow_text
    )
    assert "actions/upload-artifact" in workflow_text


def test_python_metadata_matches_supported_versions():
    """Python version metadata should match the tested support policy."""
    metadata = tomllib.loads((PROJECT_ROOT / "pyproject.toml").read_text())
    classifiers = set(metadata["project"]["classifiers"])

    assert metadata["project"]["requires-python"] == REQUIRES_PYTHON
    for version in SUPPORTED_PYTHON_VERSIONS:
        assert f"Programming Language :: Python :: {version}" in classifiers


def test_release_metadata_describes_pre_one_beta():
    """Release metadata should match the documented pre-1.0 policy."""
    import surface_pd

    metadata = tomllib.loads((PROJECT_ROOT / "pyproject.toml").read_text())
    project = metadata["project"]

    assert "version" not in project
    assert "version" in project["dynamic"]
    assert metadata["tool"]["setuptools"]["dynamic"]["version"] == {
        "attr": "surface_pd._version.__version__"
    }
    assert surface_pd.__version__ == "0.1.0"
    assert "Development Status :: 4 - Beta" in project["classifiers"]
    assert "Development Status :: 5 - Production/Stable" not in (
        project["classifiers"]
    )


def test_release_policy_and_sphinx_use_authoritative_version():
    """Policy and docs configuration should not duplicate release literals."""
    policy = (PROJECT_ROOT / "docs" / "source" / "release_policy.rst")
    sphinx_config = (PROJECT_ROOT / "docs" / "source" / "conf.py").read_text()

    assert policy.is_file()
    assert "pre-1.0 releases may introduce breaking api changes" in (
        policy.read_text().lower()
    )
    assert 'release = "1.0.0"' not in sphinx_config
    assert "surface_pd/_version.py" in sphinx_config


def test_version_file_is_not_duplicated():
    """The package should not ship a second manually maintained version."""
    assert not (PROJECT_ROOT / "surface_pd" / "VERSION").exists()
