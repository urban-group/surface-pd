"""Tests for project packaging and dependency metadata."""

import json
import re
import subprocess
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
SOFTWARE_CITATION_YEAR = "2024"
USER_DOCUMENTATION = (
    PROJECT_ROOT / "README.md",
    PROJECT_ROOT / "examples" / "README.md",
    *(PROJECT_ROOT / "docs" / "source").glob("*.rst"),
)


def _documented_autoclass_targets():
    """Return every fully qualified class target in the API reference."""
    targets = set()
    for filename in (
        "core.rst",
        "plot.rst",
        "thermodynamics.rst",
        "error.rst",
    ):
        text = (PROJECT_ROOT / "docs" / "source" / filename).read_text()
        targets.update(re.findall(r"^\.\. autoclass:: (\S+)$", text, re.M))
    return targets


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


def test_readthedocs_installs_project_with_docs_dependencies():
    """Read the Docs should install the package and its docs dependencies."""
    rtd_config = (PROJECT_ROOT / ".readthedocs.yaml").read_text()

    assert re.search(r'^    python: "3\.14"$', rtd_config, re.M)
    assert re.search(r"^    - method: pip$", rtd_config, re.M)
    assert re.search(r"^      path: \.$", rtd_config, re.M)
    assert re.search(r"^      extra_requirements:$", rtd_config, re.M)
    assert re.search(r"^        - docs$", rtd_config, re.M)
    assert (PROJECT_ROOT / "docs" / "requirements.txt").is_file()
    assert not (PROJECT_ROOT / "requirements.txt").exists()


def test_readthedocs_fails_on_sphinx_warnings():
    """Hosted documentation builds should reject Sphinx warnings."""
    rtd_config = (PROJECT_ROOT / ".readthedocs.yaml").read_text()

    assert re.search(r"^  fail_on_warning: true$", rtd_config, re.M)


def test_readme_exposes_latest_readthedocs_build_status():
    """The README should expose and link the hosted latest-build status."""
    readme = (PROJECT_ROOT / "README.md").read_text()

    assert (
        "https://readthedocs.org/projects/surface-pd/badge/?version=latest"
        in readme
    )
    assert (
        "https://surface-pd.readthedocs.io/en/latest/?badge=latest"
        in readme
    )


def test_readthedocs_maintainer_runbook_covers_recovery():
    """Maintainers should have a complete RTD verification runbook."""
    runbook_path = PROJECT_ROOT / ".github" / "READTHEDOCS.md"

    assert runbook_path.is_file()
    runbook = runbook_path.read_text()
    for required_text in (
        "urban-group/surface-pd",
        "Admin → Integrations",
        "Recent Activity",
        "Build version",
        "latest",
        "GitHub Actions",
        "Read the Docs",
        "API token",
    ):
        assert required_text in runbook


def test_repository_has_no_custom_readthedocs_deployment_trigger():
    """RTD publication should rely only on its standard Git integration."""
    workflow_text = "\n".join(
        path.read_text()
        for path in (PROJECT_ROOT / ".github" / "workflows").glob("*.y*ml")
    )

    assert "READTHEDOCS_TOKEN" not in workflow_text
    assert "api/v2/webhook" not in workflow_text
    assert "api/v3/projects" not in workflow_text


def test_repository_does_not_track_os_generated_metadata():
    """Git should ignore and never track macOS metadata files."""
    tracked_files = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.split("\0")
    tracked_ds_store = [
        path for path in tracked_files if Path(path).name == ".DS_Store"
    ]
    ignore_rules = (PROJECT_ROOT / ".gitignore").read_text().splitlines()

    assert ".DS_Store" in ignore_rules
    assert tracked_ds_store == []


def test_documented_installation_commands_match_metadata():
    """Documentation should use current Python and console-script metadata."""
    metadata = tomllib.loads((PROJECT_ROOT / "pyproject.toml").read_text())
    docs_text = "\n".join(path.read_text() for path in USER_DOCUMENTATION)

    assert metadata["project"]["requires-python"] == REQUIRES_PYTHON
    assert "Python >= 3.11, < 3.15" in docs_text
    for script_name in metadata["project"]["scripts"]:
        assert script_name in docs_text
    for stale_command in STALE_DOC_COMMANDS:
        assert stale_command not in docs_text
    assert "pip install ***" not in docs_text


def test_documented_example_paths_exist():
    """Committed enumeration configurations should resolve their inputs."""
    input_directory = (
        PROJECT_ROOT / "examples" / "enumeration-examples" / "input"
    )

    for input_path in input_directory.glob("*.json"):
        configuration = json.loads(input_path.read_text())
        target = PROJECT_ROOT / configuration["target_slab_path"]
        assert target.is_file(), f"{input_path}: {target}"

    documentation = "\n".join(path.read_text() for path in USER_DOCUMENTATION)
    assert "example/" not in documentation


def test_user_documentation_has_one_surface_pd_citation():
    """All citation locations should agree on canonical software metadata."""
    import surface_pd

    citation_paths = (
        PROJECT_ROOT / "README.md",
        PROJECT_ROOT / "examples" / "README.md",
        PROJECT_ROOT / "docs" / "source" / "how_to_cite.rst",
    )
    required_text = (
        "Li, Xinhao and Urban, Alexander",
        "surface-pd: Surface Phase Diagram Generator",
        f"version = {{{surface_pd.__version__}}}",
        f"year = {{{SOFTWARE_CITATION_YEAR}}}",
        "https://github.com/urban-group/surface-pd",
    )

    for path in citation_paths:
        citation = path.read_text()
        for text in required_text:
            assert text in citation, f"{path}: {text}"
        assert "InterPhon" not in citation


def test_readme_installation_matches_beta_repository_release():
    """Installation prose should not claim an unavailable stable release."""
    readme = (PROJECT_ROOT / "README.md").read_text()
    installation = (
        PROJECT_ROOT / "docs" / "source" / "installation.rst"
    ).read_text()

    assert "stable version" not in readme
    assert "$ pip install surface-pd" not in readme
    assert "pip install git+https://github.com/urban-group/surface-pd.git" in (
        readme
    )
    assert "git clone https://github.com/urban-group/surface-pd.git" in (
        installation
    )
    assert "git@github.com:urban-group/surface-pd.git" not in installation


def test_package_exports_resolve_to_defined_names():
    """Package-level public exports should refer to importable symbols."""
    import surface_pd.analysis as analysis
    import surface_pd.error as error
    import surface_pd.plot as plot
    import surface_pd.thermodynamics as thermodynamics
    import surface_pd.util as util

    for module in (analysis, error, plot, thermodynamics, util):
        assert all(hasattr(module, name) for name in module.__all__)


def test_ci_workflow_runs_verification_commands():
    """The CI workflow should mirror the documented local checks."""
    workflow = PROJECT_ROOT / ".github" / "workflows" / "ci.yml"
    workflow_text = workflow.read_text()
    normalized_workflow = " ".join(workflow_text.split())

    matrix_versions = ", ".join(f'"{version}"' for version in (
        SUPPORTED_PYTHON_VERSIONS
    ))
    assert f"python-version: [{matrix_versions}]" in workflow_text
    assert 'python -m pip install -e ".[dev]"' in workflow_text
    assert "ruff check ." in workflow_text
    assert "python -m pytest" in workflow_text
    assert (
        "python -m sphinx -W --keep-going -b html "
        "docs/source docs/_build/html"
    ) in normalized_workflow
    assert (
        "python -m sphinx -W --keep-going -b doctest "
        "docs/source docs/_build/doctest"
    ) in normalized_workflow
    assert "python -m pip wheel --no-deps --wheel-dir dist ." in (
        workflow_text
    )
    assert "actions/upload-artifact" in workflow_text


def test_sphinx_executes_two_principal_workflow_examples():
    """Tutorials and CI should retain two named executable examples."""
    tutorials_path = PROJECT_ROOT / "docs" / "source" / "tutorials.rst"
    tutorials = tutorials_path.read_text()
    examples_path = PROJECT_ROOT / "docs" / "source" / "python_examples.rst"

    assert "python_examples" in tutorials
    assert examples_path.is_file()
    examples = examples_path.read_text()
    testcode_groups = re.findall(
        r"^\.\. testcode:: ([\w-]+)$", examples, re.M
    )
    testoutput_groups = re.findall(
        r"^\.\. testoutput:: ([\w-]+)$", examples, re.M
    )
    expected_groups = ["enumeration", "phase-diagram"]

    assert testcode_groups == expected_groups
    assert testoutput_groups == expected_groups


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
    assert "from surface_pd import __version__" in sphinx_config


def test_sphinx_documents_exact_canonical_public_api():
    """Autodoc should cover supported classes through installed paths."""
    from surface_pd import core, error, plot, thermodynamics

    expected = {
        *(f"surface_pd.core.{name}" for name in core.__all__),
        *(f"surface_pd.plot.{name}" for name in plot.__all__),
        *(
            f"surface_pd.thermodynamics.{name}"
            for name in thermodynamics.__all__
        ),
        *(f"surface_pd.error.{name}" for name in error.__all__),
    }

    assert _documented_autoclass_targets() == expected


def test_sphinx_imports_installed_package_without_dependency_mocks():
    """Autodoc should exercise real package imports and dependencies."""
    sphinx_config = (PROJECT_ROOT / "docs" / "source" / "conf.py").read_text()

    assert "autodoc_mock_imports" not in sphinx_config
    assert "sys.path" not in sphinx_config
    assert 'os.path.abspath("../../surface_pd")' not in sphinx_config


def test_version_file_is_not_duplicated():
    """The package should not ship a second manually maintained version."""
    assert not (PROJECT_ROOT / "surface_pd" / "VERSION").exists()
