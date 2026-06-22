"""Tests for project packaging and dependency metadata."""

import tomllib
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _dependency_names(dependencies):
    """Return normalized package names from dependency specifiers."""
    return {dependency.split(">=", maxsplit=1)[0] for dependency in dependencies}


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
