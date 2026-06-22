"""Tests for command-line entry point wiring."""

import importlib
import subprocess
import sys
import tomllib
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
EXPECTED_ENTRY_POINTS = {
    "surface-enumeration": "surface_pd.cli.surface_enumeration:main",
    "surface-pd-plot": "surface_pd.cli.surface_pd_plot:main",
    "discharge-pd-gene": "surface_pd.cli.discharge_pd_gene:main",
}


def test_project_scripts_target_cli_modules():
    """Console scripts should target importable package modules."""
    metadata = tomllib.loads((PROJECT_ROOT / "pyproject.toml").read_text())

    assert metadata["project"]["scripts"] == EXPECTED_ENTRY_POINTS


@pytest.mark.parametrize("target", EXPECTED_ENTRY_POINTS.values())
def test_cli_main_callables_are_importable(target):
    """Each configured entry point should resolve to a callable main."""
    module_name, function_name = target.split(":")

    module = importlib.import_module(module_name)

    assert callable(getattr(module, function_name))


@pytest.mark.parametrize(
    "module_name",
    [
        "surface_pd.cli.surface_enumeration",
        "surface_pd.cli.surface_pd_plot",
        "surface_pd.cli.discharge_pd_gene",
    ],
)
def test_cli_modules_show_help(module_name):
    """CLI modules should be runnable and expose argparse help."""
    result = subprocess.run(
        [sys.executable, "-m", module_name, "--help"],
        cwd=PROJECT_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0
    assert "usage:" in result.stdout
