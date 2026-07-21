"""Tests for command-line entry point wiring."""

import importlib
import subprocess
import sys
import tomllib
import warnings
from pathlib import Path

import pytest

from surface_pd.cli.discharge_pd_gene import create_discharge_pd
from surface_pd.plot._phase_data_io import _read_phase_diagram_file

PROJECT_ROOT = Path(__file__).resolve().parents[1]
EXPECTED_ENTRY_POINTS = {
    "surface-enumeration": "surface_pd.cli.surface_enumeration:main",
    "surface-pd-plot": "surface_pd.cli.surface_pd_plot:main",
    "generate-discharge-pd": "surface_pd.cli.discharge_pd_gene:main",
}


def test_surface_enumeration_import_does_not_suppress_warnings(monkeypatch):
    """Importing the enumeration CLI should not suppress all warnings."""
    module_name = "surface_pd.cli.surface_enumeration"
    previous_module = sys.modules.pop(module_name, None)
    global_ignore_filter = ("ignore", None, Warning, None, 0)
    baseline_filters = [("default", None, Warning, None, 0)]
    monkeypatch.setattr(warnings, "filters", list(baseline_filters))

    try:
        importlib.import_module(module_name)

        assert global_ignore_filter not in warnings.filters
    finally:
        sys.modules.pop(module_name, None)
        if previous_module is not None:
            sys.modules[module_name] = previous_module


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


def test_discharge_generator_preserves_reference_metadata(
    tmp_path,
    monkeypatch,
):
    """Generated phase data should retain its scientific provenance."""
    charge_path = tmp_path / "charge.dat"
    charge_path.write_text(
        "# method = custom; U_Ni=? eV\n"
        "# reference_li_ev_per_atom = -1.0\n"
        "# reference_o2_raw_ev_per_molecule = -2.0\n"
        "# reference_o2_correction_ev_per_molecule = 0.25\n"
        "# reference_bulk_litmo2_ev_per_formula_unit = -3.0\n"
        "structure Li Ni O E a b gamma\n"
        "p1 1 1 2 -4 2 3 90\n"
        "p2 0 1 1 -2 2 3 90\n"
    )
    monkeypatch.chdir(tmp_path)

    create_discharge_pd(
        [str(charge_path)],
        charge_pd_end_composition=1.0,
        lithium_like_species="Li",
        oxygen_like_species="O",
        save=True,
    )

    _, input_references = _read_phase_diagram_file(charge_path)
    output_dataframe, output_references = _read_phase_diagram_file(
        tmp_path / "discharge-data1.dat"
    )
    assert output_references == input_references
    assert list(output_dataframe["dft_energy"]) == [-4, -2]
