"""Tests for command-line entry point wiring."""

import importlib
import subprocess
import sys
import tomllib
import warnings
from pathlib import Path

import numpy as np
import pytest

from surface_pd.cli.discharge_pd_gene import create_discharge_pd
from surface_pd.cli.surface_pd_plot import _prepare_surface_pd_data
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


def test_surface_plot_cli_uses_file_reference_metadata():
    """CLI help should not offer a fragile method-name selector."""
    result = subprocess.run(
        [sys.executable, "-m", "surface_pd.cli.surface_pd_plot", "--help"],
        cwd=PROJECT_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert "--functional" not in result.stdout
    assert "reference" in result.stdout.lower()


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


def _write_phase_fixture(path, *, li, tm, oxygen, energy, bulk=-3.0):
    """Write one self-contained phase row for CLI preparation tests."""
    path.write_text(
        "# method = test; U_Ni=? eV\n"
        "# reference_li_ev_per_atom = -1.0\n"
        "# reference_o2_raw_ev_per_molecule = -2.0\n"
        "# reference_o2_correction_ev_per_molecule = 0.0\n"
        f"# reference_bulk_litmo2_ev_per_formula_unit = {bulk}\n"
        "structure Li Ni O E a b gamma\n"
        f"phase {li} {tm} {oxygen} {energy} 2 3 90\n"
    )


def test_prepare_surface_pd_data_handles_single_file(tmp_path):
    """Single-file preparation should preserve phase and mesh axes."""
    path = tmp_path / "single.dat"
    _write_phase_fixture(path, li=0, tm=1, oxygen=2, energy=-3)
    mesh = np.zeros((2, 3))

    dataframe, energies, groups = _prepare_surface_pd_data(
        [path], "Li", "O", mesh, np.full_like(mesh, 300.0)
    )

    assert len(dataframe) == 1
    assert energies.shape == (1, 2, 3)
    assert groups == 1


def test_prepare_surface_pd_data_aligns_two_files(tmp_path):
    """Two-file preparation should align and concatenate on the phase axis."""
    first_path = tmp_path / "first.dat"
    second_path = tmp_path / "second.dat"
    _write_phase_fixture(first_path, li=0, tm=1, oxygen=2, energy=-3)
    _write_phase_fixture(second_path, li=1, tm=2, oxygen=4, energy=-7)
    mesh = np.zeros((2, 2))

    dataframe, energies, groups = _prepare_surface_pd_data(
        [first_path, second_path],
        "Li",
        "O",
        mesh,
        np.full_like(mesh, 300.0),
    )

    assert len(dataframe) == 2
    assert energies.shape == (2, 2, 2)
    assert groups == 2


def test_prepare_surface_pd_data_aligns_committed_lno001_pair():
    """The maintained complementary facet example should execute end to end."""
    example_root = PROJECT_ROOT / "examples" / "plotting-examples" / "LNO-001"
    data_files = [
        example_root / "LNO-001-Li-SCAN-rVV10-U-spd-charge.dat",
        example_root / "LNO-001-Ni-SCAN-rVV10-U-spd-charge.dat",
    ]

    dataframe, energies, groups = _prepare_surface_pd_data(
        data_files,
        "Li",
        "O",
        np.zeros((1, 1)),
        np.full((1, 1), 300.0),
    )

    assert len(dataframe) == 50
    assert energies.shape == (50, 1, 1)
    assert np.isfinite(energies).all()
    assert groups == 2


@pytest.mark.parametrize("number_of_files", [0, 3])
def test_prepare_surface_pd_data_rejects_invalid_file_count(number_of_files):
    """Only the scientifically defined one- and two-file paths are valid."""
    with pytest.raises(ValueError, match="Exactly one or two"):
        _prepare_surface_pd_data(
            ["unused.dat"] * number_of_files,
            "Li",
            "O",
            np.zeros(1),
            np.full(1, 300.0),
        )


def test_surface_plot_cli_rejects_three_data_files():
    """Unsupported file counts should produce an argparse usage error."""
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "surface_pd.cli.surface_pd_plot",
            "one.dat",
            "two.dat",
            "three.dat",
        ],
        cwd=PROJECT_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 2
    assert "exactly one or two DATA_FILE" in result.stderr
