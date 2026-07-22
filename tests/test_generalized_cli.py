"""End-to-end tests for the generalized plotting command."""

import json
import subprocess
import sys
from pathlib import Path

import matplotlib
import pytest

from surface_pd.cli import surface_pd_plot
from surface_pd.configuration import PhaseDiagramConfiguration
from surface_pd.plot import CompositionColoring

matplotlib.use("Agg")

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _configuration_data():
    """Return a configuration exercising all four built-in models."""
    return {
        "schema_version": 1,
        "calculation_method": "test DFT method",
        "components": ["Li", "O", "X", "C"],
        "independent_chemical_potentials": {
            "Li": {
                "model": "intercalation_voltage",
                "reference_energy_ev_per_component": -2.0,
                "electrons_per_component": 1,
                "voltage_reference": "Li/Li+",
            },
            "O": {
                "model": "fixed_pressure_oxygen",
                "raw_o2_energy_ev_per_molecule": -10.0,
                "correction_ev_per_molecule": 1.0,
            },
            "X": {
                "model": "direct",
                "state_variable": "voltage",
                "reference_energy_ev_per_component": -1.0,
            },
            "C": {"model": "constant", "value_ev_per_component": 0.5},
        },
        "reference_phases": [],
        "diagram": {
            "x_axis": {
                "state_variable": "voltage",
                "coordinates": {
                    "kind": "values",
                    "values": [0.0, 1.0],
                },
                "label": "Potential",
                "unit": "V",
            },
            "y_axis": {
                "state_variable": "temperature",
                "coordinates": {
                    "kind": "values",
                    "values": [300.0, 600.0],
                },
                "label": "Temperature",
                "unit": "K",
            },
            "fixed_conditions": {},
        },
        "datasets": [
            {
                "dataset_id": "surface",
                "path": "phases.dat",
                "number_of_surfaces": 2,
                "column_overrides": {
                    "phase_id": "name",
                    "composition": {
                        "Li": "nLi",
                        "O": "nO",
                        "X": "nX",
                        "C": "nC",
                    },
                    "dft_energy_ev": "energy",
                    "surface_area_angstrom2": "area",
                },
            }
        ],
        "alignments": [],
        "rendering": {
            "coloring": {"mode": "phase_identity"},
            "colormap": "tab20",
            "invert_x_axis": False,
            "invert_y_axis": False,
        },
    }


def _write_inputs(tmp_path, data=None):
    """Write one configuration and compatible phase table."""
    data = data or _configuration_data()
    (tmp_path / "phases.dat").write_text(
        "name nLi nO nX nC energy area note\n"
        "alpha 1 2 1 0 -20 10 ignored\n"
        "beta 0 1 2 1 -15 10 ignored\n"
    )
    path = tmp_path / "config.json"
    path.write_text(json.dumps(data))
    return path


def test_help_exposes_only_generalized_input_and_explicit_actions():
    """Help should contain no chemistry-specific legacy abstraction."""
    result = subprocess.run(
        [sys.executable, "-m", "surface_pd.cli.surface_pd_plot", "--help"],
        cwd=PROJECT_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0
    assert "CONFIG" in result.stdout
    assert "--output" in result.stdout
    assert "--show" in result.stdout
    legacy_terms = ("lithium", "oxygen-like", "--low-T", "--save")
    assert not any(term in result.stdout for term in legacy_terms)


def test_cli_requires_an_explicit_output_or_show_action(tmp_path):
    """Reading a configuration alone should not imply a side effect."""
    path = _write_inputs(tmp_path)

    with pytest.raises(SystemExit, match="2"):
        surface_pd_plot.main([str(path)])


def test_prepare_phase_diagram_covers_all_models_and_identity_rendering(
    tmp_path,
):
    """The internal workflow should evaluate and render without I/O actions."""
    configuration = PhaseDiagramConfiguration.read_json(
        _write_inputs(tmp_path)
    )

    result, figure, axes, colorbar = surface_pd_plot._prepare_phase_diagram(
        configuration
    )

    assert result.phase_ids == ("surface:alpha", "surface:beta")
    assert result.surface_grand_potential_ev_per_angstrom2.shape == (2, 2, 2)
    assert axes.get_xlabel() == "Potential (V)"
    assert axes.get_ylabel() == "Temperature (K)"
    assert colorbar.ax.get_ylabel() == "Stable phase"
    assert not axes.xaxis_inverted()
    assert not axes.yaxis_inverted()
    figure.clear()


@pytest.mark.parametrize(
    ("coloring", "expected_normalization"),
    [
        (
            {
                "mode": "atomic_fraction",
                "component": "Li",
                "label": "Lithium fraction",
                "unit": "1",
            },
            "atomic_fraction",
        ),
        (
            {
                "mode": "component_ratio",
                "component": "Li",
                "reference_component": "O",
                "label": "Li/O ratio",
                "unit": "1",
            },
            "component_ratio",
        ),
    ],
)
def test_configuration_constructs_composition_coloring(
    coloring, expected_normalization
):
    """JSON rendering choices should become the reviewed coloring object."""
    data = _configuration_data()
    data["rendering"]["coloring"] = coloring

    configuration = PhaseDiagramConfiguration(data)

    actual = configuration.create_coloring()
    assert isinstance(actual, CompositionColoring)
    assert actual.normalization == expected_normalization


def test_prepare_phase_diagram_applies_rendering_options(tmp_path):
    """Colormap, composition coloring, and inversion should reach rendering."""
    data = _configuration_data()
    data["rendering"] = {
        "coloring": {
            "mode": "atomic_fraction",
            "component": "Li",
            "label": "Lithium fraction",
            "unit": "1",
        },
        "colormap": "plasma",
        "invert_x_axis": True,
        "invert_y_axis": True,
    }
    configuration = PhaseDiagramConfiguration.read_json(
        _write_inputs(tmp_path, data)
    )

    _, figure, axes, colorbar = surface_pd_plot._prepare_phase_diagram(
        configuration
    )

    assert axes.xaxis_inverted()
    assert axes.yaxis_inverted()
    assert colorbar.ax.get_ylabel() == "Lithium fraction (1)"
    assert colorbar.mappable.cmap.name == "plasma"
    figure.clear()


def test_cli_writes_output_without_showing(tmp_path, monkeypatch):
    """Output-only operation should save exactly the requested file."""
    path = _write_inputs(tmp_path)
    output = tmp_path / "diagram.png"
    shown = False

    def record_show():
        nonlocal shown
        shown = True

    monkeypatch.setattr(surface_pd_plot.plt, "show", record_show)

    surface_pd_plot.main([str(path), "--output", str(output)])

    assert output.is_file()
    assert output.stat().st_size > 0
    assert not shown


def test_cli_show_action_does_not_create_an_output(tmp_path, monkeypatch):
    """Display-only operation should call show without inventing a filename."""
    path = _write_inputs(tmp_path)
    calls = []
    monkeypatch.setattr(surface_pd_plot.plt, "show", lambda: calls.append(1))

    surface_pd_plot.main([str(path), "--show"])

    assert calls == [1]
    assert not (tmp_path / "diagram.pdf").exists()


def test_cli_can_save_and_show_in_one_invocation(tmp_path, monkeypatch):
    """The two explicit actions should compose without changing evaluation."""
    path = _write_inputs(tmp_path)
    output = tmp_path / "diagram.svg"
    calls = []
    monkeypatch.setattr(surface_pd_plot.plt, "show", lambda: calls.append(1))

    surface_pd_plot.main(
        [str(path), "--output", str(output), "--show"]
    )

    assert output.is_file()
    assert calls == [1]


def test_invalid_configuration_has_no_rendering_side_effects(
    tmp_path, monkeypatch
):
    """Validation failure should occur before saving or displaying."""
    path = tmp_path / "invalid.json"
    path.write_text("{}")
    output = tmp_path / "must-not-exist.png"
    calls = []
    monkeypatch.setattr(surface_pd_plot.plt, "show", lambda: calls.append(1))

    with pytest.raises(SystemExit, match="2"):
        surface_pd_plot.main(
            [str(path), "--output", str(output), "--show"]
        )

    assert not output.exists()
    assert calls == []


def test_cli_reports_table_context_without_a_traceback(tmp_path):
    """Expected input errors should remain concise command-line diagnostics."""
    path = _write_inputs(tmp_path)
    (tmp_path / "phases.dat").write_text(
        "name nLi nO nX nC energy area\n"
        "alpha 1 invalid 1 0 -20 10\n"
    )
    output = tmp_path / "must-not-exist.png"

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "surface_pd.cli.surface_pd_plot",
            str(path),
            "--output",
            str(output),
        ],
        cwd=PROJECT_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 2
    assert "dataset 'surface' row 2 column 'nO'" in result.stderr
    assert "Traceback" not in result.stderr
    assert not output.exists()


def test_cli_evaluates_configured_alignment(tmp_path):
    """The command workflow should accept direct aligned dataset views."""
    data = _configuration_data()
    data["components"] = ["A", "B", "C"]
    data["independent_chemical_potentials"] = {
        "A": {
            "model": "direct",
            "state_variable": "mu_A",
            "reference_energy_ev_per_component": 0.0,
        },
        "B": {
            "model": "direct",
            "state_variable": "mu_B",
            "reference_energy_ev_per_component": 0.0,
        },
    }
    data["reference_phases"] = [
        {
            "reference_id": "bulk-C",
            "composition": {"A": 0, "B": 0, "C": 1},
            "energy_ev_per_formula_unit": -5.0,
        }
    ]
    data["diagram"]["x_axis"]["state_variable"] = "mu_A"
    data["diagram"]["y_axis"]["state_variable"] = "mu_B"
    column_overrides = {
        "phase_id": "name",
        "composition": {"A": "nA", "B": "nB", "C": "nC"},
        "dft_energy_ev": "energy",
        "surface_area_angstrom2": "area",
    }
    data["datasets"] = [
        {
            "dataset_id": "root",
            "path": "root.dat",
            "number_of_surfaces": 2,
            "column_overrides": column_overrides,
        },
        {
            "dataset_id": "target",
            "path": "target.dat",
            "number_of_surfaces": 2,
            "column_overrides": column_overrides,
        },
    ]
    data["alignments"] = [
        {
            "reference_dataset_id": "root",
            "target_dataset_id": "target",
            "reference_anchor_phase_id": "anchor",
            "target_anchor_phase_id": "anchor",
            "bulk_reference_id": "bulk-C",
        }
    ]
    (tmp_path / "root.dat").write_text(
        "name nA nB nC energy area\nanchor 1 1 1 10 10\n"
    )
    (tmp_path / "target.dat").write_text(
        "name nA nB nC energy area\nanchor 1 1 2 30 10\n"
    )
    path = tmp_path / "aligned.json"
    path.write_text(json.dumps(data))

    result, figure, _, _ = surface_pd_plot._prepare_phase_diagram(
        PhaseDiagramConfiguration.read_json(path)
    )

    assert result.datasets[1].energy_offset_ev == pytest.approx(-25.0)
    figure.clear()
