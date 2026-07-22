"""End-to-end tests for the generalized plotting command."""

import json
import subprocess
import sys
from pathlib import Path

import matplotlib
import pytest
from matplotlib.colorbar import Colorbar
from matplotlib.legend import Legend

from surface_pd.cli import surface_pd_plot
from surface_pd.configuration import PhaseDiagramConfiguration

matplotlib.use("Agg")

PROJECT_ROOT = Path(__file__).resolve().parents[1]
_EVALUATION = {
    "x_range": (0.0, 1.0),
    "y_range": (300.0, 600.0),
    "mesh_points": 2,
}


def _cli_arguments(path):
    """Return one explicit quick-inspection evaluation domain."""
    return [
        str(path),
        "--x-range",
        "0",
        "1",
        "--y-range",
        "300",
        "600",
    ]


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
                "label": "Potential",
                "unit": "V",
            },
            "y_axis": {
                "state_variable": "temperature",
                "label": "Temperature",
                "unit": "K",
            },
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
    assert "--x-range" in result.stdout
    assert "--y-range" in result.stdout
    assert "--mesh-points" in result.stdout
    assert "--condition" in result.stdout
    assert "--color-component" in result.stdout
    assert "--phase-identity-colors" in result.stdout
    assert "--colormap" not in result.stdout
    assert "--boundary-linewidth" not in result.stdout
    legacy_terms = ("lithium", "oxygen-like", "--low-T", "--save")
    assert not any(term in result.stdout for term in legacy_terms)


def test_cli_requires_an_explicit_output_or_show_action(tmp_path):
    """Reading a configuration alone should not imply a side effect."""
    path = _write_inputs(tmp_path)

    with pytest.raises(SystemExit, match="2"):
        surface_pd_plot.main(_cli_arguments(path))


def test_prepare_phase_diagram_defaults_to_first_independent_component(
    tmp_path,
):
    """The internal workflow should evaluate and render without I/O actions."""
    configuration = PhaseDiagramConfiguration.read_json(
        _write_inputs(tmp_path)
    )

    result, figure, axes, colorbar = surface_pd_plot._prepare_phase_diagram(
        configuration, **_EVALUATION
    )

    assert result.phase_ids == ("surface:alpha", "surface:beta")
    assert result.surface_grand_potential_ev_per_angstrom2.shape == (2, 2, 2)
    assert axes.get_xlabel() == "Potential (V)"
    assert axes.get_ylabel() == "Temperature (K)"
    assert colorbar.ax.get_ylabel() == "Li atomic fraction (1)"
    assert not axes.xaxis_inverted()
    assert not axes.yaxis_inverted()
    figure.clear()


def test_prepare_phase_diagram_uses_explicit_range_and_mesh_density(tmp_path):
    """CLI evaluation controls should determine the numerical mesh."""
    configuration = PhaseDiagramConfiguration.read_json(
        _write_inputs(tmp_path)
    )

    result, figure, _, _ = surface_pd_plot._prepare_phase_diagram(
        configuration,
        x_range=(-2.0, 3.0),
        y_range=(250.0, 1250.0),
        mesh_points=5,
    )

    assert result.mesh_shape == (5, 5)
    assert result.specification.x_axis.values.tolist() == [
        -2.0,
        -0.75,
        0.5,
        1.75,
        3.0,
    ]
    assert result.specification.y_axis.values[[0, -1]].tolist() == [
        250.0,
        1250.0,
    ]
    figure.clear()


def test_cli_fixed_condition_parser_is_explicit():
    """Additional thermodynamic slices should use NAME=VALUE syntax."""
    assert surface_pd_plot._fixed_condition("delta_mu=-0.25") == (
        "delta_mu",
        -0.25,
    )
    with pytest.raises(SystemExit):
        surface_pd_plot._parser().parse_args(
            [
                "config.json",
                "--x-range",
                "0",
                "1",
                "--y-range",
                "0",
                "1",
                "--condition",
                "invalid",
            ]
        )


def test_prepare_phase_diagram_accepts_simple_cli_coloring_choices(tmp_path):
    """Quick inspection may select a component or discrete phase identity."""
    configuration = PhaseDiagramConfiguration.read_json(
        _write_inputs(tmp_path)
    )

    _, figure, _, colorbar = surface_pd_plot._prepare_phase_diagram(
        configuration, color_component="O", **_EVALUATION
    )
    _, identity_figure, _, identity_legend = (
        surface_pd_plot._prepare_phase_diagram(
            configuration, phase_identity_colors=True, **_EVALUATION
        )
    )

    assert isinstance(colorbar, Colorbar)
    assert colorbar.ax.get_ylabel() == "O atomic fraction (1)"
    assert isinstance(identity_legend, Legend)
    assert identity_legend.get_title().get_text() == "Stable phase"
    figure.clear()
    identity_figure.clear()


def test_prepare_rejects_nonindependent_color_component(tmp_path):
    """CLI composition selection should reflect model provenance."""
    configuration = PhaseDiagramConfiguration.read_json(
        _write_inputs(tmp_path)
    )

    with pytest.raises(ValueError, match="independent component"):
        surface_pd_plot._prepare_phase_diagram(
            configuration, color_component="not-present", **_EVALUATION
        )


def test_cli_writes_output_without_showing(tmp_path, monkeypatch):
    """Output-only operation should save exactly the requested file."""
    path = _write_inputs(tmp_path)
    output = tmp_path / "diagram.png"
    shown = False

    def record_show():
        nonlocal shown
        shown = True

    monkeypatch.setattr(surface_pd_plot.plt, "show", record_show)

    surface_pd_plot.main(
        [*_cli_arguments(path), "--output", str(output)]
    )

    assert output.is_file()
    assert output.stat().st_size > 0
    assert not shown


def test_cli_show_action_does_not_create_an_output(tmp_path, monkeypatch):
    """Display-only operation should call show without inventing a filename."""
    path = _write_inputs(tmp_path)
    calls = []
    monkeypatch.setattr(surface_pd_plot.plt, "show", lambda: calls.append(1))

    surface_pd_plot.main([*_cli_arguments(path), "--show"])

    assert calls == [1]
    assert not (tmp_path / "diagram.pdf").exists()


def test_cli_can_save_and_show_in_one_invocation(tmp_path, monkeypatch):
    """The two explicit actions should compose without changing evaluation."""
    path = _write_inputs(tmp_path)
    output = tmp_path / "diagram.svg"
    calls = []
    monkeypatch.setattr(surface_pd_plot.plt, "show", lambda: calls.append(1))

    surface_pd_plot.main(
        [*_cli_arguments(path), "--output", str(output), "--show"]
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
            [*_cli_arguments(path), "--output", str(output), "--show"]
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
            *_cli_arguments(path),
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
        PhaseDiagramConfiguration.read_json(path), **_EVALUATION
    )

    assert result.datasets[1].energy_offset_ev == pytest.approx(-25.0)
    figure.clear()
