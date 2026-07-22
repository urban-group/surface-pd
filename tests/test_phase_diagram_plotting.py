"""Tests for chemistry-independent phase-diagram rendering."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.collections import LineCollection, QuadMesh

from surface_pd.configuration import PhaseDiagramConfiguration
from surface_pd.plot import CompositionColoring, plot_phase_diagram
from surface_pd.thermodynamics import (
    DiagramAxis,
    DirectChemicalPotential,
    GrandPotentialModel,
    Phase,
    PhaseDataset,
    PhaseDiagramSpecification,
)

_PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _diagram_result(energies=(0.0, 1.0)):
    """Return a small result whose representative changes along x."""
    model = GrandPotentialModel(
        ("A", "B"),
        {
            "A": DirectChemicalPotential("x", 0.0),
            "B": DirectChemicalPotential("y", 0.0),
        },
        (),
    )
    dataset = PhaseDataset(
        "surfaces",
        ("A", "B"),
        (
            Phase("equal", {"A": 1, "B": 1}, energies[0], 1.0, 1),
            Phase("A-rich", {"A": 2, "B": 1}, energies[1], 1.0, 1),
        ),
        "plot test method",
    )
    specification = PhaseDiagramSpecification(
        DiagramAxis("x", (0.0, 2.0), "Horizontal condition", "X-unit"),
        DiagramAxis("y", (0.0, 1.0), "Vertical condition", "Y-unit"),
        {},
    )
    return specification.evaluate(model, (dataset,))


def _quad_mesh(axes):
    """Return the single colored mesh created by the renderer."""
    meshes = [
        collection
        for collection in axes.collections
        if isinstance(collection, QuadMesh)
    ]
    assert len(meshes) == 1
    return meshes[0]


def _boundary_lines(axes):
    """Return the phase-boundary collection created by the renderer."""
    return [
        collection
        for collection in axes.collections
        if isinstance(collection, LineCollection)
    ]


def test_identity_renderer_uses_qualified_phase_labels_without_side_effects(
    monkeypatch,
):
    """Rendering should return objects without external side effects."""
    result = _diagram_result()
    stable_before = result.stable_phase_mask.copy()
    energies_before = result.surface_grand_potential_ev_per_angstrom2.copy()
    representatives_before = result.representative_phase_indices.copy()
    monkeypatch.setattr(
        plt,
        "show",
        lambda: pytest.fail("renderer must not call pyplot.show"),
    )

    figure, axes, colorbar = plot_phase_diagram(result)

    assert figure is axes.figure
    assert colorbar.ax.figure is figure
    assert axes.get_xlabel() == "Horizontal condition (X-unit)"
    assert axes.get_ylabel() == "Vertical condition (Y-unit)"
    assert [label.get_text() for label in colorbar.ax.get_yticklabels()] == [
        "surfaces:equal",
        "surfaces:A-rich",
    ]
    np.testing.assert_array_equal(
        _quad_mesh(axes).get_array(), [[0, 1], [0, 1]]
    )
    np.testing.assert_array_equal(result.stable_phase_mask, stable_before)
    np.testing.assert_array_equal(
        result.surface_grand_potential_ev_per_angstrom2, energies_before
    )
    np.testing.assert_array_equal(
        result.representative_phase_indices, representatives_before
    )
    assert len(_boundary_lines(axes)) == 1
    plt.close(figure)


def test_boundary_overlay_matches_nearest_mesh_cell_edges():
    """The visual boundary should coincide with the colored-cell boundary."""
    result = _diagram_result()

    figure, axes, _ = plot_phase_diagram(
        result,
        boundary_color="magenta",
        boundary_linewidth=2.5,
    )

    boundaries = _boundary_lines(axes)
    assert len(boundaries) == 1
    np.testing.assert_allclose(
        boundaries[0].get_segments(),
        [
            [[1.0, -0.5], [1.0, 0.5]],
            [[1.0, 0.5], [1.0, 1.5]],
        ],
    )
    np.testing.assert_allclose(boundaries[0].get_linewidths(), [2.5])
    np.testing.assert_allclose(
        boundaries[0].get_colors(), [[1.0, 0.0, 1.0, 1.0]]
    )
    plt.close(figure)


def test_uniform_and_disabled_boundaries_add_no_line_collection():
    """Uniform fields and an explicit opt-out should have no internal lines."""
    uniform = _diagram_result(energies=(0.0, 10.0))
    uniform_figure, uniform_axes, _ = plot_phase_diagram(uniform)
    disabled_figure, disabled_axes, _ = plot_phase_diagram(
        _diagram_result(), boundary_color=None
    )

    assert _boundary_lines(uniform_axes) == []
    assert _boundary_lines(disabled_axes) == []
    plt.close(uniform_figure)
    plt.close(disabled_figure)


def test_lno_001_boundary_layout_regression():
    """The maintained voltage-temperature example should retain its layout."""
    configuration = PhaseDiagramConfiguration.read_json(
        _PROJECT_ROOT
        / "examples/plotting-examples/lno-001-scan/charge.json"
    )
    result = configuration.diagram_specification.evaluate(
        configuration.model, configuration.load_datasets()
    )

    figure, axes, _ = plot_phase_diagram(result)

    assert np.unique(result.representative_phase_indices).size == 13
    boundaries = _boundary_lines(axes)
    assert len(boundaries) == 1
    assert len(boundaries[0].get_segments()) == 1215
    plt.close(figure)


def test_boundary_color_must_be_understood_by_matplotlib():
    """Invalid boundary colors should fail before creating a figure."""
    with pytest.raises(ValueError, match="boundary_color"):
        plot_phase_diagram(_diagram_result(), boundary_color="not-a-color")


@pytest.mark.parametrize("boundary_linewidth", [0.0, -1.0, np.inf])
def test_boundary_linewidth_must_be_finite_and_positive(boundary_linewidth):
    """Enabled boundaries require a visible finite line width."""
    with pytest.raises(ValueError, match="boundary_linewidth"):
        plot_phase_diagram(
            _diagram_result(), boundary_linewidth=boundary_linewidth
        )


def test_renderer_uses_first_phase_for_ties_without_losing_tie_mask():
    """A raster representative should not overwrite numerical co-stability."""
    result = _diagram_result(energies=(0.0, 0.0))
    assert result.stable_phase_ids_at(0, 0) == (
        "surfaces:equal",
        "surfaces:A-rich",
    )

    figure, axes, _ = plot_phase_diagram(result)

    assert _quad_mesh(axes).get_array()[0, 0] == 0
    assert result.stable_phase_ids_at(0, 0) == (
        "surfaces:equal",
        "surfaces:A-rich",
    )
    plt.close(figure)


def test_atomic_fraction_coloring_is_explicit_and_component_neutral():
    """Atomic fractions should be calculated from absolute composition."""
    result = _diagram_result()
    coloring = CompositionColoring(
        component="A",
        normalization="atomic_fraction",
        reference_component=None,
        label="A atomic fraction",
        unit="fraction",
    )

    np.testing.assert_allclose(coloring.phase_values(result), [0.5, 2 / 3])
    figure, axes, colorbar = plot_phase_diagram(
        result, coloring=coloring, cmap="viridis"
    )

    np.testing.assert_allclose(
        _quad_mesh(axes).get_array(),
        [[0.5, 2 / 3], [0.5, 2 / 3]],
    )
    assert colorbar.ax.get_ylabel() == "A atomic fraction (fraction)"
    plt.close(figure)


def test_component_ratio_coloring_uses_declared_denominator():
    """Component ratios should never infer a host or transition metal."""
    result = _diagram_result()
    coloring = CompositionColoring(
        component="A",
        normalization="component_ratio",
        reference_component="B",
        label="A per B",
        unit="A/B",
    )

    np.testing.assert_allclose(coloring.phase_values(result), [1.0, 2.0])


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"normalization": "unknown"}, "normalization"),
        ({"normalization": "component_ratio"}, "reference_component"),
        (
            {
                "normalization": "atomic_fraction",
                "reference_component": "B",
            },
            "reference_component",
        ),
    ],
)
def test_composition_coloring_rejects_ambiguous_configuration(
    kwargs, message
):
    """Normalization meaning should be complete at construction."""
    values = {
        "component": "A",
        "normalization": "atomic_fraction",
        "reference_component": None,
        "label": "A content",
        "unit": "fraction",
    }
    values.update(kwargs)

    with pytest.raises(ValueError, match=message):
        CompositionColoring(**values)


def test_composition_coloring_rejects_missing_or_zero_components():
    """Every plotted phase must define a finite requested composition ratio."""
    result = _diagram_result()
    missing = CompositionColoring(
        "C", "atomic_fraction", None, "C fraction", "fraction"
    )
    zero_denominator = CompositionColoring(
        "A", "component_ratio", "C", "A per C", "A/C"
    )

    with pytest.raises(ValueError, match="component.*C"):
        missing.phase_values(result)
    with pytest.raises(ValueError, match="reference component.*C"):
        zero_denominator.phase_values(result)


def test_axis_inversion_and_existing_axes_are_presentation_only():
    """Presentation controls should not alter numerical results."""
    result = _diagram_result()
    x_before = result.x_mesh.copy()
    figure, supplied_axes = plt.subplots()

    returned_figure, returned_axes, _ = plot_phase_diagram(
        result,
        ax=supplied_axes,
        invert_x_axis=True,
        invert_y_axis=True,
    )

    assert returned_figure is figure
    assert returned_axes is supplied_axes
    assert returned_axes.xaxis_inverted()
    assert returned_axes.yaxis_inverted()
    np.testing.assert_array_equal(result.x_mesh, x_before)
    plt.close(figure)
