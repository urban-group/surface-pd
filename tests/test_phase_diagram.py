"""Tests for generalized two-dimensional numerical phase diagrams."""

from pathlib import Path

import numpy as np
import pytest

from surface_pd.thermodynamics import (
    DatasetAlignment,
    DiagramAxis,
    DirectChemicalPotential,
    GrandPotentialModel,
    Phase,
    PhaseDataset,
    PhaseDiagramResult,
    PhaseDiagramSpecification,
    ReferencePhase,
)

_METHOD = "diagram test method"


def test_numerical_phase_diagram_module_does_not_import_matplotlib():
    """Numerical evaluation should remain independent of rendering."""
    module_path = (
        Path(__file__).resolve().parents[1]
        / "surface_pd"
        / "thermodynamics"
        / "phase_diagram.py"
    )

    assert "matplotlib" not in module_path.read_text().lower()


def _axis(variable, values=(0.0, 1.0), label=None, unit="arb."):
    """Return a concise numerical diagram axis."""
    display_label = variable.title() if label is None else label
    return DiagramAxis(variable, values, display_label, unit)


def _model(variable_names):
    """Return a model with one direct potential per named state variable."""
    components = tuple(f"C{index}" for index in range(len(variable_names)))
    models = {
        component: DirectChemicalPotential(variable, 0.0)
        for component, variable in zip(
            components, variable_names, strict=True
        )
    }
    return GrandPotentialModel(components, models, ())


def _dataset(components, dataset_id="dataset", energies=(0.0, 1.0)):
    """Return phases with identical composition and selected energies."""
    composition = dict.fromkeys(components, 1)
    phases = tuple(
        Phase(f"p{index}", composition, energy, 2.0, 1)
        for index, energy in enumerate(energies)
    )
    return PhaseDataset(dataset_id, components, phases, _METHOD)


@pytest.mark.parametrize(
    "values",
    [(), (1.0,), (0.0, 0.0), (1.0, 0.0), (0.0, np.nan), [[0.0, 1.0]]],
)
def test_diagram_axis_requires_increasing_finite_coordinates(values):
    """Numerical axes should be directly suitable for conventional plots."""
    with pytest.raises((TypeError, ValueError), match="values"):
        _axis("voltage", values)


@pytest.mark.parametrize(
    ("variable", "label", "unit"),
    [("bad name", "Voltage", "V"), ("voltage", "", "V"),
     ("voltage", "Voltage", "")],
)
def test_diagram_axis_requires_explicit_identity_and_labels(
    variable, label, unit
):
    """State identity, display label, and units should not be inferred."""
    with pytest.raises((TypeError, ValueError)):
        _axis(variable, label=label, unit=unit)


def test_diagram_axis_owns_read_only_values():
    """Caller mutation should not change a configured axis."""
    values = np.array([0.0, 1.0])
    axis = _axis("voltage", values)
    values[0] = 99.0

    np.testing.assert_array_equal(axis.values, [0.0, 1.0])
    assert not axis.values.flags.writeable


@pytest.mark.parametrize(
    ("x_name", "y_name", "fixed_name"),
    [
        ("voltage", "temperature", "delta_mu"),
        ("voltage", "delta_mu", "temperature"),
        ("temperature", "delta_mu", "voltage"),
    ],
)
def test_required_axis_pairs_preserve_xy_mesh_orientation(
    x_name, y_name, fixed_name
):
    """All initially supported axis pairs should use shape (n_y, n_x)."""
    model = _model((x_name, y_name, fixed_name))
    dataset = _dataset(model.components)
    specification = PhaseDiagramSpecification(
        _axis(x_name, (1.0, 2.0, 3.0)),
        _axis(y_name, (10.0, 20.0)),
        {fixed_name: 7.0},
    )

    result = specification.evaluate(model, (dataset,))

    assert isinstance(result, PhaseDiagramResult)
    assert result.mesh_shape == (2, 3)
    np.testing.assert_array_equal(
        result.x_mesh,
        [[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]],
    )
    np.testing.assert_array_equal(
        result.y_mesh,
        [[10.0, 10.0, 10.0], [20.0, 20.0, 20.0]],
    )
    potentials = result.dataset_results[0].chemical_potentials_ev
    np.testing.assert_array_equal(potentials["C0"], result.x_mesh)
    np.testing.assert_array_equal(potentials["C1"], result.y_mesh)
    np.testing.assert_array_equal(potentials["C2"], np.full((2, 3), 7.0))


def test_specification_requires_distinct_axes_and_disjoint_fixed_conditions():
    """Each state variable should have exactly one source in a diagram."""
    voltage = _axis("voltage")
    with pytest.raises(ValueError, match="distinct"):
        PhaseDiagramSpecification(voltage, voltage, {})
    with pytest.raises(ValueError, match="fixed_conditions.*axis"):
        PhaseDiagramSpecification(
            voltage, _axis("temperature"), {"voltage": 1.0}
        )


@pytest.mark.parametrize(
    ("fixed", "message"),
    [({}, "missing.*delta_mu"),
     ({"delta_mu": 0.0, "unused": 1.0}, "unused.*unused")],
)
def test_specification_requires_exact_model_state_variables(fixed, message):
    """Missing conditions and likely configuration typos should fail early."""
    model = _model(("voltage", "temperature", "delta_mu"))
    specification = PhaseDiagramSpecification(
        _axis("voltage"), _axis("temperature"), fixed
    )

    with pytest.raises(ValueError, match=message):
        specification.evaluate(model, (_dataset(model.components),))


@pytest.mark.parametrize("value", [True, [1.0], np.nan, np.inf])
def test_fixed_conditions_must_be_finite_scalars(value):
    """A fixed condition should not introduce a hidden third mesh axis."""
    with pytest.raises((TypeError, ValueError), match="fixed_conditions"):
        PhaseDiagramSpecification(
            _axis("voltage"), _axis("temperature"), {"delta_mu": value}
        )


def test_result_preserves_ties_and_complete_energy_tensor():
    """Numerical co-stability should not be collapsed for rendering."""
    model = _model(("voltage", "temperature"))
    dataset = _dataset(
        model.components,
        energies=(0.0, 1.0e-10, 3.0e-10),
    )
    result = PhaseDiagramSpecification(
        _axis("voltage"), _axis("temperature"), {}
    ).evaluate(model, (dataset,))

    assert result.independent_components == ("C0", "C1")
    assert result.surface_grand_potential_ev_per_angstrom2.shape == (3, 2, 2)
    assert result.stable_phase_mask.shape == (3, 2, 2)
    assert np.all(result.stable_phase_mask[0])
    assert np.all(result.stable_phase_mask[1])
    assert not np.any(result.stable_phase_mask[2])
    np.testing.assert_array_equal(result.representative_phase_indices, 0)
    assert result.stable_phase_ids_at(0, 0) == (
        "dataset:p0",
        "dataset:p1",
    )
    assert not result.surface_grand_potential_ev_per_angstrom2.flags.writeable
    assert not result.stable_phase_mask.flags.writeable


def test_multiple_datasets_preserve_input_and_phase_order():
    """Qualified identities should make multi-dataset order deterministic."""
    model = _model(("voltage", "temperature"))
    first = _dataset(model.components, "first", (0.0,))
    second = _dataset(model.components, "second", (1.0,))

    result = PhaseDiagramSpecification(
        _axis("voltage"), _axis("temperature"), {}
    ).evaluate(model, (first, second))

    assert result.datasets == (first, second)
    assert result.phase_ids == ("first:p0", "second:p0")
    assert len(result.dataset_results) == 2


def test_multiple_datasets_require_unique_ids():
    """Qualified phase identity is ambiguous when dataset IDs repeat."""
    model = _model(("voltage", "temperature"))
    first = _dataset(model.components, "duplicate", (0.0,))
    second = _dataset(model.components, "duplicate", (1.0,))
    specification = PhaseDiagramSpecification(
        _axis("voltage"), _axis("temperature"), {}
    )

    with pytest.raises(ValueError, match="dataset_id.*unique"):
        specification.evaluate(model, (first, second))


def test_aligned_dataset_requires_its_reference_in_same_diagram():
    """Aligned target energies should be compared in their declared gauge."""
    model = _model(("voltage", "temperature"))
    components = model.components
    reference = _dataset(components, "reference", (0.0,))
    target_composition = {
        components[0]: 2,
        components[1]: 1,
    }
    target = PhaseDataset(
        "target",
        components,
        (Phase("p0", target_composition, 3.0, 2.0, 1),),
        _METHOD,
    )
    bulk = ReferencePhase(
        "bulk",
        {components[0]: 1, components[1]: 0},
        2.0,
        _METHOD,
    )
    aligned = DatasetAlignment(
        reference, target, "p0", "p0", bulk
    ).create_aligned_dataset()
    specification = PhaseDiagramSpecification(
        _axis("voltage"), _axis("temperature"), {}
    )

    with pytest.raises(ValueError, match="reference.*reference"):
        specification.evaluate(model, (aligned,))

    result = specification.evaluate(model, (reference, aligned))
    assert result.phase_ids == ("reference:p0", "target:p0")
