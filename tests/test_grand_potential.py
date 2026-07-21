"""Tests for constrained chemical potentials and grand potentials."""

from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from surface_pd.thermodynamics import (
    DirectChemicalPotential,
    GrandPotentialModel,
    Phase,
    PhaseDataset,
    ReferencePhase,
    ThermodynamicState,
)

_METHOD = "test DFT method; U_X=? eV"


def _reference(reference_id, composition, energy, method=_METHOD):
    """Return a concise reference equality for solver tests."""
    return ReferencePhase(
        reference_id=reference_id,
        composition=composition,
        energy_ev_per_formula_unit=energy,
        calculation_method=method,
    )


def _analytical_model():
    """Return a model where mu_B = (7 - mu_A) / 2."""
    return GrandPotentialModel(
        components=("A", "B"),
        independent_chemical_potentials={
            "A": DirectChemicalPotential("delta_mu_a", 1.0)
        },
        reference_phases=(
            _reference("AB2", {"A": 1, "B": 2}, 7.0),
        ),
    )


def _analytical_dataset(method=_METHOD):
    """Return two phases with hand-computable grand potentials."""
    return PhaseDataset(
        dataset_id="analytical",
        components=("A", "B"),
        phases=(
            Phase("p1", {"A": 1, "B": 1}, 10.0, 2.0, 1),
            Phase("p2", {"A": 2, "B": 2}, 20.0, 5.0, 2),
        ),
        calculation_method=method,
    )


def test_model_solves_dependent_chemical_potential_over_state():
    """The implementation should directly solve the reference equality."""
    model = _analytical_model()
    state = ThermodynamicState({"delta_mu_a": [0.0, 2.0]})

    chemical_potentials = model.chemical_potentials(state)

    assert tuple(chemical_potentials) == ("A", "B")
    np.testing.assert_allclose(chemical_potentials["A"], [1.0, 3.0])
    np.testing.assert_allclose(chemical_potentials["B"], [3.0, 2.0])
    assert model.dependent_components == ("B",)
    assert model.required_state_variables == frozenset({"delta_mu_a"})
    with pytest.raises(ValueError, match="read-only"):
        chemical_potentials["A"][0] = 0.0
    with pytest.raises(TypeError):
        chemical_potentials["A"] = np.array([0.0, 0.0])


def test_model_evaluates_total_cell_and_area_grand_potentials():
    """Every output normalization should follow its documented equation."""
    state = ThermodynamicState({"delta_mu_a": [0.0, 2.0]})

    result = _analytical_model().evaluate(_analytical_dataset(), state)

    assert result.dataset_id == "analytical"
    assert result.phase_ids == ("analytical:p1", "analytical:p2")
    assert result.state_shape == (2,)
    np.testing.assert_allclose(
        result.total_grand_potential_ev,
        [[6.0, 5.0], [12.0, 10.0]],
    )
    np.testing.assert_allclose(
        result.grand_potential_ev_per_surface_cell,
        [[6.0, 5.0], [6.0, 5.0]],
    )
    np.testing.assert_allclose(
        result.surface_grand_potential_ev_per_angstrom2,
        [[3.0, 2.5], [1.2, 1.0]],
    )
    assert not result.total_grand_potential_ev.flags.writeable
    assert not result.chemical_potentials_ev["B"].flags.writeable
    with pytest.raises(FrozenInstanceError):
        result.dataset_id = "changed"


def test_scalar_state_produces_one_value_per_phase():
    """Scalar state variables should not add a trailing array dimension."""
    result = _analytical_model().evaluate(
        _analytical_dataset(),
        ThermodynamicState({"delta_mu_a": 0.0}),
    )

    assert result.state_shape == ()
    assert result.total_grand_potential_ev.shape == (2,)


def test_all_components_may_be_independent_without_references():
    """A reference-free model is valid when every potential is supplied."""
    model = GrandPotentialModel(
        components=("A",),
        independent_chemical_potentials={
            "A": DirectChemicalPotential("mu_a", 0.0)
        },
        reference_phases=(),
    )

    potentials = model.chemical_potentials(
        ThermodynamicState({"mu_a": -2.0})
    )

    assert potentials["A"] == pytest.approx(-2.0)
    assert model.dependent_components == ()


def test_all_components_may_be_solved_from_references():
    """A fully constrained system needs no independent reservoir model."""
    model = GrandPotentialModel(
        components=("A", "B"),
        independent_chemical_potentials={},
        reference_phases=(
            _reference("A", {"A": 1}, 2.0),
            _reference("B", {"B": 1}, 3.0),
        ),
    )

    potentials = model.chemical_potentials(ThermodynamicState({}))

    assert potentials["A"] == pytest.approx(2.0)
    assert potentials["B"] == pytest.approx(3.0)


@pytest.mark.parametrize(
    ("components", "independent", "references", "message"),
    [
        (
            ("A", "B", "C"),
            {"A": DirectChemicalPotential("mu_a", 0.0)},
            (_reference("ABC", {"A": 1, "B": 1, "C": 1}, 0.0),),
            "number of reference phases",
        ),
        (
            ("A", "B"),
            {"A": DirectChemicalPotential("mu_a", 0.0)},
            (
                _reference("AB", {"A": 1, "B": 1}, 0.0),
                _reference("B", {"B": 1}, 0.0),
            ),
            "number of reference phases",
        ),
        (
            ("A", "B", "C"),
            {"A": DirectChemicalPotential("mu_a", 0.0)},
            (
                _reference("AB", {"A": 1, "B": 1}, 0.0),
                _reference("A2B2", {"A": 2, "B": 2}, 0.0),
            ),
            "singular",
        ),
    ],
)
def test_model_rejects_unsolvable_dependent_systems(
    components, independent, references, message
):
    """Only square, full-rank dependent systems are supported initially."""
    with pytest.raises(ValueError, match=message):
        GrandPotentialModel(components, independent, references)


def test_model_rejects_unknown_components_and_reference_provenance_mixing():
    """Component and DFT reference conventions should be explicit."""
    with pytest.raises(ValueError, match="independent.*outside"):
        GrandPotentialModel(
            ("A",),
            {"B": DirectChemicalPotential("mu_b", 0.0)},
            (),
        )
    with pytest.raises(ValueError, match="reference.*outside"):
        GrandPotentialModel(
            ("A", "B"),
            {"A": DirectChemicalPotential("mu_a", 0.0)},
            (_reference("BC", {"B": 1, "C": 1}, 0.0),),
        )
    with pytest.raises(ValueError, match="calculation_method"):
        GrandPotentialModel(
            ("A", "B", "C"),
            {"A": DirectChemicalPotential("mu_a", 0.0)},
            (
                _reference("AB", {"A": 1, "B": 1}, 0.0),
                _reference(
                    "AC", {"A": 1, "C": 1}, 0.0, "other method"
                ),
            ),
        )


def test_model_rejects_incompatible_dataset_before_evaluation():
    """The evaluator should not silently reorder or mix physical inputs."""
    model = _analytical_model()
    state = ThermodynamicState({"delta_mu_a": 0.0})
    reversed_basis = PhaseDataset(
        "reversed",
        ("B", "A"),
        (Phase("p", {"A": 1, "B": 1}, 0.0, 1.0, 1),),
        _METHOD,
    )

    with pytest.raises(ValueError, match="component basis"):
        model.evaluate(reversed_basis, state)
    with pytest.raises(ValueError, match="calculation_method"):
        model.evaluate(_analytical_dataset("other method"), state)


class _InvalidResultModel:
    """Chemical-potential model returning invalid numerical data."""

    required_state_variables = frozenset()

    def __init__(self, result):
        self.result = result

    def evaluate(self, state):
        """Return the intentionally invalid test value."""
        return self.result


@pytest.mark.parametrize("result", [np.nan, [1.0, 2.0]])
def test_model_validates_independent_model_results(result):
    """Reservoir models must return finite arrays with the state shape."""
    model = GrandPotentialModel(("A",), {"A": _InvalidResultModel(result)}, ())

    with pytest.raises(ValueError, match="chemical potential.*A"):
        model.chemical_potentials(ThermodynamicState({}))
