"""Tests for generalized thermodynamic state and reservoir models."""

from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from surface_pd.thermodynamics import (
    ChemicalPotentialModel,
    ConstantChemicalPotential,
    DirectChemicalPotential,
    FixedPressureOxygenChemicalPotential,
    IntercalationChemicalPotential,
    ThermodynamicState,
)


def test_empty_thermodynamic_state_is_scalar():
    """An empty state should provide the scalar broadcast shape."""
    state = ThermodynamicState({})

    assert state.shape == ()
    assert dict(state) == {}


def test_thermodynamic_state_requires_mapping():
    """A positional sequence should not be accepted as named state data."""
    with pytest.raises(TypeError, match="mapping"):
        ThermodynamicState([("temperature", 298.15)])


def test_thermodynamic_state_broadcasts_and_owns_values():
    """State arrays should be compatible, copied, and read-only."""
    temperature = np.array([[298.15], [500.0]])
    voltage = np.array([[0.0, 1.0, 2.0]])

    state = ThermodynamicState(
        {"temperature": temperature, "voltage": voltage}
    )
    temperature[0, 0] = 999.0

    assert state.shape == (2, 3)
    assert len(state) == 2
    assert state["temperature"].shape == (2, 3)
    assert state["voltage"].shape == (2, 3)
    assert state["temperature"][0, 0] == pytest.approx(298.15)
    with pytest.raises(ValueError, match="read-only"):
        state["temperature"][0, 0] = 400.0


@pytest.mark.parametrize("name", ["", " temperature", "temperature ", "T-K"])
def test_thermodynamic_state_rejects_invalid_variable_names(name):
    """State-variable names should be unambiguous identifiers."""
    with pytest.raises(ValueError, match="variable name"):
        ThermodynamicState({name: 1.0})


@pytest.mark.parametrize(
    "value",
    [True, 1 + 2j, "298.15", np.nan, np.inf, np.array([1.0, np.nan])],
)
def test_thermodynamic_state_rejects_invalid_values(value):
    """State values should contain only finite real numbers."""
    with pytest.raises((TypeError, ValueError), match="finite real"):
        ThermodynamicState({"temperature": value})


def test_thermodynamic_state_rejects_incompatible_shapes():
    """Every state value should share one NumPy broadcast shape."""
    with pytest.raises(ValueError, match="broadcast"):
        ThermodynamicState(
            {
                "temperature": np.ones((2, 2)),
                "voltage": np.ones(3),
            }
        )


def test_constant_chemical_potential_follows_state_shape():
    """A fixed reservoir should broadcast over the full state mesh."""
    state = ThermodynamicState({"voltage": np.array([0.0, 1.0, 2.0])})
    model = ConstantChemicalPotential(value_ev_per_component=-3.5)

    result = model.evaluate(state)

    assert isinstance(model, ChemicalPotentialModel)
    assert model.required_state_variables == frozenset()
    assert result.shape == state.shape
    assert result.dtype == float
    assert np.array_equal(result, [-3.5, -3.5, -3.5])
    with pytest.raises(FrozenInstanceError):
        model.value_ev_per_component = 0.0


@pytest.mark.parametrize("value", [True, "-3.5", np.nan, np.inf])
def test_constant_chemical_potential_rejects_invalid_energy(value):
    """Constant chemical potentials should be finite real energies."""
    with pytest.raises((TypeError, ValueError), match="finite real"):
        ConstantChemicalPotential(value_ev_per_component=value)


def test_direct_chemical_potential_uses_named_offset():
    """A direct model should add its named offset to its reference."""
    state = ThermodynamicState(
        {
            "delta_mu_O": np.array([-1.0, 0.0, 0.5]),
            "delta_mu_H": np.array([0.1, 0.2, 0.3]),
        }
    )
    model = DirectChemicalPotential(
        state_variable="delta_mu_O",
        reference_energy_ev_per_component=-4.5,
    )

    result = model.evaluate(state)

    assert model.required_state_variables == frozenset({"delta_mu_O"})
    assert np.array_equal(result, [-5.5, -4.5, -4.0])


def test_direct_chemical_potential_reports_missing_offset():
    """Missing model inputs should identify the required state variable."""
    model = DirectChemicalPotential(
        state_variable="delta_mu_O",
        reference_energy_ev_per_component=0.0,
    )

    with pytest.raises(ValueError, match="delta_mu_O"):
        model.evaluate(ThermodynamicState({"voltage": 3.0}))


@pytest.mark.parametrize("state_variable", ["", "delta mu", " delta_mu"])
def test_direct_chemical_potential_rejects_invalid_variable_name(
    state_variable,
):
    """Direct-coordinate variable names should follow the state contract."""
    with pytest.raises(ValueError, match="variable name"):
        DirectChemicalPotential(
            state_variable=state_variable,
            reference_energy_ev_per_component=0.0,
        )


def test_intercalation_chemical_potential_supports_valence():
    """The voltage law should scale by electrons transferred per component."""
    state = ThermodynamicState({"voltage": np.array([0.0, 1.5, 3.0])})
    model = IntercalationChemicalPotential(
        reference_energy_ev_per_component=-1.0,
        electrons_per_component=2,
        voltage_reference="Mg/Mg2+",
    )

    result = model.evaluate(state)

    assert model.required_state_variables == frozenset({"voltage"})
    assert np.array_equal(result, [-1.0, -4.0, -7.0])


@pytest.mark.parametrize("electrons", [True, 0, -1, 1.5])
def test_intercalation_chemical_potential_rejects_invalid_electron_count(
    electrons,
):
    """The voltage law should require a positive integer electron count."""
    with pytest.raises((TypeError, ValueError), match="positive integer"):
        IntercalationChemicalPotential(
            reference_energy_ev_per_component=-2.0,
            electrons_per_component=electrons,
            voltage_reference="Li/Li+",
        )


@pytest.mark.parametrize("reference", ["", "  ", "Li/Li+\nSHE"])
def test_intercalation_chemical_potential_rejects_invalid_reference(reference):
    """Voltage-reference provenance should be nonempty and single-line."""
    with pytest.raises(ValueError, match="voltage_reference"):
        IntercalationChemicalPotential(
            reference_energy_ev_per_component=-2.0,
            electrons_per_component=1,
            voltage_reference=reference,
        )


@pytest.mark.parametrize(
    ("raw_o2", "correction", "temperature", "expected_ev_per_oxygen"),
    [
        (-9.86018, 1.36, 298.15, -4.522051912226518),
        (-12.00701, 0.0, 298.15, -6.275466912226517),
        (-9.86018, 1.36, 1000.0, -5.344828641479072),
    ],
)
def test_fixed_pressure_oxygen_model_reproduces_legacy_values(
    raw_o2,
    correction,
    temperature,
    expected_ev_per_oxygen,
):
    """The extracted approximation should reproduce corrected legacy values."""
    model = FixedPressureOxygenChemicalPotential(
        raw_o2_energy_ev_per_molecule=raw_o2,
        correction_ev_per_molecule=correction,
    )

    result = model.evaluate(ThermodynamicState({"temperature": temperature}))

    assert model.required_state_variables == frozenset({"temperature"})
    assert result.shape == ()
    assert result.item() == pytest.approx(expected_ev_per_oxygen)


@pytest.mark.parametrize(
    "temperature",
    [0.0, -1.0, np.array([298.15, 0.0]), np.array([-1.0, 298.15])],
)
def test_fixed_pressure_oxygen_model_rejects_nonpositive_temperature(
    temperature,
):
    """The analytical entropy expression requires positive temperature."""
    model = FixedPressureOxygenChemicalPotential(
        raw_o2_energy_ev_per_molecule=-9.86018,
        correction_ev_per_molecule=1.36,
    )

    with pytest.raises(ValueError, match="Temperature must be positive"):
        model.evaluate(ThermodynamicState({"temperature": temperature}))


def test_models_require_thermodynamic_state():
    """Model evaluation should not silently accept an arbitrary mapping."""
    model = ConstantChemicalPotential(value_ev_per_component=-1.0)

    with pytest.raises(TypeError, match="ThermodynamicState"):
        model.evaluate({"temperature": 298.15})
