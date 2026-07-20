"""Regression tests for the publication reference-energy dataset."""

import pytest

from surface_pd.plot._reference_energies import (
    _BULK_ENERGY_BY_FUNCTIONAL,
    _LI_ENERGY_BY_FUNCTIONAL,
    _O2_ENERGY_BY_FUNCTIONAL,
    _get_bulk_energy,
    _get_reference_energies,
)


def test_historic_elemental_reference_energies_are_exact():
    """Li and O2 values should match Xinhao Li's published calculations."""
    assert dict(_LI_ENERGY_BY_FUNCTIONAL) == {
        "PBE+U": -1.89965,
        "SCAN+rVV10+U": -2.33333,
        "r2SCAN+rVV10+U": -2.32338,
    }
    assert dict(_O2_ENERGY_BY_FUNCTIONAL) == {
        "PBE+U": -8.50018,
        "SCAN+rVV10+U": -12.00701,
        "r2SCAN+rVV10+U": -11.54833,
    }


def test_historic_bulk_reference_energies_are_exact():
    """Only supplied, non-placeholder LiTMO2 values should be supported."""
    assert {
        species: dict(energies)
        for species, energies in _BULK_ENERGY_BY_FUNCTIONAL.items()
    } == {
        "Ni": {
            "PBE+U": -19.92283375,
            "SCAN+rVV10+U": -36.8133525,
        },
        "Co": {
            "PBE+U": -22.69242,
            "SCAN+rVV10+U": -37.2001966667,
            "r2SCAN+rVV10+U": -32.5698933333,
        },
        "Mn": {
            "PBE+U": -26.319605,
            "r2SCAN+rVV10+U": -36.4717,
        },
    }


@pytest.mark.parametrize(
    ("species", "functional"),
    [("Mn", "SCAN+rVV10+U"), ("Ni", "r2SCAN+rVV10+U")],
)
def test_missing_or_placeholder_bulk_references_are_rejected(
    species,
    functional,
):
    """Incomplete protocols must not silently enter a calculation."""
    with pytest.raises(ValueError, match="No Li.*O2 reference energy"):
        _get_bulk_energy(species, functional)


@pytest.mark.parametrize("functional", ["PBE", "SCAN", "unknown"])
def test_shorthand_and_unknown_functionals_are_rejected(functional):
    """Functional labels must identify the complete calculation protocol."""
    with pytest.raises(ValueError, match="Unsupported DFT functional"):
        _get_reference_energies("Ni", functional)


def test_unknown_transition_metal_is_rejected():
    """Only transition metals with confirmed bulk references are valid."""
    with pytest.raises(ValueError, match="Unsupported transition metal"):
        _get_reference_energies("Fe", "PBE+U")
