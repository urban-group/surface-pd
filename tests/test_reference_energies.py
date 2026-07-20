"""Tests for user-supplied scientific reference energies."""

from dataclasses import FrozenInstanceError

import pytest

from surface_pd.plot import ReferenceEnergies


def test_reference_energies_expose_effective_oxygen_energy():
    """The auditable raw value and correction should produce the used value."""
    references = ReferenceEnergies(
        method="PBE+U; U_Ni=? eV",
        li_ev_per_atom=-1.89965,
        o2_raw_ev_per_molecule=-9.86018,
        o2_correction_ev_per_molecule=1.36,
        bulk_litmo2_ev_per_formula_unit=-19.92283375,
    )

    assert references.o2_ev_per_molecule == pytest.approx(-8.50018)


def test_reference_energies_are_immutable():
    """A validated scientific input must not change during calculation."""
    references = ReferenceEnergies(
        method="SCAN+rVV10+U; U_Ni=? eV",
        li_ev_per_atom=-2.33333,
        o2_raw_ev_per_molecule=-12.00701,
        o2_correction_ev_per_molecule=0.0,
        bulk_litmo2_ev_per_formula_unit=-36.8133525,
    )

    with pytest.raises(FrozenInstanceError):
        references.li_ev_per_atom = 0.0


@pytest.mark.parametrize("method", ["", "   "])
def test_reference_energies_reject_empty_method(method):
    """Method provenance must be present even when details are unknown."""
    with pytest.raises(ValueError, match="method must not be empty"):
        ReferenceEnergies(method, -1.0, -2.0, 0.0, -3.0)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("li_ev_per_atom", float("nan")),
        ("o2_raw_ev_per_molecule", float("inf")),
        ("o2_correction_ev_per_molecule", float("-inf")),
        ("bulk_litmo2_ev_per_formula_unit", True),
    ],
)
def test_reference_energies_reject_invalid_numbers(field, value):
    """All energy fields must be finite real numbers, excluding booleans."""
    values = {
        "method": "custom method",
        "li_ev_per_atom": -1.0,
        "o2_raw_ev_per_molecule": -2.0,
        "o2_correction_ev_per_molecule": 0.0,
        "bulk_litmo2_ev_per_formula_unit": -3.0,
    }
    values[field] = value

    with pytest.raises(ValueError, match=field):
        ReferenceEnergies(**values)
