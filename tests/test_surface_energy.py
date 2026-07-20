"""Tests for dimensionally consistent surface-energy calculations."""

import numpy as np
import pytest

from surface_pd.plot import SurfaceEnergy


def _surface_energy(temperature, functional="PBE+U", tm_species="Ni"):
    """Return a minimal calculator for oxygen-potential tests."""
    return SurfaceEnergy(
        V=0.0,
        T=temperature,
        nLi=1,
        nTM=1,
        nO=2,
        dft_energy=0.0,
        a=1.0,
        b=1.0,
        gamma=90.0,
        TM_species=tm_species,
        functional=functional,
    )


@pytest.mark.parametrize(
    ("functional", "temperature", "expected_ev_per_oxygen"),
    [
        ("PBE+U", 298.15, -4.522051912226518),
        ("SCAN+rVV10+U", 298.15, -6.275466912226517),
        ("PBE+U", 1000.0, -5.344828641479072),
    ],
)
def test_oxygen_chemical_potential_uses_consistent_units(
    functional,
    temperature,
    expected_ev_per_oxygen,
):
    """Only the molar thermal correction should be converted to eV."""
    result = _surface_energy(temperature, functional).g_oxygen()

    assert result == pytest.approx(expected_ev_per_oxygen)


def test_oxygen_chemical_potential_preserves_temperature_shape():
    """Array temperatures should produce chemical potentials of equal shape."""
    temperature = np.array([[298.15, 500.0], [750.0, 1000.0]])

    result = _surface_energy(temperature).g_oxygen()

    assert result.shape == temperature.shape


@pytest.mark.parametrize(
    "temperature",
    [0.0, -1.0, np.array([298.15, 0.0]), np.array([-1.0, 298.15])],
)
def test_oxygen_chemical_potential_rejects_nonpositive_temperature(
    temperature,
):
    """Absolute temperature must be positive for the entropy correction."""
    with pytest.raises(ValueError, match="Temperature must be positive"):
        _surface_energy(temperature).g_oxygen()


def test_surface_energy_uses_confirmed_reference_set():
    """Surface energy should use the matching Li, O2, and LiTMO2 values."""
    calculator = SurfaceEnergy(
        V=3.0,
        T=298.15,
        nLi=0,
        nTM=1,
        nO=1,
        dft_energy=-10.0,
        a=2.0,
        b=3.0,
        gamma=90.0,
        TM_species="Co",
        functional="PBE+U",
    )

    expected = (
        -10.0
        + (-1.89965 - 3.0)
        + (-4.522051912226518)
        - (-22.69242)
    ) / 12.0
    assert calculator.get_gibbs_free_energy() == pytest.approx(expected)


def test_surface_energy_rejects_unsupported_reference_combination():
    """Construction should fail before an unsupported calculation begins."""
    with pytest.raises(ValueError, match="No LiMnO2 reference energy"):
        _surface_energy(298.15, "SCAN+rVV10+U", tm_species="Mn")


def test_surface_energy_has_no_duplicate_lithium_reference():
    """Lithium energy should have one functional-dependent source of truth."""
    assert not hasattr(_surface_energy(298.15), "e_li")
