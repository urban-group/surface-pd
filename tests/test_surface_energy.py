"""Tests for dimensionally consistent surface-energy calculations."""

import numpy as np
import pytest

from surface_pd.plot import SurfaceEnergy


def _surface_energy(temperature, functional="PBE"):
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
        TM_species="Ni",
        functional=functional,
    )


@pytest.mark.parametrize(
    ("functional", "temperature", "expected_ev_per_oxygen"),
    [
        ("PBE", 298.15, -5.201961912226517),
        ("SCAN", 298.15, -5.496961912226517),
        ("PBE", 1000.0, -6.024738641479072),
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
