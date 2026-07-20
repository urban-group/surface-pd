"""Tests for dimensionally consistent surface-energy calculations."""

import numpy as np
import pytest

from surface_pd.plot import ReferenceEnergies, SurfaceEnergy


def _references(
    *,
    li=-1.89965,
    o2_raw=-9.86018,
    o2_correction=1.36,
    bulk=-19.92283375,
):
    """Return an explicit reference set for surface-energy tests."""
    return ReferenceEnergies(
        method="test method; U_Ni=? eV",
        li_ev_per_atom=li,
        o2_raw_ev_per_molecule=o2_raw,
        o2_correction_ev_per_molecule=o2_correction,
        bulk_litmo2_ev_per_formula_unit=bulk,
    )


def _surface_energy(temperature, references=None):
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
        reference_energies=references or _references(),
    )


@pytest.mark.parametrize(
    ("references", "temperature", "expected_ev_per_oxygen"),
    [
        (_references(), 298.15, -4.522051912226518),
        (
            _references(
                li=-2.33333,
                o2_raw=-12.00701,
                o2_correction=0.0,
                bulk=-36.8133525,
            ),
            298.15,
            -6.275466912226517,
        ),
        (_references(), 1000.0, -5.344828641479072),
    ],
)
def test_oxygen_chemical_potential_uses_consistent_units(
    references,
    temperature,
    expected_ev_per_oxygen,
):
    """Only the molar thermal correction should be converted to eV."""
    result = _surface_energy(temperature, references).g_oxygen()

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
        reference_energies=_references(bulk=-22.69242),
    )

    expected = (
        -10.0
        + (-1.89965 - 3.0)
        + (-4.522051912226518)
        - (-22.69242)
    ) / 12.0
    assert calculator.get_gibbs_free_energy() == pytest.approx(expected)


def test_surface_energy_requires_reference_value_object():
    """Arbitrary method labels or mappings must not replace validated data."""
    with pytest.raises(TypeError, match="reference_energies"):
        _surface_energy(298.15, references="PBE+U")
