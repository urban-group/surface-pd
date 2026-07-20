"""Confirmed DFT reference energies used by surface phase calculations.

The values were provided by Xinhao Li and restored from commit ``c0279c7``.
Energies are in eV: Li per atom, O2 per molecule, and bulk LiTMO2 per
LiTMO2 formula unit. Functional labels describe the complete computational
protocol and must not be shortened.
"""

from types import MappingProxyType

_SUPPORTED_FUNCTIONALS = (
    "PBE+U",
    "SCAN+rVV10+U",
    "r2SCAN+rVV10+U",
)

_LI_ENERGY_BY_FUNCTIONAL = MappingProxyType(
    {
        "PBE+U": -1.89965,
        "SCAN+rVV10+U": -2.33333,
        "r2SCAN+rVV10+U": -2.32338,
    }
)

# Xinhao's PBE+U O2 protocol applies a +1.36 eV correction to the raw
# isolated-molecule energy of -9.86018 eV.
_PBE_U_O2_RAW_ENERGY_EV = -9.86018
_PBE_U_O2_CORRECTION_EV = 1.36
_O2_ENERGY_BY_FUNCTIONAL = MappingProxyType(
    {
        "PBE+U": _PBE_U_O2_RAW_ENERGY_EV + _PBE_U_O2_CORRECTION_EV,
        "SCAN+rVV10+U": -12.00701,
        "r2SCAN+rVV10+U": -11.54833,
    }
)

# The historic Mn SCAN+rVV10+U value was zero (a placeholder), and no Ni
# r2SCAN+rVV10+U value was supplied. Neither combination is supported.
_BULK_ENERGY_BY_FUNCTIONAL = MappingProxyType(
    {
        "Ni": MappingProxyType(
            {
                "PBE+U": -19.92283375,
                "SCAN+rVV10+U": -36.8133525,
            }
        ),
        "Co": MappingProxyType(
            {
                "PBE+U": -22.69242,
                "SCAN+rVV10+U": -37.2001966667,
                "r2SCAN+rVV10+U": -32.5698933333,
            }
        ),
        "Mn": MappingProxyType(
            {
                "PBE+U": -26.319605,
                "r2SCAN+rVV10+U": -36.4717,
            }
        ),
    }
)


def _validate_functional(functional: str) -> None:
    """Raise ``ValueError`` unless *functional* is a complete protocol."""
    if functional not in _SUPPORTED_FUNCTIONALS:
        supported = ", ".join(_SUPPORTED_FUNCTIONALS)
        raise ValueError(
            f"Unsupported DFT functional {functional!r}. "
            f"Supported protocols are: {supported}."
        )


def _get_bulk_energy(tm_species: str, functional: str) -> float:
    """Return a confirmed bulk LiTMO2 energy in eV per formula unit."""
    _validate_functional(functional)
    if tm_species not in _BULK_ENERGY_BY_FUNCTIONAL:
        supported = ", ".join(_BULK_ENERGY_BY_FUNCTIONAL)
        raise ValueError(
            f"Unsupported transition metal {tm_species!r}. "
            f"Confirmed references are available for: {supported}."
        )
    try:
        return _BULK_ENERGY_BY_FUNCTIONAL[tm_species][functional]
    except KeyError as error:
        raise ValueError(
            f"No Li{tm_species}O2 reference energy is available for "
            f"{functional!r}."
        ) from error


def _get_reference_energies(
    tm_species: str,
    functional: str,
) -> tuple[float, float, float]:
    """Return Li, O2, and bulk energies for one confirmed protocol."""
    bulk_energy = _get_bulk_energy(tm_species, functional)
    return (
        _LI_ENERGY_BY_FUNCTIONAL[functional],
        _O2_ENERGY_BY_FUNCTIONAL[functional],
        bulk_energy,
    )
