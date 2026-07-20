"""User-provided reference energies for surface phase calculations."""

import math
from dataclasses import dataclass
from numbers import Real


@dataclass(frozen=True)
class ReferenceEnergies:
    """Scientific reference energies supplied with phase-diagram data.

    Parameters
    ----------
    method : str
        Free-text calculation-method provenance. The text is retained for
        auditing and is never interpreted as a lookup key.
    li_ev_per_atom : float
        Energy of the lithium reference in eV per atom.
    o2_raw_ev_per_molecule : float
        Raw energy of the oxygen reference in eV per O2 molecule.
    o2_correction_ev_per_molecule : float
        Explicit correction added to the raw oxygen energy, in eV per O2
        molecule. Use zero when no correction is applied.
    bulk_litmo2_ev_per_formula_unit : float
        Bulk LiTMO2 energy in eV per formula unit.

    Raises
    ------
    ValueError
        If ``method`` is empty or contains a newline, or if an energy is not a
        finite real number. Boolean values are not accepted as energies.
    """

    method: str
    li_ev_per_atom: float
    o2_raw_ev_per_molecule: float
    o2_correction_ev_per_molecule: float
    bulk_litmo2_ev_per_formula_unit: float

    def __post_init__(self) -> None:
        """Normalize and validate all supplied scientific metadata."""
        if not isinstance(self.method, str) or not self.method.strip():
            raise ValueError("method must not be empty")
        if "\n" in self.method or "\r" in self.method:
            raise ValueError("method must be a single line")
        object.__setattr__(self, "method", self.method.strip())

        energy_fields = (
            "li_ev_per_atom",
            "o2_raw_ev_per_molecule",
            "o2_correction_ev_per_molecule",
            "bulk_litmo2_ev_per_formula_unit",
        )
        for field_name in energy_fields:
            value = getattr(self, field_name)
            if (
                isinstance(value, bool)
                or not isinstance(value, Real)
                or not math.isfinite(value)
            ):
                raise ValueError(
                    f"{field_name} must be a finite real number"
                )
            object.__setattr__(self, field_name, float(value))

    @property
    def o2_ev_per_molecule(self) -> float:
        """Return the corrected oxygen energy in eV per O2 molecule."""
        return (
            self.o2_raw_ev_per_molecule
            + self.o2_correction_ev_per_molecule
        )
