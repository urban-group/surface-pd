"""Immutable phase and reference data for grand-potential models."""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from types import MappingProxyType

from ._validation import (
    validate_composition,
    validate_finite_real,
    validate_identifier,
    validate_positive_integer,
    validate_provenance,
    validate_variable_name,
)


@dataclass(frozen=True)
class Phase:
    """Store one calculated surface phase in canonical internal units.

    Parameters
    ----------
    phase_id : str
        Dataset-local phase identifier. It must be nonempty, single-line, and
        contain no colon.
    composition : mapping of str to int
        Absolute, nonnegative elemental counts in the calculated cell. At
        least one count must be positive; counts are not normalized.
    dft_energy_ev : float
        Finite total DFT energy of the calculated cell in eV.
    surface_area_angstrom2 : float
        Positive surface-unit-cell area in square angstroms.
    number_of_surfaces : int
        Positive number of equivalent surfaces represented by the cell.
    """

    phase_id: str
    composition: Mapping[str, int]
    dft_energy_ev: float
    surface_area_angstrom2: float
    number_of_surfaces: int

    def __post_init__(self) -> None:
        """Validate inputs and take ownership of the composition."""
        object.__setattr__(
            self, "phase_id", validate_identifier(self.phase_id, "phase_id")
        )
        object.__setattr__(
            self, "composition", validate_composition(self.composition)
        )
        object.__setattr__(
            self,
            "dft_energy_ev",
            validate_finite_real(self.dft_energy_ev, "dft_energy_ev"),
        )
        area = validate_finite_real(
            self.surface_area_angstrom2, "surface_area_angstrom2"
        )
        if area <= 0:
            raise ValueError("surface_area_angstrom2 must be positive")
        object.__setattr__(self, "surface_area_angstrom2", area)
        object.__setattr__(
            self,
            "number_of_surfaces",
            validate_positive_integer(
                self.number_of_surfaces, "number_of_surfaces"
            ),
        )


@dataclass(frozen=True)
class ReferencePhase:
    """Store a bulk reference equality per formula unit.

    Parameters
    ----------
    reference_id : str
        Nonempty, single-line reference identifier containing no colon.
    composition : mapping of str to int
        Absolute, nonnegative elemental counts per formula unit.
    energy_ev_per_formula_unit : float
        Finite reference energy in eV per formula unit.
    calculation_method : str
        Nonempty, single-line free text identifying the calculation method.
    """

    reference_id: str
    composition: Mapping[str, int]
    energy_ev_per_formula_unit: float
    calculation_method: str

    def __post_init__(self) -> None:
        """Validate inputs and take ownership of the composition."""
        object.__setattr__(
            self,
            "reference_id",
            validate_identifier(self.reference_id, "reference_id"),
        )
        object.__setattr__(
            self, "composition", validate_composition(self.composition)
        )
        object.__setattr__(
            self,
            "energy_ev_per_formula_unit",
            validate_finite_real(
                self.energy_ev_per_formula_unit,
                "energy_ev_per_formula_unit",
            ),
        )
        object.__setattr__(
            self,
            "calculation_method",
            validate_provenance(
                self.calculation_method, "calculation_method"
            ),
        )


@dataclass(frozen=True)
class PhaseDataset:
    """Group compatible phases with an explicit ordered component basis.

    Parameters
    ----------
    dataset_id : str
        Nonempty, single-line dataset identifier containing no colon.
    components : sequence of str
        Nonempty ordered component basis. Names must be unique and
        identifier-like.
    phases : sequence of Phase
        Nonempty collection of phases with unique local identifiers. Each
        phase composition must exactly match ``components``.
    calculation_method : str
        Nonempty, single-line free text identifying the calculation method.
    """

    dataset_id: str
    components: Sequence[str]
    phases: Sequence[Phase]
    calculation_method: str
    _phase_by_id: Mapping[str, Phase] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        """Validate the dataset and own its component and phase sequences."""
        dataset_id = validate_identifier(self.dataset_id, "dataset_id")
        if isinstance(self.components, (str, bytes)) or not isinstance(
            self.components, Sequence
        ):
            raise TypeError("components must be a nonempty sequence")
        components = tuple(self.components)
        if not components:
            raise ValueError("components must be a nonempty sequence")
        try:
            components = tuple(
                validate_variable_name(component) for component in components
            )
        except ValueError as error:
            raise ValueError(
                "components must contain identifier-like strings"
            ) from error
        if len(set(components)) != len(components):
            raise ValueError("components must be unique")

        if isinstance(self.phases, (str, bytes)) or not isinstance(
            self.phases, Sequence
        ):
            raise TypeError("phases must be a nonempty sequence of Phase")
        phases = tuple(self.phases)
        if not phases:
            raise ValueError("dataset must contain at least one phase")
        phase_by_id = {}
        component_set = set(components)
        for phase in phases:
            if not isinstance(phase, Phase):
                raise TypeError("phases must contain only Phase objects")
            if phase.phase_id in phase_by_id:
                raise ValueError(f"duplicate phase_id {phase.phase_id!r}")
            if set(phase.composition) != component_set:
                raise ValueError(
                    f"phase {phase.phase_id!r} does not match component basis"
                )
            phase_by_id[phase.phase_id] = phase

        object.__setattr__(self, "dataset_id", dataset_id)
        object.__setattr__(self, "components", components)
        object.__setattr__(self, "phases", phases)
        object.__setattr__(
            self,
            "calculation_method",
            validate_provenance(
                self.calculation_method, "calculation_method"
            ),
        )
        object.__setattr__(
            self, "_phase_by_id", MappingProxyType(phase_by_id)
        )

    @property
    def qualified_phase_ids(self) -> tuple[str, ...]:
        """Return deterministic ``dataset_id:phase_id`` identifiers."""
        return tuple(
            self.qualified_phase_id(phase.phase_id) for phase in self.phases
        )

    def get_phase(self, phase_id: str) -> Phase:
        """Return a phase by its dataset-local identifier.

        Parameters
        ----------
        phase_id : str
            Dataset-local phase identifier.
        """
        return self._phase_by_id[phase_id]

    def qualified_phase_id(self, phase_id: str) -> str:
        """Return the globally unambiguous identifier for one phase.

        Parameters
        ----------
        phase_id : str
            Dataset-local identifier of a phase in this dataset.
        """
        phase = self.get_phase(phase_id)
        return f"{self.dataset_id}:{phase.phase_id}"
