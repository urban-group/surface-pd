"""Explicit energy alignment between compatible surface datasets."""

from dataclasses import dataclass, field

import numpy as np

from ._validation import validate_identifier
from .phase_data import Phase, PhaseDataset, ReferencePhase

_AREA_RELATIVE_TOLERANCE = 5e-3
_AREA_ABSOLUTE_TOLERANCE_ANGSTROM2 = 1e-8


@dataclass(frozen=True)
class DatasetAlignment:
    r"""Align a target dataset to an explicit reference-dataset anchor.

    Parameters
    ----------
    reference_dataset : PhaseDataset
        Root dataset whose energy convention remains unchanged.
    target_dataset : PhaseDataset
        Dataset receiving the calculated constant total-energy offset.
    reference_anchor_phase_id : str
        Local identifier of the explicit anchor in ``reference_dataset``.
    target_anchor_phase_id : str
        Local identifier of the explicit anchor in ``target_dataset``.
    bulk_reference : ReferencePhase
        Single bulk formula unit relating the two anchor compositions.

    Raises
    ------
    TypeError
        If datasets are not ordinary :class:`PhaseDataset` objects or the bulk
        reference has the wrong type.
    ValueError
        If identities, component bases, calculation methods, surface geometry,
        or anchor stoichiometry are incompatible.

    Notes
    -----
    The integer :math:`k` is defined by

    .. math::

        \mathbf{n}_\mathrm{target} - \mathbf{n}_\mathrm{reference}
        = k\mathbf{n}_\mathrm{bulk}.

    The total-energy offset added to every target phase is

    .. math::

        \Delta E_\mathrm{target} = E_\mathrm{reference}
        - E_\mathrm{target} + kG_\mathrm{bulk}.

    Alignment is direct to one reference dataset. Aligned views cannot be used
    as alignment inputs, so chained and cyclic alignment are unsupported.
    """

    reference_dataset: PhaseDataset
    target_dataset: PhaseDataset
    reference_anchor_phase_id: str
    target_anchor_phase_id: str
    bulk_reference: ReferencePhase
    _bulk_unit_count: int = field(init=False, repr=False)
    _energy_offset_ev: float = field(init=False, repr=False)

    def __post_init__(self) -> None:
        """Validate compatibility and calculate the target-energy offset."""
        if not isinstance(
            self.reference_dataset, PhaseDataset
        ) or not isinstance(
            self.target_dataset,
            PhaseDataset,
        ):
            raise TypeError(
                "alignment inputs must be ordinary PhaseDataset objects"
            )
        if not isinstance(self.bulk_reference, ReferencePhase):
            raise TypeError("bulk_reference must be a ReferencePhase")
        if (
            self.reference_dataset is self.target_dataset
            or self.reference_dataset.dataset_id
            == self.target_dataset.dataset_id
        ):
            raise ValueError("alignment requires two distinct datasets")

        reference_anchor_id = validate_identifier(
            self.reference_anchor_phase_id,
            "reference_anchor_phase_id",
        )
        target_anchor_id = validate_identifier(
            self.target_anchor_phase_id,
            "target_anchor_phase_id",
        )
        reference_anchor = self._get_anchor(
            self.reference_dataset,
            reference_anchor_id,
            "reference_anchor_phase_id",
        )
        target_anchor = self._get_anchor(
            self.target_dataset,
            target_anchor_id,
            "target_anchor_phase_id",
        )

        if self.reference_dataset.components != self.target_dataset.components:
            raise ValueError(
                "alignment datasets must use identical component bases"
            )
        components = self.reference_dataset.components
        if set(self.bulk_reference.composition) != set(components):
            raise ValueError(
                "bulk reference composition must match the dataset component "
                "basis"
            )
        method = self.reference_dataset.calculation_method
        if self.target_dataset.calculation_method != method:
            raise ValueError(
                "alignment datasets must use one calculation_method"
            )
        if self.bulk_reference.calculation_method != method:
            raise ValueError(
                "bulk reference calculation_method must match the datasets"
            )

        if reference_anchor.number_of_surfaces != (
            target_anchor.number_of_surfaces
        ):
            raise ValueError(
                "alignment anchor numbers of surfaces must match"
            )
        if not np.isclose(
            reference_anchor.surface_area_angstrom2,
            target_anchor.surface_area_angstrom2,
            rtol=_AREA_RELATIVE_TOLERANCE,
            atol=_AREA_ABSOLUTE_TOLERANCE_ANGSTROM2,
        ):
            raise ValueError(
                "alignment anchor surface areas must match within 0.5%"
            )

        composition_difference = tuple(
            target_anchor.composition[component]
            - reference_anchor.composition[component]
            for component in components
        )
        bulk_composition = tuple(
            self.bulk_reference.composition[component]
            for component in components
        )
        bulk_unit_count = self._integer_bulk_unit_count(
            composition_difference, bulk_composition
        )
        energy_offset = (
            reference_anchor.dft_energy_ev
            - target_anchor.dft_energy_ev
            + bulk_unit_count
            * self.bulk_reference.energy_ev_per_formula_unit
        )
        if not np.isfinite(energy_offset):
            raise ValueError(
                "calculated alignment energy offset must be finite"
            )

        object.__setattr__(
            self, "reference_anchor_phase_id", reference_anchor_id
        )
        object.__setattr__(self, "target_anchor_phase_id", target_anchor_id)
        object.__setattr__(self, "_bulk_unit_count", bulk_unit_count)
        object.__setattr__(self, "_energy_offset_ev", float(energy_offset))

    @staticmethod
    def _get_anchor(
        dataset: PhaseDataset,
        phase_id: str,
        field_name: str,
    ) -> Phase:
        """Return one explicitly named anchor with a field-specific error."""
        try:
            return dataset.get_phase(phase_id)
        except KeyError as error:
            raise ValueError(
                f"{field_name} {phase_id!r} is not present in dataset "
                f"{dataset.dataset_id!r}"
            ) from error

    @staticmethod
    def _integer_bulk_unit_count(
        difference: tuple[int, ...],
        bulk_composition: tuple[int, ...],
    ) -> int:
        """Return exact integer ``k`` in ``difference = k * bulk``."""
        candidate_counts = []
        for difference_count, bulk_count in zip(
            difference, bulk_composition, strict=True
        ):
            if bulk_count == 0:
                if difference_count != 0:
                    raise ValueError(
                        "anchor compositions must differ by an integer number "
                        "of bulk reference units"
                    )
                continue
            quotient, remainder = divmod(difference_count, bulk_count)
            if remainder:
                raise ValueError(
                    "anchor compositions must differ by an integer number of "
                    "bulk reference units"
                )
            candidate_counts.append(quotient)

        if not candidate_counts or len(set(candidate_counts)) != 1:
            raise ValueError(
                "anchor compositions must differ by an integer number of "
                "bulk reference units"
            )
        return candidate_counts[0]

    @property
    def bulk_unit_count(self) -> int:
        """Return signed integer bulk units in target minus reference."""
        return self._bulk_unit_count

    @property
    def energy_offset_ev(self) -> float:
        """Return the total-energy offset added to every target phase in eV."""
        return self._energy_offset_ev

    @property
    def reference_anchor_id(self) -> str:
        """Return the qualified reference-dataset anchor identity."""
        return self.reference_dataset.qualified_phase_id(
            self.reference_anchor_phase_id
        )

    @property
    def target_anchor_id(self) -> str:
        """Return the qualified target-dataset anchor identity."""
        return self.target_dataset.qualified_phase_id(
            self.target_anchor_phase_id
        )

    def create_aligned_dataset(self) -> "AlignedPhaseDataset":
        """Return an immutable target-energy view backed by this alignment."""
        return AlignedPhaseDataset(self)


@dataclass(frozen=True)
class AlignedPhaseDataset:
    """Expose a target dataset through one validated energy alignment.

    Parameters
    ----------
    alignment : DatasetAlignment
        Validated direct-to-reference alignment defining the source dataset,
        energy offset, anchors, bulk formula unit, and alignment direction.

    Notes
    -----
    Raw :class:`Phase` energies remain unchanged. A
    :class:`GrandPotentialModel` recognizes this view and adds
    ``energy_offset_ev`` immediately before grand-potential evaluation.
    """

    alignment: DatasetAlignment

    def __post_init__(self) -> None:
        """Require a validated alignment as the sole source of an offset."""
        if not isinstance(self.alignment, DatasetAlignment):
            raise TypeError("alignment must be a DatasetAlignment")

    @property
    def source_dataset(self) -> PhaseDataset:
        """Return the unchanged target dataset underlying this view."""
        return self.alignment.target_dataset

    @property
    def reference_dataset_id(self) -> str:
        """Return the dataset ID defining the reference energy convention."""
        return self.alignment.reference_dataset.dataset_id

    @property
    def energy_offset_ev(self) -> float:
        """Return the validated total-energy offset in eV."""
        return self.alignment.energy_offset_ev

    @property
    def dataset_id(self) -> str:
        """Return the source dataset identifier."""
        return self.source_dataset.dataset_id

    @property
    def components(self) -> tuple[str, ...]:
        """Return the source dataset's ordered component basis."""
        return self.source_dataset.components

    @property
    def phases(self) -> tuple[Phase, ...]:
        """Return the source dataset's immutable phase sequence."""
        return self.source_dataset.phases

    @property
    def calculation_method(self) -> str:
        """Return the source dataset's calculation-method provenance."""
        return self.source_dataset.calculation_method

    @property
    def qualified_phase_ids(self) -> tuple[str, ...]:
        """Return the source dataset's qualified phase identities."""
        return self.source_dataset.qualified_phase_ids

    def get_phase(self, phase_id: str) -> Phase:
        """Return a source phase by its local identifier.

        Parameters
        ----------
        phase_id : str
            Dataset-local phase identifier.
        """
        return self.source_dataset.get_phase(phase_id)

    def qualified_phase_id(self, phase_id: str) -> str:
        """Return the qualified identity of one source phase.

        Parameters
        ----------
        phase_id : str
            Dataset-local phase identifier.
        """
        return self.source_dataset.qualified_phase_id(phase_id)
