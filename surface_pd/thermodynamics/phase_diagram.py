"""Two-dimensional numerical phase-diagram evaluation."""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from types import MappingProxyType

import numpy as np

from ._validation import (
    validate_finite_real,
    validate_provenance,
    validate_variable_name,
)
from .alignment import AlignedPhaseDataset
from .grand_potential import GrandPotentialModel, GrandPotentialResult
from .phase_data import PhaseDataset
from .state import ThermodynamicState

_STABILITY_ABSOLUTE_TOLERANCE_EV_PER_ANGSTROM2 = 1e-10


def _read_only_array(
    value: object, dtype: type, field_name: str
) -> np.ndarray:
    """Return an owned, read-only NumPy array."""
    try:
        result = np.array(value, dtype=dtype, copy=True)
    except (TypeError, ValueError) as error:
        raise TypeError(
            f"{field_name} must be numerical array data"
        ) from error
    result.setflags(write=False)
    return result


@dataclass(frozen=True)
class DiagramAxis:
    """Define one numerical and display axis of a phase diagram.

    Parameters
    ----------
    state_variable : str
        Identifier of the thermodynamic state variable varied on this axis.
    values : array-like
        At least two finite, strictly increasing one-dimensional coordinates.
    label : str
        Nonempty, single-line display label without units.
    unit : str
        Nonempty, single-line display unit.
    """

    state_variable: str
    values: np.ndarray
    label: str
    unit: str

    def __post_init__(self) -> None:
        """Validate identity and display metadata and own the coordinates."""
        variable = validate_variable_name(self.state_variable)
        values = np.asarray(self.values)
        if (
            values.ndim != 1
            or values.size < 2
            or np.issubdtype(values.dtype, np.bool_)
            or np.issubdtype(values.dtype, np.complexfloating)
            or not np.issubdtype(values.dtype, np.number)
        ):
            raise ValueError(
                "values must be finite, one-dimensional numerical coordinates "
                "with at least two points"
            )
        values = np.array(values, dtype=float, copy=True)
        if not np.isfinite(values).all() or not np.all(np.diff(values) > 0):
            raise ValueError("values must be finite and strictly increasing")
        values.setflags(write=False)

        object.__setattr__(self, "state_variable", variable)
        object.__setattr__(self, "values", values)
        object.__setattr__(
            self, "label", validate_provenance(self.label, "label")
        )
        object.__setattr__(
            self, "unit", validate_provenance(self.unit, "unit")
        )


@dataclass(frozen=True)
class PhaseDiagramResult:
    """Store a complete, tie-preserving numerical phase diagram.

    Parameters
    ----------
    specification : PhaseDiagramSpecification
        Axes and fixed conditions used for evaluation.
    datasets : sequence of PhaseDataset or AlignedPhaseDataset
        Evaluated datasets in declared order.
    dataset_results : sequence of GrandPotentialResult
        Per-dataset results retaining resolved chemical potentials and every
        grand-potential normalization.
    independent_components : sequence of str
        Ordered components whose chemical potentials are defined directly by
        the model rather than solved from reference constraints.
    phase_ids : sequence of str
        Qualified phase identities in energy-tensor order.
    surface_grand_potential_ev_per_angstrom2 : array-like
        Complete intensive energy tensor with shape
        ``(number_of_phases, number_of_y_values, number_of_x_values)``.
    stable_phase_mask : array-like of bool
        Co-stability mask with the same shape as the energy tensor.
    representative_phase_indices : array-like of int
        First stable phase index at every point, with shape
        ``(number_of_y_values, number_of_x_values)``.
    x_mesh, y_mesh : array-like
        Coordinate meshes with shape
        ``(number_of_y_values, number_of_x_values)``.
    stability_tolerance_ev_per_angstrom2 : float
        Absolute tie tolerance in eV per square angstrom; relative tolerance is
        always zero.
    """

    specification: "PhaseDiagramSpecification"
    datasets: Sequence[PhaseDataset | AlignedPhaseDataset]
    dataset_results: Sequence[GrandPotentialResult]
    independent_components: Sequence[str]
    phase_ids: Sequence[str]
    surface_grand_potential_ev_per_angstrom2: np.ndarray
    stable_phase_mask: np.ndarray
    representative_phase_indices: np.ndarray
    x_mesh: np.ndarray
    y_mesh: np.ndarray
    stability_tolerance_ev_per_angstrom2: float = (
        _STABILITY_ABSOLUTE_TOLERANCE_EV_PER_ANGSTROM2
    )

    def __post_init__(self) -> None:
        """Validate tensor relationships and own every numerical result."""
        if not isinstance(self.specification, PhaseDiagramSpecification):
            raise TypeError(
                "specification must be a PhaseDiagramSpecification"
            )
        datasets = tuple(self.datasets)
        dataset_results = tuple(self.dataset_results)
        independent_components = tuple(
            validate_variable_name(component)
            for component in self.independent_components
        )
        phase_ids = tuple(self.phase_ids)
        if not datasets or len(datasets) != len(dataset_results):
            raise ValueError(
                "datasets and dataset_results must have equal nonzero length"
            )
        if not phase_ids or len(set(phase_ids)) != len(phase_ids):
            raise ValueError("phase_ids must be nonempty and unique")
        if not independent_components or len(
            set(independent_components)
        ) != len(independent_components):
            raise ValueError(
                "independent_components must be nonempty and unique"
            )
        if any(
            component not in dataset.components
            for dataset in datasets
            for component in independent_components
        ):
            raise ValueError(
                "independent_components must be present in every dataset"
            )

        mesh_shape = (
            self.specification.y_axis.values.size,
            self.specification.x_axis.values.size,
        )
        energy_shape = (len(phase_ids), *mesh_shape)
        energies = _read_only_array(
            self.surface_grand_potential_ev_per_angstrom2,
            float,
            "surface_grand_potential_ev_per_angstrom2",
        )
        if energies.shape != energy_shape or not np.isfinite(energies).all():
            raise ValueError(
                "surface grand potentials must be finite and have shape "
                f"{energy_shape}"
            )
        stable_mask = _read_only_array(
            self.stable_phase_mask, bool, "stable_phase_mask"
        )
        if (
            stable_mask.shape != energy_shape
            or not stable_mask.any(axis=0).all()
        ):
            raise ValueError(
                "stable_phase_mask must identify at least one phase per point"
            )
        representatives = _read_only_array(
            self.representative_phase_indices,
            int,
            "representative_phase_indices",
        )
        if (
            representatives.shape != mesh_shape
            or np.any(representatives < 0)
            or np.any(representatives >= len(phase_ids))
        ):
            raise ValueError(
                "representative_phase_indices must contain valid phase indices"
            )
        x_mesh = _read_only_array(self.x_mesh, float, "x_mesh")
        y_mesh = _read_only_array(self.y_mesh, float, "y_mesh")
        if x_mesh.shape != mesh_shape or y_mesh.shape != mesh_shape:
            raise ValueError(f"coordinate meshes must have shape {mesh_shape}")
        tolerance = validate_finite_real(
            self.stability_tolerance_ev_per_angstrom2,
            "stability_tolerance_ev_per_angstrom2",
        )
        if tolerance < 0:
            raise ValueError(
                "stability_tolerance_ev_per_angstrom2 must be nonnegative"
            )

        object.__setattr__(self, "datasets", datasets)
        object.__setattr__(self, "dataset_results", dataset_results)
        object.__setattr__(
            self, "independent_components", independent_components
        )
        object.__setattr__(self, "phase_ids", phase_ids)
        object.__setattr__(
            self, "surface_grand_potential_ev_per_angstrom2", energies
        )
        object.__setattr__(self, "stable_phase_mask", stable_mask)
        object.__setattr__(
            self, "representative_phase_indices", representatives
        )
        object.__setattr__(self, "x_mesh", x_mesh)
        object.__setattr__(self, "y_mesh", y_mesh)
        object.__setattr__(
            self, "stability_tolerance_ev_per_angstrom2", tolerance
        )

    @property
    def mesh_shape(self) -> tuple[int, int]:
        """Return plotting shape ``(number_of_y, number_of_x)``."""
        return self.y_mesh.shape

    def stable_phase_ids_at(
        self, y_index: int, x_index: int
    ) -> tuple[str, ...]:
        """Return all co-stable phase identities at one mesh point.

        Parameters
        ----------
        y_index : int
            Index along the y-axis coordinates.
        x_index : int
            Index along the x-axis coordinates.
        """
        stable = self.stable_phase_mask[:, y_index, x_index]
        return tuple(
            phase_id
            for phase_id, is_stable in zip(
                self.phase_ids, stable, strict=True
            )
            if is_stable
        )


@dataclass(frozen=True)
class PhaseDiagramSpecification:
    """Define and evaluate a two-dimensional thermodynamic diagram.

    Parameters
    ----------
    x_axis : DiagramAxis
        Horizontal state coordinate.
    y_axis : DiagramAxis
        Vertical state coordinate.
    fixed_conditions : mapping of str to float
        Finite scalar values for every required state variable not used as an
        axis.
    """

    x_axis: DiagramAxis
    y_axis: DiagramAxis
    fixed_conditions: Mapping[str, float]

    def __post_init__(self) -> None:
        """Validate axes and defensively own scalar fixed conditions."""
        if not isinstance(self.x_axis, DiagramAxis) or not isinstance(
            self.y_axis, DiagramAxis
        ):
            raise TypeError("x_axis and y_axis must be DiagramAxis objects")
        if self.x_axis.state_variable == self.y_axis.state_variable:
            raise ValueError(
                "x_axis and y_axis state variables must be distinct"
            )
        if not isinstance(self.fixed_conditions, Mapping):
            raise TypeError("fixed_conditions must be a mapping")

        fixed = {}
        axis_variables = {
            self.x_axis.state_variable,
            self.y_axis.state_variable,
        }
        for variable, value in self.fixed_conditions.items():
            try:
                valid_variable = validate_variable_name(variable)
                fixed_value = validate_finite_real(
                    value, f"fixed_conditions[{variable!r}]"
                )
            except (TypeError, ValueError) as error:
                raise type(error)(
                    "fixed_conditions must contain finite scalar values: "
                    f"{error}"
                ) from error
            if valid_variable in axis_variables:
                raise ValueError(
                    "fixed_conditions must not duplicate an axis variable"
                )
            fixed[valid_variable] = fixed_value
        object.__setattr__(
            self, "fixed_conditions", MappingProxyType(fixed)
        )

    def evaluate(
        self,
        model: GrandPotentialModel,
        datasets: Sequence[PhaseDataset | AlignedPhaseDataset],
    ) -> PhaseDiagramResult:
        """Evaluate every phase and retain tolerance-aware co-stability.

        Parameters
        ----------
        model : GrandPotentialModel
            General thermodynamic model used for every dataset.
        datasets : sequence of PhaseDataset or AlignedPhaseDataset
            Nonempty ordered candidate datasets. IDs must be unique, and every
        aligned view must be accompanied by its ordinary reference dataset.

        Returns
        -------
        PhaseDiagramResult
            Complete energies, per-dataset results, and stable-phase data.
        """
        if not isinstance(model, GrandPotentialModel):
            raise TypeError("model must be a GrandPotentialModel")
        if isinstance(datasets, (str, bytes)) or not isinstance(
            datasets, Sequence
        ):
            raise TypeError("datasets must be a nonempty sequence")
        datasets = tuple(datasets)
        if not datasets:
            raise ValueError("datasets must be a nonempty sequence")
        if any(
            not isinstance(item, (PhaseDataset, AlignedPhaseDataset))
            for item in datasets
        ):
            raise TypeError(
                "datasets must contain PhaseDataset or AlignedPhaseDataset"
            )
        dataset_ids = [dataset.dataset_id for dataset in datasets]
        if len(set(dataset_ids)) != len(dataset_ids):
            raise ValueError("dataset_id values must be unique")
        ordinary_ids = {
            dataset.dataset_id
            for dataset in datasets
            if isinstance(dataset, PhaseDataset)
        }
        for dataset in datasets:
            if (
                isinstance(dataset, AlignedPhaseDataset)
                and dataset.reference_dataset_id not in ordinary_ids
            ):
                raise ValueError(
                    f"aligned dataset {dataset.dataset_id!r} requires "
                    "reference "
                    f"dataset {dataset.reference_dataset_id!r}"
                )

        supplied_variables = {
            self.x_axis.state_variable,
            self.y_axis.state_variable,
            *self.fixed_conditions,
        }
        required_variables = set(model.required_state_variables)
        missing = sorted(required_variables - supplied_variables)
        unused = sorted(supplied_variables - required_variables)
        if missing:
            raise ValueError(
                "diagram is missing required state variables: "
                + ", ".join(missing)
            )
        if unused:
            raise ValueError(
                "diagram contains unused state variables: "
                + ", ".join(unused)
            )

        x_mesh, y_mesh = np.meshgrid(
            self.x_axis.values,
            self.y_axis.values,
            indexing="xy",
        )
        state_values = dict(self.fixed_conditions)
        state_values[self.x_axis.state_variable] = x_mesh
        state_values[self.y_axis.state_variable] = y_mesh
        state = ThermodynamicState(state_values)

        dataset_results = tuple(
            model.evaluate(dataset, state) for dataset in datasets
        )
        phase_ids = tuple(
            phase_id
            for result in dataset_results
            for phase_id in result.phase_ids
        )
        surface_energies = np.concatenate(
            [
                result.surface_grand_potential_ev_per_angstrom2
                for result in dataset_results
            ],
            axis=0,
        )

        minimum_energy = np.min(surface_energies, axis=0)
        stable_mask = np.isclose(
            surface_energies,
            minimum_energy[np.newaxis, ...],
            rtol=0.0,
            atol=_STABILITY_ABSOLUTE_TOLERANCE_EV_PER_ANGSTROM2,
        )
        # argmax returns the first True value, preserving declared phase order
        # solely for deterministic single-color rendering.
        representatives = np.argmax(stable_mask, axis=0)
        return PhaseDiagramResult(
            specification=self,
            datasets=datasets,
            dataset_results=dataset_results,
            independent_components=tuple(
                model.independent_chemical_potentials
            ),
            phase_ids=phase_ids,
            surface_grand_potential_ev_per_angstrom2=surface_energies,
            stable_phase_mask=stable_mask,
            representative_phase_indices=representatives,
            x_mesh=x_mesh,
            y_mesh=y_mesh,
        )
