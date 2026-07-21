"""Constrained chemical potentials and surface grand potentials."""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from numbers import Integral
from types import MappingProxyType

import numpy as np

from ._validation import validate_identifier, validate_variable_name
from .alignment import AlignedPhaseDataset
from .chemical_potential import ChemicalPotentialModel
from .phase_data import PhaseDataset, ReferencePhase
from .state import ThermodynamicState


def _read_only_float_array(
    value: object,
    expected_shape: tuple[int, ...],
    field_name: str,
) -> np.ndarray:
    """Return a finite, owned, read-only float array of a required shape."""
    array = np.asarray(value)
    if (
        np.issubdtype(array.dtype, np.bool_)
        or np.issubdtype(array.dtype, np.complexfloating)
        or not np.issubdtype(array.dtype, np.number)
    ):
        raise TypeError(f"{field_name} must contain finite real numbers")
    result = np.array(array, dtype=float, copy=True)
    if result.shape != expected_shape or not np.isfinite(result).all():
        raise ValueError(
            f"{field_name} must be finite and have shape {expected_shape}"
        )
    result.setflags(write=False)
    return result


@dataclass(frozen=True)
class GrandPotentialResult:
    """Store immutable output from a grand-potential evaluation.

    Parameters
    ----------
    dataset_id : str
        Identifier of the evaluated phase dataset.
    phase_ids : sequence of str
        Ordered, dataset-qualified phase identifiers.
    chemical_potentials_ev : mapping of str to array-like
        Resolved chemical potentials in eV per component. Each array has
        ``state_shape``.
    total_grand_potential_ev : array-like
        Total grand potential in eV with shape
        ``(number_of_phases, *state_shape)``.
    grand_potential_ev_per_surface_cell : array-like
        Grand potential in eV per represented surface cell, with the same
        shape as ``total_grand_potential_ev``.
    surface_grand_potential_ev_per_angstrom2 : array-like
        Surface grand potential in eV per square angstrom, with the same shape
        as ``total_grand_potential_ev``.
    state_shape : sequence of int
        Common shape of the thermodynamic state variables.

    Notes
    -----
    Instances are normally returned by :meth:`GrandPotentialModel.evaluate`.
    All arrays and mappings are defensively copied and read-only.
    """

    dataset_id: str
    phase_ids: Sequence[str]
    chemical_potentials_ev: Mapping[str, np.ndarray]
    total_grand_potential_ev: np.ndarray
    grand_potential_ev_per_surface_cell: np.ndarray
    surface_grand_potential_ev_per_angstrom2: np.ndarray
    state_shape: Sequence[int]

    def __post_init__(self) -> None:
        """Validate result dimensions and take ownership of numerical data."""
        dataset_id = validate_identifier(self.dataset_id, "dataset_id")

        if isinstance(self.phase_ids, (str, bytes)) or not isinstance(
            self.phase_ids, Sequence
        ):
            raise TypeError("phase_ids must be a nonempty sequence")
        phase_ids = tuple(self.phase_ids)
        if not phase_ids:
            raise ValueError("phase_ids must be a nonempty sequence")
        if any(
            not isinstance(phase_id, str)
            or not phase_id.strip()
            or "\n" in phase_id
            or "\r" in phase_id
            for phase_id in phase_ids
        ):
            raise ValueError(
                "phase_ids must contain nonempty single-line text"
            )
        if len(set(phase_ids)) != len(phase_ids):
            raise ValueError("phase_ids must be unique")

        if isinstance(self.state_shape, (str, bytes)) or not isinstance(
            self.state_shape, Sequence
        ):
            raise TypeError("state_shape must be a sequence of dimensions")
        state_shape = tuple(self.state_shape)
        if any(
            isinstance(size, bool)
            or not isinstance(size, Integral)
            or size < 0
            for size in state_shape
        ):
            raise ValueError(
                "state_shape dimensions must be nonnegative integers"
            )
        state_shape = tuple(int(size) for size in state_shape)

        if not isinstance(self.chemical_potentials_ev, Mapping):
            raise TypeError("chemical_potentials_ev must be a mapping")
        chemical_potentials = {}
        for component, values in self.chemical_potentials_ev.items():
            valid_component = validate_variable_name(component)
            chemical_potentials[valid_component] = _read_only_float_array(
                values,
                state_shape,
                f"chemical potential for {valid_component}",
            )
        if not chemical_potentials:
            raise ValueError("chemical_potentials_ev must not be empty")

        energy_shape = (len(phase_ids), *state_shape)
        object.__setattr__(self, "dataset_id", dataset_id)
        object.__setattr__(self, "phase_ids", phase_ids)
        object.__setattr__(self, "state_shape", state_shape)
        object.__setattr__(
            self,
            "chemical_potentials_ev",
            MappingProxyType(chemical_potentials),
        )
        for field_name in (
            "total_grand_potential_ev",
            "grand_potential_ev_per_surface_cell",
            "surface_grand_potential_ev_per_angstrom2",
        ):
            object.__setattr__(
                self,
                field_name,
                _read_only_float_array(
                    getattr(self, field_name), energy_shape, field_name
                ),
            )


@dataclass(frozen=True)
class GrandPotentialModel:
    r"""Resolve constrained chemical potentials and evaluate surface phases.

    Parameters
    ----------
    components : sequence of str
        Unique component names in deterministic matrix-column order.
    independent_chemical_potentials : mapping of str to ChemicalPotentialModel
        Models for user-selected independent components. Every other component
        is dependent and is solved from ``reference_phases``.
    reference_phases : sequence of ReferencePhase
        Constant reference equalities used to solve dependent potentials.

    Raises
    ------
    ValueError
        If the dependent system is not square and full-rank, component names
        are incompatible, or reference calculation methods differ.

    Notes
    -----
    If :math:`C_D` and :math:`C_I` contain the dependent and independent
    columns of the reference-composition matrix, respectively, this class
    solves

    .. math::

        C_D\mu_D(\mathbf{x}) =
        \mathbf{g}_\mathrm{ref} - C_I\mu_I(\mathbf{x}).

    The independent/dependent partition is explicit; the class does not choose
    it automatically and does not use a least-squares solution.
    """

    components: Sequence[str]
    independent_chemical_potentials: Mapping[str, ChemicalPotentialModel]
    reference_phases: Sequence[ReferencePhase]
    _dependent_components: tuple[str, ...] = field(init=False, repr=False)
    _dependent_matrix: np.ndarray = field(init=False, repr=False)
    _independent_matrix: np.ndarray = field(init=False, repr=False)
    _reference_energies: np.ndarray = field(init=False, repr=False)

    def __post_init__(self) -> None:
        """Validate and assemble the constant reference matrices."""
        components = self._validate_components(self.components)
        independent = self._validate_independent_models(
            self.independent_chemical_potentials, components
        )
        references = self._validate_references(
            self.reference_phases, components
        )
        dependent = tuple(
            component
            for component in components
            if component not in independent
        )
        if len(references) != len(dependent):
            raise ValueError(
                "number of reference phases must equal number of dependent "
                "components"
            )

        # C follows the declared component order, making the partition and
        # every subsequent result deterministic and easy to inspect.
        reference_matrix = np.array(
            [
                [
                    reference.composition.get(component, 0)
                    for component in components
                ]
                for reference in references
            ],
            dtype=float,
        ).reshape(len(references), len(components))
        dependent_indices = [components.index(name) for name in dependent]
        independent_indices = [components.index(name) for name in independent]
        dependent_matrix = reference_matrix[:, dependent_indices]
        independent_matrix = reference_matrix[:, independent_indices]

        if (
            dependent
            and np.linalg.matrix_rank(dependent_matrix) != len(dependent)
        ):
            raise ValueError("dependent reference system is singular")

        reference_energies = np.array(
            [reference.energy_ev_per_formula_unit for reference in references],
            dtype=float,
        )
        for array in (
            dependent_matrix,
            independent_matrix,
            reference_energies,
        ):
            array.setflags(write=False)

        object.__setattr__(self, "components", components)
        object.__setattr__(
            self,
            "independent_chemical_potentials",
            MappingProxyType(independent),
        )
        object.__setattr__(self, "reference_phases", references)
        object.__setattr__(self, "_dependent_components", dependent)
        object.__setattr__(self, "_dependent_matrix", dependent_matrix)
        object.__setattr__(self, "_independent_matrix", independent_matrix)
        object.__setattr__(self, "_reference_energies", reference_energies)

    @staticmethod
    def _validate_components(components: object) -> tuple[str, ...]:
        """Return a nonempty, unique, ordered component basis."""
        if isinstance(components, (str, bytes)) or not isinstance(
            components, Sequence
        ):
            raise TypeError("components must be a nonempty sequence")
        result = tuple(validate_variable_name(name) for name in components)
        if not result:
            raise ValueError("components must be a nonempty sequence")
        if len(set(result)) != len(result):
            raise ValueError("components must be unique")
        return result

    @staticmethod
    def _validate_independent_models(
        models: object, components: tuple[str, ...]
    ) -> dict[str, ChemicalPotentialModel]:
        """Return independent models ordered by the component basis."""
        if not isinstance(models, Mapping):
            raise TypeError(
                "independent_chemical_potentials must be a mapping"
            )
        unknown = set(models) - set(components)
        if unknown:
            raise ValueError(
                "independent components outside model basis: "
                f"{sorted(unknown)!r}"
            )
        result = {}
        for component in components:
            if component not in models:
                continue
            model = models[component]
            if not isinstance(model, ChemicalPotentialModel):
                raise TypeError(
                    f"chemical-potential model for {component!r} does not "
                    "satisfy ChemicalPotentialModel"
                )
            result[component] = model
        return result

    @staticmethod
    def _validate_references(
        references: object, components: tuple[str, ...]
    ) -> tuple[ReferencePhase, ...]:
        """Return compatible reference equalities with stable ordering."""
        if isinstance(references, (str, bytes)) or not isinstance(
            references, Sequence
        ):
            raise TypeError("reference_phases must be a sequence")
        result = tuple(references)
        if any(not isinstance(item, ReferencePhase) for item in result):
            raise TypeError(
                "reference_phases must contain only ReferencePhase objects"
            )
        reference_ids = [item.reference_id for item in result]
        if len(set(reference_ids)) != len(reference_ids):
            raise ValueError("reference phase identifiers must be unique")
        component_set = set(components)
        for reference in result:
            if not set(reference.composition).issubset(component_set):
                raise ValueError(
                    f"reference {reference.reference_id!r} contains "
                    "components outside model basis"
                )
        methods = {reference.calculation_method for reference in result}
        if len(methods) > 1:
            raise ValueError(
                "reference phases must use one calculation_method"
            )
        return result

    @property
    def dependent_components(self) -> tuple[str, ...]:
        """Return components solved from reference equalities."""
        return self._dependent_components

    @property
    def required_state_variables(self) -> frozenset[str]:
        """Return the union required by all independent reservoir models."""
        required = set()
        for model in self.independent_chemical_potentials.values():
            required.update(model.required_state_variables)
        return frozenset(required)

    def chemical_potentials(
        self, state: ThermodynamicState
    ) -> Mapping[str, np.ndarray]:
        """Resolve every component chemical potential over a state.

        Parameters
        ----------
        state : ThermodynamicState
            Named scalar or array-valued thermodynamic conditions.

        Returns
        -------
        mapping of str to numpy.ndarray
            Read-only potentials in component order and eV per component. Each
            array has ``state.shape``.
        """
        if not isinstance(state, ThermodynamicState):
            raise TypeError("state must be a ThermodynamicState")

        independent_values = {}
        for component, model in self.independent_chemical_potentials.items():
            independent_values[component] = _read_only_float_array(
                model.evaluate(state),
                state.shape,
                f"chemical potential for {component}",
            )

        dependent_values = self._solve_dependent_potentials(
            independent_values, state.shape
        )
        ordered = {}
        for component in self.components:
            values = independent_values.get(component)
            if values is None:
                values = dependent_values[component]
            ordered[component] = values
        return MappingProxyType(ordered)

    def _solve_dependent_potentials(
        self,
        independent_values: Mapping[str, np.ndarray],
        state_shape: tuple[int, ...],
    ) -> dict[str, np.ndarray]:
        """Solve ``C_D mu_D = g - C_I mu_I`` over the full state."""
        if not self.dependent_components:
            return {}

        if independent_values:
            independent_stack = np.stack(
                [
                    independent_values[name]
                    for name in self.independent_chemical_potentials
                ]
            )
            independent_contribution = np.tensordot(
                self._independent_matrix,
                independent_stack,
                axes=(1, 0),
            )
        else:
            independent_contribution = np.zeros(
                (len(self.reference_phases), *state_shape)
            )
        energy_shape = (len(self.reference_phases),) + (1,) * len(state_shape)
        right_hand_side = self._reference_energies.reshape(
            energy_shape
        ) - independent_contribution

        # NumPy solves columns of right-hand sides, so flatten only the state
        # dimensions and restore them immediately after the exact solve.
        flat_solution = np.linalg.solve(
            self._dependent_matrix,
            right_hand_side.reshape(len(self.dependent_components), -1),
        )
        solution = flat_solution.reshape(
            (len(self.dependent_components), *state_shape)
        )
        dependent_values = {}
        for index, component in enumerate(self.dependent_components):
            values = np.array(solution[index], copy=True)
            values.setflags(write=False)
            dependent_values[component] = values
        return dependent_values

    def evaluate(
        self,
        dataset: PhaseDataset | AlignedPhaseDataset,
        state: ThermodynamicState,
    ) -> GrandPotentialResult:
        r"""Evaluate all phases without selecting stable phases.

        Parameters
        ----------
        dataset : PhaseDataset or AlignedPhaseDataset
            Absolute phase energies, compositions, areas, and numbers of
            represented surfaces.
        state : ThermodynamicState
            Named scalar or array-valued thermodynamic conditions.

        Returns
        -------
        GrandPotentialResult
            Complete phase-by-state grand-potential tensors.

        Notes
        -----
        For phase :math:`s`, the evaluated quantities are

        .. math::

            \Omega_s = E_s - \sum_i n_{s,i}\mu_i,

        .. math::

            \Omega_s^\mathrm{cell} =
            \Omega_s/N_{s,\mathrm{surfaces}},
            \qquad
            \gamma_s =
            \Omega_s/(N_{s,\mathrm{surfaces}} A_s).
        """
        if not isinstance(dataset, (PhaseDataset, AlignedPhaseDataset)):
            raise TypeError(
                "dataset must be a PhaseDataset or AlignedPhaseDataset"
            )
        if tuple(dataset.components) != tuple(self.components):
            raise ValueError("dataset component basis must match model order")
        if self.reference_phases and dataset.calculation_method != (
            self.reference_phases[0].calculation_method
        ):
            raise ValueError(
                "dataset and reference phases must use one calculation_method"
            )

        chemical_potentials = self.chemical_potentials(state)
        composition_matrix = np.array(
            [
                [phase.composition[component] for component in self.components]
                for phase in dataset.phases
            ],
            dtype=float,
        )
        potential_stack = np.stack(
            [chemical_potentials[component] for component in self.components]
        )
        reservoir_energy = np.tensordot(
            composition_matrix, potential_stack, axes=(1, 0)
        )

        phase_axis_shape = (len(dataset.phases),) + (1,) * len(state.shape)
        dft_energies = np.array(
            [phase.dft_energy_ev for phase in dataset.phases]
        ).reshape(phase_axis_shape)
        if isinstance(dataset, AlignedPhaseDataset):
            dft_energies = dft_energies + dataset.energy_offset_ev
        numbers_of_surfaces = np.array(
            [phase.number_of_surfaces for phase in dataset.phases]
        ).reshape(phase_axis_shape)
        areas = np.array(
            [phase.surface_area_angstrom2 for phase in dataset.phases]
        ).reshape(phase_axis_shape)

        total = dft_energies - reservoir_energy
        per_surface_cell = total / numbers_of_surfaces
        per_area = per_surface_cell / areas
        return GrandPotentialResult(
            dataset_id=dataset.dataset_id,
            phase_ids=dataset.qualified_phase_ids,
            chemical_potentials_ev=chemical_potentials,
            total_grand_potential_ev=total,
            grand_potential_ev_per_surface_cell=per_surface_cell,
            surface_grand_potential_ev_per_angstrom2=per_area,
            state_shape=state.shape,
        )
