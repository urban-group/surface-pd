"""Surface-constrained systematic composition enumeration."""

import math
from collections.abc import Mapping
from dataclasses import dataclass
from numbers import Integral, Real
from types import MappingProxyType

import numpy as np
from pymatgen.core import DummySpecies
from pymatgen.core.periodic_table import get_el_sp
from pymatgen.transformations.advanced_transformations import (
    EnumerateStructureTransformation,
)
from pymatgen.transformations.standard_transformations import (
    SubstitutionTransformation,
)

from surface_pd.core.enumeration_slab import EnumerationSlab
from surface_pd.core.pre_check import PreCheck
from surface_pd.error import IncompatibleSymmError, NoInversionSymmetryError


@dataclass(frozen=True)
class SurfaceEnumerationMetadata:
    """Provenance for one finalized surface-enumeration result.

    Parameters
    ----------
    transformation_matrix : tuple of tuple of int
        Integer transformation from the parent lattice to the result lattice.
    area_multiplier : int
        In-plane surface-cell area multiplier.
    raw_candidate_rank : int
        Zero-based position in the ranked raw pymatgen result list.
    symmetric : bool
        Whether symmetric two-surface finalization was applied.
    """

    transformation_matrix: tuple[tuple[int, int, int], ...]
    area_multiplier: int
    raw_candidate_rank: int
    symmetric: bool


class SurfaceEnumerator:
    """Enumerate ordered surface structures at a defined composition.

    Parameters
    ----------
    replacements : dict
        Enumerated parent species mapped to replacement-species fractional
        occupancies. A sum below one represents vacancies, for example
        ``{"Li": {"Li": 0.5}}``. Multiple replacement species represent
        substitutional orderings.
    min_cell_size : int, default=1
        Minimum in-plane surface-cell area multiplier.
    max_cell_size : int, default=1
        Maximum in-plane surface-cell area multiplier.
    enum_precision_parameter : float, default=1e-5
        Positive finite-coordinate tolerance passed to enumlib.

    Notes
    -----
    Enumeration requires a working enumlib installation discoverable by
    pymatgen. Temporary marker species are managed internally. Enumlib ranks
    symmetrically distinct raw candidates. Accepted candidates retain that
    order, and surface-pd does not apply a second structure-matching or
    deduplication pass after surface filtering or symmetric finalization.
    Callers that require a different, tolerance-dependent equivalence policy
    can explicitly compare selected results with pymatgen's
    :class:`~pymatgen.analysis.structure_matcher.StructureMatcher`.
    """

    def __init__(
        self,
        replacements: dict,
        min_cell_size: int = 1,
        max_cell_size: int = 1,
        enum_precision_parameter: float = 1e-5,
    ):
        self.replacements = _validate_replacements(replacements)
        self.min_cell_size = _validate_positive_integer(
            "min_cell_size", min_cell_size
        )
        self.max_cell_size = _validate_positive_integer(
            "max_cell_size", max_cell_size
        )
        if self.min_cell_size > self.max_cell_size:
            raise ValueError("min_cell_size must not exceed max_cell_size")
        if (
            isinstance(enum_precision_parameter, bool)
            or not isinstance(enum_precision_parameter, Real)
            or not math.isfinite(enum_precision_parameter)
        ):
            raise TypeError(
                "enum_precision_parameter must be a finite real number"
            )
        if enum_precision_parameter <= 0:
            raise ValueError("enum_precision_parameter must be positive")
        self.enum_precision_parameter = float(enum_precision_parameter)

    def apply_enumeration(
        self, structure: EnumerationSlab, max_structures: int = 2000
    ):
        """Enumerate ordered in-plane derivatives of a parent slab.

        Parameters
        ----------
        structure : EnumerationSlab
            Configured parent slab. The slab is not mutated.
        max_structures : int, default=2000
            Positive limit on ranked raw candidates requested from pymatgen.
            Surface filtering and finalization happen afterward, so fewer
            structures may be returned. Enumeration is not repeated to fill
            this limit with finalized structures.

        Returns
        -------
        list of EnumerationSlab
            Finalized surface slabs whose supercells expand only in plane.
            Each result exposes immutable ``enumeration_metadata`` containing
            its transformation, area multiplier, raw rank, and symmetry mode.

        Raises
        ------
        RuntimeError
            If pymatgen cannot execute enumlib.
        """
        if not isinstance(structure, EnumerationSlab):
            raise TypeError("structure must be an EnumerationSlab")
        max_structures = _validate_positive_integer(
            "max_structures", max_structures
        )
        if structure.enumerated_species is None:
            raise ValueError("structure.enumerated_species must be configured")
        if set(self.replacements) != set(structure.enumerated_species):
            raise ValueError(
                "replacements must contain the same species as "
                "structure.enumerated_species"
            )
        analysis = structure.analyze()
        self._validate_realizable_composition(structure, analysis)
        if structure.symmetric:
            _validate_symmetric_selective_dynamics(structure, analysis)
            pre_check = PreCheck(structure, analysis=analysis)
            if not pre_check.relax_both_surfaces():
                raise IncompatibleSymmError
            if not pre_check.has_inversion_symmetry():
                raise NoInversionSymmetryError

        markers = [
            DummySpecies(f"Xsurface{i + 1}")
            for i in range(len(structure.enumerated_species))
        ]
        marked = structure._surface_substitute(markers, analysis=analysis)
        substitutions = {
            marker: dict(self.replacements[species])
            for marker, species in zip(
                markers, structure.enumerated_species, strict=True
            )
        }
        raw_results = _apply_raw_enumeration(
            marked,
            substitutions,
            self.min_cell_size,
            self.max_cell_size,
            self.enum_precision_parameter,
            max_structures,
        )

        results = []
        for raw_candidate_rank, result in enumerate(raw_results):
            candidate = result["structure"]
            transformation = _in_plane_transformation(
                structure,
                candidate,
                self.min_cell_size,
                self.max_cell_size,
            )
            if transformation is None:
                continue
            finalized = EnumerationSlab.from_structure(
                candidate,
                direction=structure.direction,
                layer_tolerance_angstrom=(
                    structure.layer_tolerance_angstrom
                ),
                enumerated_species=structure.enumerated_species,
                num_enumerated_layers=structure.num_enumerated_layers,
                symmetric=structure.symmetric,
            )
            _complete_selective_dynamics(
                finalized,
                analysis.fixed_region_bounds_angstrom,
            )
            if structure.symmetric:
                finalized = finalized.symmetrize_top_base()
                finalized = EnumerationSlab.from_structure(
                    finalized,
                    direction=structure.direction,
                    layer_tolerance_angstrom=(
                        structure.layer_tolerance_angstrom
                    ),
                    enumerated_species=structure.enumerated_species,
                    num_enumerated_layers=structure.num_enumerated_layers,
                    symmetric=True,
                )
            area_multiplier = abs(round(np.linalg.det(transformation)))
            finalized._set_enumeration_metadata(
                SurfaceEnumerationMetadata(
                    transformation_matrix=tuple(
                        tuple(int(value) for value in row)
                        for row in transformation
                    ),
                    area_multiplier=int(area_multiplier),
                    raw_candidate_rank=raw_candidate_rank,
                    symmetric=bool(structure.symmetric),
                )
            )
            results.append(finalized)
        return results

    def _validate_realizable_composition(self, structure, analysis):
        """Require one allowed area multiplier to realize all occupancies."""
        top_indices = structure._select_enumerated_site_indices(
            analysis.layers, only_top=True
        )
        for multiplier in range(self.min_cell_size, self.max_cell_size + 1):
            realizable = True
            for species, replacements in self.replacements.items():
                site_count = len(top_indices[species])
                for occupancy in replacements.values():
                    population = site_count * multiplier * occupancy
                    if not math.isclose(
                        population, round(population), abs_tol=1e-8
                    ):
                        realizable = False
                        break
                if not realizable:
                    break
            if realizable:
                return
        raise ValueError(
            "requested occupancies cannot be realized from the selected "
            f"surface sites with cell sizes {self.min_cell_size} through "
            f"{self.max_cell_size}: "
            + "; ".join(
                f"{species} has {len(top_indices[species])} selected sites "
                f"with occupancies {dict(replacements)}"
                for species, replacements in self.replacements.items()
            )
        )


def _validate_positive_integer(name, value):
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError(f"{name} must be an integer")
    if value <= 0:
        raise ValueError(f"{name} must be positive")
    return int(value)


def _validate_replacements(replacements):
    if not isinstance(replacements, Mapping):
        raise TypeError("replacements must be a mapping")
    if not replacements:
        raise ValueError("replacements must not be empty")
    validated = {}
    for parent, targets in replacements.items():
        if not isinstance(parent, str) or not parent.strip():
            raise ValueError("parent species must be nonempty strings")
        get_el_sp(parent)
        if not isinstance(targets, Mapping) or not targets:
            raise ValueError("each parent species needs replacement species")
        validated_targets = {}
        total = 0.0
        for target, occupancy in targets.items():
            if not isinstance(target, str) or not target.strip():
                raise ValueError(
                    "replacement species must be nonempty strings"
                )
            get_el_sp(target)
            if (
                isinstance(occupancy, bool)
                or not isinstance(occupancy, Real)
                or not math.isfinite(occupancy)
            ):
                raise TypeError("occupancies must be finite real numbers")
            if not 0 <= occupancy <= 1:
                raise ValueError("occupancies must be between zero and one")
            validated_targets[target] = float(occupancy)
            total += float(occupancy)
        if total > 1.0 + 1e-8:
            raise ValueError("replacement occupancies must not sum above one")
        validated[parent] = MappingProxyType(validated_targets)
    return MappingProxyType(validated)


def _apply_raw_enumeration(
    structure,
    substitutions,
    min_cell_size,
    max_cell_size,
    enum_precision_parameter,
    max_structures,
):
    """Apply pymatgen's unrestricted three-dimensional enumerator."""
    partial = SubstitutionTransformation(substitutions).apply_transformation(
        structure
    )
    transformation = EnumerateStructureTransformation(
        min_cell_size=min_cell_size,
        max_cell_size=max_cell_size,
        enum_precision_parameter=enum_precision_parameter,
    )
    return transformation.apply_transformation(
        partial, return_ranked_list=max_structures
    )


def _in_plane_transformation(parent, child, min_cell_size, max_cell_size):
    """Return a strict in-plane integer transformation, or ``None``."""
    matrix = child.lattice.matrix @ np.linalg.inv(parent.lattice.matrix)
    integer_matrix = np.rint(matrix)
    if not np.allclose(matrix, integer_matrix, atol=1e-6, rtol=0):
        return None
    integer_matrix = integer_matrix.astype(int)
    direction = parent.direction
    expected_normal = np.zeros(3, dtype=int)
    expected_normal[direction] = 1
    if not np.array_equal(integer_matrix[direction], expected_normal):
        return None
    if not np.array_equal(integer_matrix[:, direction], expected_normal):
        return None
    multiplier = abs(round(np.linalg.det(integer_matrix)))
    if not min_cell_size <= multiplier <= max_cell_size:
        return None
    return integer_matrix


def _complete_selective_dynamics(structure, fixed_region_bounds_angstrom):
    """Restore flags omitted from sites created by pymatgen enumeration."""
    if fixed_region_bounds_angstrom is not None:
        lower_limit, upper_limit = fixed_region_bounds_angstrom
    coordinates = structure._site_coordinates_angstrom()
    for index, site in enumerate(structure):
        flags = site.properties.get("selective_dynamics")
        if flags is not None:
            continue
        coordinate = coordinates[index]
        fixed = (
            fixed_region_bounds_angstrom is not None
            and lower_limit - 1e-8 <= coordinate <= upper_limit + 1e-8
        )
        site.properties["selective_dynamics"] = [not fixed] * 3


def _validate_symmetric_selective_dynamics(structure, analysis):
    """Validate complete flags and relaxed outer layers for both surfaces."""
    for index, site in enumerate(structure):
        flags = site.properties.get("selective_dynamics")
        if (
            flags is None
            or len(flags) != 3
            or any(not isinstance(value, (bool, np.bool_)) for value in flags)
        ):
            raise ValueError(
                "symmetric enumeration requires three Boolean "
                f"selective_dynamics flags on site {index}"
            )
    if not analysis.fixed_site_indices:
        raise ValueError(
            "symmetric enumeration requires a fixed region identified by "
            "selective_dynamics"
        )
    outer_indices = (
        analysis.layers[0].site_indices + analysis.layers[-1].site_indices
    )
    if any(
        not all(structure[index].properties["selective_dynamics"])
        for index in outer_indices
    ):
        raise IncompatibleSymmError
