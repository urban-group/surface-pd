"""Surface-constrained systematic composition enumeration."""

import numpy as np
from pymatgen.core import DummySpecies
from pymatgen.transformations.advanced_transformations import (
    EnumerateStructureTransformation,
)
from pymatgen.transformations.standard_transformations import (
    SubstitutionTransformation,
)

from surface_pd.core.enumeration_slab import EnumerationSlab


class EnumWithComposition:
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
    pymatgen. Temporary marker species are managed internally.
    """

    def __init__(
        self,
        replacements: dict,
        min_cell_size: int = 1,
        max_cell_size: int = 1,
        enum_precision_parameter: float = 1e-5,
    ):
        self.replacements = replacements
        self.min_cell_size = min_cell_size
        self.max_cell_size = max_cell_size
        self.enum_precision_parameter = enum_precision_parameter

    def apply_enumeration(
        self, structure: EnumerationSlab, max_structures: int = 2000
    ):
        """Enumerate ordered in-plane derivatives of a parent slab.

        Parameters
        ----------
        structure : EnumerationSlab
            Configured parent slab. The slab is not mutated.
        max_structures : int, default=2000
            Maximum number of ranked raw structures requested from pymatgen.

        Returns
        -------
        list of EnumerationSlab
            Finalized surface slabs whose supercells expand only in plane.

        Raises
        ------
        RuntimeError
            If pymatgen cannot execute enumlib.
        """
        if not isinstance(structure, EnumerationSlab):
            raise TypeError("structure must be an EnumerationSlab")
        if structure.enumerated_species is None:
            raise ValueError("structure.enumerated_species must be configured")
        if set(self.replacements) != set(structure.enumerated_species):
            raise ValueError(
                "replacements must contain the same species as "
                "structure.enumerated_species"
            )

        markers = [
            DummySpecies(f"Xsurface{i + 1}")
            for i in range(len(structure.enumerated_species))
        ]
        marked = structure._surface_substitute(markers)
        substitutions = {
            marker: self.replacements[species]
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
        for result in raw_results:
            candidate = result["structure"]
            if not _is_in_plane_derivative(
                structure,
                candidate,
                self.min_cell_size,
                self.max_cell_size,
            ):
                continue
            finalized = EnumerationSlab.from_structure(
                candidate,
                direction=structure.direction,
                tolerance=structure.tolerance,
                enumerated_species=structure.enumerated_species,
                num_enumerated_layers=structure.num_enumerated_layers,
                symmetric=structure.symmetric,
            )
            if structure.symmetric:
                finalized = finalized.symmetrize_top_base()
                finalized = EnumerationSlab.from_structure(
                    finalized,
                    direction=structure.direction,
                    tolerance=structure.tolerance,
                    enumerated_species=structure.enumerated_species,
                    num_enumerated_layers=structure.num_enumerated_layers,
                    symmetric=True,
                )
            results.append(finalized)
        return results


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


def _is_in_plane_derivative(parent, child, min_cell_size, max_cell_size):
    """Return whether a child lattice is a strict in-plane supercell."""
    matrix = child.lattice.matrix @ np.linalg.inv(parent.lattice.matrix)
    integer_matrix = np.rint(matrix)
    if not np.allclose(matrix, integer_matrix, atol=1e-6, rtol=0):
        return False
    integer_matrix = integer_matrix.astype(int)
    direction = parent.direction
    expected_normal = np.zeros(3, dtype=int)
    expected_normal[direction] = 1
    if not np.array_equal(integer_matrix[direction], expected_normal):
        return False
    if not np.array_equal(integer_matrix[:, direction], expected_normal):
        return False
    multiplier = abs(round(np.linalg.det(integer_matrix)))
    return min_cell_size <= multiplier <= max_cell_size
