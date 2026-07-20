"""
Enumeration module for systematic surface composition generation.

This module provides the EnumWithComposition class that wraps pymatgen's
enumeration functionality to systematically generate surface structures
with different compositions.
"""

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
    subs_dict : dict
        Mapping accepted by pymatgen's ``SubstitutionTransformation``. Each
        key identifies a species in the parent slab and each value maps
        replacement species to fractional occupancies, for example
        ``{"X": {"Li": 0.5, "O": 0.5}}``.
    min_cell_size : int, default=1
        Minimum multiple of the parent cell passed to enumlib. Must be
        positive and no larger than ``max_cell_size``.
    max_cell_size : int, default=1
        Maximum multiple of the parent cell passed to enumlib. Must be
        positive and no smaller than ``min_cell_size``.
    enum_precision_parameter : float, default=1e-5
        Positive finite-coordinate tolerance passed to enumlib.

    Notes
    -----
    Enumeration requires a working enumlib installation discoverable by
    pymatgen. Validation of substitutions and cell-size bounds is delegated to
    pymatgen.
    """

    def __init__(
        self,
        subs_dict: dict,
        min_cell_size: int = 1,
        max_cell_size: int = 1,
        enum_precision_parameter: float = 1e-5,
    ):
        self.subs_dict = subs_dict
        self.min_cell_size = min_cell_size
        self.max_cell_size = max_cell_size
        self.enum_precision_parameter = enum_precision_parameter

    def apply_enumeration(
        self, structure: EnumerationSlab, max_structures: int = 2000
    ):
        """Enumerate ordered derivatives of a parent slab.

        Parameters
        ----------
        structure : EnumerationSlab
            Parent slab to substitute and enumerate. The slab is not mutated.
        max_structures : int, default=2000
            Maximum number of ranked structures requested from pymatgen. Must
            be positive.

        Returns
        -------
        list of dict
            Pymatgen ranked-result dictionaries. Each dictionary contains a
            ``structure`` entry and ranking metadata supplied by
            ``EnumerateStructureTransformation``.

        Raises
        ------
        RuntimeError
            If pymatgen cannot execute enumlib.

        Notes
        -----
        Validation errors raised by pymatgen are propagated unchanged.
        """
        subs = SubstitutionTransformation(self.subs_dict)
        surface_structure_partial = subs.apply_transformation(structure)
        enum = EnumerateStructureTransformation(
            min_cell_size=self.min_cell_size,
            max_cell_size=self.max_cell_size,
            enum_precision_parameter=self.enum_precision_parameter,
        )
        enumerated_structures = enum.apply_transformation(
            surface_structure_partial, return_ranked_list=max_structures
        )
        return enumerated_structures
