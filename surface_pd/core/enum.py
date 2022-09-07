from surface_pd.core.slab import Slab

from pymatgen.transformations.standard_transformations \
    import SubstitutionTransformation
from pymatgen.transformations.advanced_transformations \
    import EnumerateStructureTransformation


class EnumWithComposition(object):
    """
    EnumWithComposition class to enumerate the parent slab model with
    defined composition.

    Args:
        subs_dict (dict): Species and occupancy dictionaries containing the
            species mapping in string-string pairs.
        min_cell_size (int): The minimum cell size. Must be an int.
            Defaults to 1.
        max_cell_size (int): The maximum cell size. Must be an int.
            Defaults to 1.
        enum_precision_parameter (float): Finite precision parameter for
            enumlib. Defaults to 1e-5.
    """

    def __init__(self,
                 subs_dict: dict,
                 min_cell_size: int = 1,
                 max_cell_size: int = 1,
                 enum_precision_parameter: float = 1e-5):
        self.subs_dict = subs_dict
        self.min_cell_size = min_cell_size
        self.max_cell_zie = max_cell_size
        self.enum_precision_parameter = enum_precision_parameter

    def apply_enumeration(self,
                          structure: Slab,
                          max_structures: int = 2000):
        """
        Apply enumeration to parent slab model.

        Args:
            structure: Parent slab model.
            max_structures: Number of structures to be returned at most for
                each composition. Defaults to 2000.

        Returns:
            A sequence of all enumerated structures.

        """

        subs = SubstitutionTransformation(self.subs_dict)
        surface_structure_partial = subs.apply_transformation(
            structure)
        enum = EnumerateStructureTransformation(
            min_cell_size=self.min_cell_size,
            max_cell_size=self.max_cell_zie,
            enum_precision_parameter=self.enum_precision_parameter
        )
        enumerated_structures = enum.apply_transformation(
            surface_structure_partial,
            return_ranked_list=max_structures)
        return enumerated_structures
