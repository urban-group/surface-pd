import collections

from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Species
from pymatgen.transformations.standard_transformations \
    import SubstitutionTransformation
from pymatgen.transformations.advanced_transformations \
    import EnumerateStructureTransformation


class EnumWithComposition(object):
    """
    EnumWithComposition class to enumerate the parent slab model with
    defined composition.

    Args:
        subs_species (list): to-be-enumerated species.
        subs_composition (list): target composition.
        min_cell_size (int): The minimum cell size. Must be an int.
            Defaults to 1.
        max_cell_size (int): The maximum cell size. Must be an int.
            Defaults to 1.
        enum_precision_parameter (float): Finite precision parameter for
            enumlib.
    """

    def __init__(self,
                 subs_species: list = None,
                 subs_composition: list = None,
                 min_cell_size: int = 1,
                 max_cell_size: int = 1,
                 enum_precision_parameter: float = 1e-5):
        self.subs_species = subs_species
        self.subs_composition = subs_composition
        self.min_cell_size = min_cell_size
        self.cell_zie = max_cell_size
        self.enum_precision_parameter = enum_precision_parameter

    def apply_enumeration(self, structure: Structure):
        """
        Apply enumeration to parent slab model.

        Args:
            structure (Structure): parent slab model.

        Returns:
            A sequence of all enumerated structures.
        """
        num_subs = len(self.subs_species)

        subs_dict = collections.defaultdict(float)
        for i in range(num_subs):
            temp_dict = collections.defaultdict(float)
            temp_dict[self.subs_species[i]] = self.subs_composition[i]
            subs_dict[self.subs_species[i]] = temp_dict

        # Questions 1: what the "subs_dict" dictionary above looks like?
        # Questions 2: When I looked at the self.surface_substitute function in Slab Class, I thought that Li and O atoms can be substituted by only one element, respectively.
        # In this EnumWithComposition Class, is there no limitations in the List length for "subs_species" and "subs_composition" arguments?
        # The List arguments with arbitrary length can be handled?
        subs = SubstitutionTransformation(subs_dict)
        surface_structure_partial = subs.apply_transformation(
            structure)
        
        enum = EnumerateStructureTransformation(
            min_cell_size=self.min_cell_size,
            max_cell_size=self.cell_zie,
            enum_precision_parameter=self.enum_precision_parameter
        )
        enumerated_structures = enum.apply_transformation(
            surface_structure_partial, return_ranked_list=2000)
        return enumerated_structures
