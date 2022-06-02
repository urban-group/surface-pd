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
    """
    def __init__(self,
                 subs_species: list=None,
                 subs_composition: list=None,
                 min_cell_size: float=1,
                 max_cell_size: int=None,
                 enum_precision_parameter: float=0.00001):
        """

        Args:
            subs_species: to-be-enumerated species
            subs_composition: target composition
            min_cell_size:
            max_cell_size:
            enum_precision_parameter:
        """
        self.subs_species = subs_species
        self.subs_composition = subs_composition
        self.min_cell_size = min_cell_size
        self.cell_zie = max_cell_size
        self.enum_precision_parameter = enum_precision_parameter

    def apply_enumeration(self, structure: Structure):
        """

        Args:
            structure: parent slab model

        Returns:

        """
        num_subs = len(self.subs_species)

        subs_dict = collections.defaultdict(float)
        for i in range(num_subs):
            temp_dict = collections.defaultdict(float)
            temp_dict[self.subs_species[i]] = self.subs_composition[i]
            subs_dict[self.subs_species[i]] = temp_dict

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
