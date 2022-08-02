from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.surface import get_slab_regions
from pymatgen.core.composition import Composition, Element

from surface_pd.error import *
from surface_pd.core import Slab


class PreCheck(object):
    """
    PreCheck class to check whether the input slab model is valid.

    Args:
        structure: Input slab model.

    """

    def __init__(self, structure: Slab):
        self.structure = structure

    def compound_check(self):
        """
        This function is used to check whether it is a Li-TM-O compound,
        if it is not Li-TM-O compound, the name of other elements will be
        internally changed into pseudo "Li-TM-O" compound.

        Returns:
            If it is a Li-TM-O system, pseudo Li element,
            pseudo TM element, and pseudo O element.
        """
        pseudo_Li, pseudo_TM, pseudo_O = '', '', ''

        num_pseudo_O = 0
        symbol_set = self.structure.symbol_set
        for element in symbol_set:
            # The slab model should always have more oxygen atoms than others.
            if (self.structure.composition[element]) > num_pseudo_O:
                num_pseudo_O = self.structure.composition[element]
                pseudo_O = element
            # The pseudo Li element should have a positive charge "+1".
            if list(Element(element).common_oxidation_states) == [1]:
                pseudo_Li = element

        symbol_set = list(symbol_set)
        for i in (pseudo_Li, pseudo_O):
            symbol_set.remove(i)
        pseudo_TM = symbol_set.copy()

        if pseudo_Li == 'Li' and pseudo_O == 'O':
            return True, pseudo_Li, pseudo_TM, pseudo_O
        else:
            return False, pseudo_Li, pseudo_TM, pseudo_O

    def pseudo_compound_generator(self):
        """
        This function is used to replace pseudo Li and O element with "Li"
        and "O" element.

        Returns:
            Li-TM-O system structure

        """
        _, pseudo_Li, pseudo_TM, pseudo_O = self.compound_check()
        self.structure.replace_species({Element(pseudo_Li): Element("Li"),
                                        Element(pseudo_O): Element("O")})
        return self.structure

    def is_not_slab(self):
        """
        Check whether the input structure is a slab.

        """
        try:
            ranges = get_slab_regions(self.structure)
        except ValueError:
            return True
        else:
            if len(ranges) == 2:
                slab_width = ranges[0][1] - ranges[1][0]
            else:
                slab_width = ranges[0][1] - ranges[0][0]
            if slab_width < 0.9:
                return False
            else:
                return True

    def is_not_cuboid(self):
        """
        Check whether the input structure is a cuboid.

        """
        if max(self.structure.lattice.abc) != self.structure.lattice.c:
            return True
        else:
            return False

    def not_all_has_selective_dynamics(self):
        """
        Check if all sites have selective dynamics as the site properties.

        """
        for site in self.structure:
            try:
                if site.properties['selective_dynamics'] == \
                        ([False, False, False] or [True, True, True]):
                    return False
            except KeyError:
                return True
            else:
                if site.properties['selective_dynamics'] == [False]:
                    return True

    def has_no_inversion_symmetry(self, symprec=0.1):
        """
        Check if the input structure has inversion symmetry.

        Args:
            symprec: Symmetry detection tolerance. Defaults to 0.1. (
                relatively loose).

        """
        # sga = SpacegroupAnalyzer(self.structure, symprec=symprec)
        # return not sga.is_laue()
        return not Slab.from_sites(self.structure). \
            is_symmetry(symprec=symprec,
                        return_isc=False)

    @staticmethod
    def has_validate_composition(num_lithium_like_atoms: int,
                                 lithium_like_composition: list,
                                 num_oxygen_like_atoms: int,
                                 oxygen_like_composition: list):
        pass_1, pass_2 = [], []
        for composition in lithium_like_composition:
            f = composition * num_lithium_like_atoms
            if f.is_integer():
                pass_1.append(True)
            else:
                pass_1.append(False)
                break
        for composition in oxygen_like_composition:
            f = composition * num_oxygen_like_atoms
            if f.is_integer():
                pass_2.append(True)
            else:
                pass_2.append(False)
                break

        if all(pass_1) and all(pass_2):
            return True
        else:
            return False
