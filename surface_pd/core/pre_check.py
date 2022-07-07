from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.surface import get_slab_regions

from surface_pd.error import *


class PreCheck(object):
    """
    PreCheck class to check whether the input slab model is valid.

    Args:
        structure: Input slab model.

    """

    def __init__(self, structure: Structure):
        self.structure = structure

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
        sga = SpacegroupAnalyzer(self.structure, symprec=symprec)
        return not sga.is_laue()
