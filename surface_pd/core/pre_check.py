from pymatgen.core.surface import get_slab_regions

from surface_pd.core import Slab
from surface_pd.util import check_int


class PreCheck(object):
    """
    PreCheck class to check whether the input slab model is valid.

    Args:
        structure: Input slab model.

    """

    def __init__(self, structure: Slab):
        self.structure = structure

    def is_slab(self):
        """
        Check whether the input structure is a slab.

        """
        try:
            ranges = get_slab_regions(self.structure)
        except ValueError:
            return False
        else:
            if len(ranges) == 2:
                slab_width = ranges[0][1] - ranges[1][0]
            else:
                slab_width = ranges[0][1] - ranges[0][0]
            if slab_width < 0.9:
                return True
            else:
                return False

    def is_cuboid(self):
        """
        Check whether the input structure is a cuboid.

        """
        if max(self.structure.lattice.abc) != self.structure.lattice.c:
            return False
        else:
            return True

    def all_has_selective_dynamics(self):
        """
        Check if all sites have selective dynamics as the site properties.

        """
        for site in self.structure:
            try:
                if site.properties['selective_dynamics'] == \
                        ([False, False, False] or [True, True, True]):
                    return True
            except KeyError:
                return False
            else:
                if site.properties['selective_dynamics'] == [False]:
                    return False

    def has_inversion_symmetry(self, symprec=0.1):
        """
        Check if the input structure has inversion symmetry.

        Args:
            symprec: Symmetry detection tolerance. Defaults to 0.1. (
                fairly strict).

        """
        return Slab.from_sites(self.structure).is_symmetry(
            symprec=symprec,
            return_isc=False)

    def relax_both_surfaces(self):
        """
        Check if the input slab model has only one side of the surface
        will be relaxed.

        """
        # Wrap any out of unit cell atoms back
        self.structure.wrap_pbc()

        fixed_atoms_c_set = []
        for site in self.structure:
            if site.properties['selective_dynamics'] == [False, False, False]:
                fixed_atoms_c_set.append(
                    site.frac_coords[self.structure.direction])
        min_fixed_atom_c = min(fixed_atoms_c_set)

        # Check the slab region
        ranges = get_slab_regions(self.structure)
        if len(ranges) == 2:
            lower_boundary = min(ranges[0][1], ranges[1][0])
        else:
            lower_boundary = min(ranges[0])

        if min_fixed_atom_c - lower_boundary < self.structure.tolerance:
            print('********************************************************')
            print(' The slab model provided has whole bottom surface fixed. '
                  '\n '
                  ' Therefore, no symmetrization is needed.')
            print('********************************************************')
            return False
        else:
            return True

    def has_validate_composition(self,
                                 replace,
                                 max_cell_size):
        """
        Check whether the user defined composition can be enumerated.

        Args:
            replace: Species and occupancy dictionaries containing the species
                mapping in string-string pairs. E.g. {'Li': {'Li': 0.5}},
                stored in the list.
            max_cell_size: Maximum number of supercells of the input slab.

        """
        is_int = []

        num_atoms = {
            key: len(value)
            for key, value in self.structure.index_extraction()[2].items()
        }
        for comp in replace:
            for s, n_s in num_atoms.items():
                if self.structure.symmetric:
                    f = comp[s][s] * n_s * max_cell_size / 2
                    f = check_int(f)
                else:
                    f = comp[s][s] * n_s * max_cell_size
                    f = check_int(f)
                if f.is_integer():
                    is_int.append(True)
                else:
                    is_int.append(False)
                    break
        if all(is_int):
            return True
        else:
            return False
