import copy

import numpy as np

from pymatgen.core.composition import Element

from surface_pd.core import Slab
from surface_pd.error import (PrimitiveStructureFinderError,
                              NoInversionSymmetryError,
                              SlabOrientationError,
                              NonCentralInversionSymmetryError)


class PostCheck(object):
    """
    PostCheck class to check whether the output slab model is valid.

    Args:
        refined_structure: After refined structure, to be checked.

    """

    def __init__(self, refined_structure: Slab):
        self.refined_structure = refined_structure

    def slab_size_check(self,
                        total_num_sites,
                        enumerated_num_sites,
                        criteria):
        """
        Check whether the after refined structure has the right geometry.

        Args:
            total_num_sites: Number of sites that should be after enumeration.
            enumerated_num_sites: Actual number of sites for the enumerated
                slab model.
            criteria: Lattice parameter perpendicular to the input(parent)
                slab model surface.

        Returns:
            Indicator and after-operated structure.

        """
        if enumerated_num_sites > total_num_sites:
            refined_prim = self.refined_structure.copy()
            refined_prim = refined_prim.get_primitive_structure() \
                .get_reduced_structure().get_sorted_structure()
            return -1, refined_prim
        elif enumerated_num_sites < total_num_sites:
            if total_num_sites % enumerated_num_sites != 0:
                raise PrimitiveStructureFinderError
            else:
                return 1, self.refined_structure
        else:
            # Handle the case that the slab is repeated in the space
            if max(self.refined_structure.lattice.abc) > criteria * 2 - 5:
                refined_prim = copy.deepcopy(self.refined_structure)
                refined_prim = refined_prim.get_primitive_structure() \
                    .get_reduced_structure().get_sorted_structure()
                return -1, refined_prim
            # Handle the rotated case
            if max(self.refined_structure.lattice.abc) != \
                    self.refined_structure.lattice.c:
                if (max(self.refined_structure.lattice.abc) ==
                        self.refined_structure.lattice.a):
                    refined_rotated = copy.deepcopy(self.refined_structure)
                    refined_rotated.make_supercell(
                        [[0, 0, 1],
                         [0, 1, 0],
                         [1, 0, 0]]
                    )
                elif (max(self.refined_structure.lattice.abc) ==
                      self.refined_structure.lattice.b):
                    refined_rotated = copy.deepcopy(self.refined_structure)
                    refined_rotated.make_supercell(
                        [[1, 0, 0],
                         [0, 0, 1],
                         [0, 1, 0]]
                    )
                return 2, refined_rotated
            else:
                return 0, self.refined_structure

    def final_check(self,
                    species,
                    composition_list,
                    index,
                    keep_symmetric,
                    criteria,
                    symprec=1e-5):
        """
        Perform final check to see whether the enumerated slab models are
        correct. \n
        1. Has the correct geometry. \n
        2. Has the inversion symmetry center even with a quite high symmetry
        detection parameter. \n
        3. Has the inversion symmetry center shifted to the origin (0,
        0, 0).

        Args:
            species: Target species that will be enumerated.
            composition_list: Compositions of enumerated species.
            index: Unique index of this slab model.
            keep_symmetric: Whether the symmetry of the slab model will be
                kept after the enumeration.
            criteria:
            symprec: Tolerance for symmetry finding. Defaults to 1e-5.

        Returns:
            If everything is fine, it will just pass. If something
            unexpected happens, the specific error report will generate.

        """

        if keep_symmetric:
            try:
                symmetric, origin, _ = Slab.from_sites(
                    self.refined_structure).is_symmetry(
                    symprec=symprec,
                    return_isc=True)
            except TypeError:
                symmetric = Slab.from_sites(
                    self.refined_structure).is_symmetry(
                    symprec, return_isc=False)
            if not symmetric:
                print("{}{} -- structure_{}".format(species,
                                                    composition_list,
                                                    index))
                self.refined_structure.to(
                    fmt='poscar',
                    filename='debug-structure-{}.vasp'.format(index))
                raise NoInversionSymmetryError
            else:
                if any(origin) != 0.:
                    print("{}{} -- structure_{}".format(species,
                                                        composition_list,
                                                        index))
                    self.refined_structure.to(
                        fmt='poscar',
                        filename='debug-structure-{}.vasp'.format(index))
                    raise NonCentralInversionSymmetryError
        else:
            pass
        if (round(criteria, 2) - 0.02) >=\
                round(self.refined_structure.lattice.abc[
                      self.refined_structure.direction], 2) or \
                round(criteria, 2) + 0.02 <= \
                round(self.refined_structure.lattice.abc[
                          self.refined_structure.direction], 3):
            print("{}{} -- structure_{}".format(species,
                                                composition_list,
                                                index))
            self.refined_structure.to(
                fmt='poscar',
                filename='debug-structure-{}.vasp'.format(index))
            raise SlabOrientationError
