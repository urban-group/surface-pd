import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from surface_pd.error import (PrimitiveStructureFinderError,
                              NoInversionSymmetryError,
                              SlabOrientationError,
                              NonCentralInversionSymmetryError)


class PostCheck(object):
    """
    PostCheck class to check whether the output slab model is valid.

    Args:
        refined_structure: After refined structure, to be checked

    """

    def __init__(self, refined_structure: Structure):
        self.refined_structure = refined_structure

    def slab_size_check(self,
                        total_num_sites,
                        enumerated_num_sites,
                        input_c):
        """
        Check whether the after refined structure has the right geometry.

        Args:
            total_num_sites: Number of sites that should be after enumeration
            enumerated_num_sites: Actual number of sites for the enumerated
                slab model
            input_c: c lattice parameter of the input(parent) slab model

        Returns:
            indicator and after-operated structure

        """
        if enumerated_num_sites > total_num_sites:
            refined_prim = self.refined_structure.copy()
            refined_prim = refined_prim.get_primitive_structure() \
                .get_reduced_structure().get_sorted_structure()
            return -1, refined_prim
        elif enumerated_num_sites < total_num_sites:
            if total_num_sites % enumerated_num_sites != 0:
                raise PrimitiveStructureFinderError
            # else:
            # multiple = total_num_sites / enumerated_num_sites
            # a, b, c = refined_structure.lattice.abc
            # scaling_matrix = define_scaling_matrix(a, b, multiple)
            # refined_super = refined_structure.copy()
            # refined_super.make_supercell(scaling_matrix)
            else:
                return 1, self.refined_structure
        else:
            if max(self.refined_structure.lattice.abc) > input_c * 2 - 5:
                refined_prim = self.refined_structure.copy()
                refined_prim = refined_prim.get_primitive_structure() \
                    .get_reduced_structure().get_sorted_structure()
                return -1, refined_prim
            if max(self.refined_structure.lattice.abc) != \
                    self.refined_structure.lattice.c:
                if (max(self.refined_structure.lattice.abc) ==
                        self.refined_structure.lattice.a):
                    refined_rotated = self.refined_structure.copy()
                    refined_rotated.make_supercell(
                        [[0, 0, 1],
                         [0, 1, 0],
                         [1, 0, 0]]
                    )
                elif (max(self.refined_structure.lattice.abc) ==
                      self.refined_structure.lattice.b):
                    refined_rotated = self.refined_structure.copy()
                    refined_rotated.make_supercell(
                        [[0, 1, 0],
                         [1, 0, 0],
                         [0, 0, 1]]
                    )
                return 2, refined_rotated
            else:
                return 0, self.refined_structure

    def final_check(self,
                    Li_composition,
                    O_composition,
                    index,
                    symprec=1e-5):
        """
        Perform final check to see whether the generated slab models are
        correct. \n
        1. Has the correct geometry. \n
        2. Has the inversion symmetry center even with a quite high symmetry
        detection parameter. \n
        3. Has the inversion symmetry center shifted to the origin (0,
        0, 0).

        Args:
            Li_composition: Surface Li composition of the slab model.
            O_composition: Surface O composition of the slab model.
            index: Index of the slab model in the dataset.
            symprec: Symmetry detection tolerance. Defaults to 1e-5 (tight
                tolerance is used here).

        Returns:
            If everything is fine, it will just pass. If something
            unexpected happens, the specific error report will generate.

        """

        sga = SpacegroupAnalyzer(self.refined_structure, symprec=symprec)
        ops = sga.get_symmetry_operations()
        inversion = ops[1]
        assert (np.all(inversion.rotation_matrix == -np.identity(3)))
        origin = inversion.translation_vector / 2
        if not sga.is_laue():
            print("{}Li{}O -- structure_{}".format(Li_composition,
                                                   O_composition,
                                                   index))
            raise NoInversionSymmetryError
        if max(self.refined_structure.lattice.abc) != \
                self.refined_structure.lattice.c:
            print("{}Li{}O -- structure_{}".format(Li_composition,
                                                   O_composition,
                                                   index))
            raise SlabOrientationError
        if any(origin) != 0.:
            print("{}Li{}O -- structure_{}".format(Li_composition,
                                                   O_composition,
                                                   index))
            raise NonCentralInversionSymmetryError
