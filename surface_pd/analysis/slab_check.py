import numpy as np

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from surface_pd.error import (PrimitiveStructureFinderError,
                              NoInversionSymmetryError,
                              SlabOrientationError,
                              NonCentralInversionSymmetryError)


def equal_lattices(lat1, lat2, dtol=0.001, atol=0.01):
    # lattice_check
    lencheck = not np.any(
        np.abs(np.array(lat1.lengths) - np.array(lat2.lengths)) > dtol)
    angcheck = not np.any(
        np.abs(np.array(lat1.angles) - np.array(lat2.angles)) > atol)
    return lencheck and angcheck


def slab_size_check(refined_structure,
                    total_num_sites,
                    enumerated_num_sites,
                    input_c):
    """
    Check whether the after refined structure has the right geometry that we
    expect
    Args:
        refined_structure: after refined slab model
        total_num_sites: number of sites taht should be after enumeration
        enumerated_num_sites: actual number of sites for the enumerated slab
        model
        input_c: c lattice parameter of the input(parent) slab model

    Returns:

    """
    if enumerated_num_sites > total_num_sites:
        refined_prim = refined_structure.copy()
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
            return 1, refined_structure
    else:
        if max(refined_structure.lattice.abc) > input_c * 2 - 5:
            refined_prim = refined_structure.copy()
            refined_prim = refined_prim.get_primitive_structure() \
                .get_reduced_structure().get_sorted_structure()
            return -1, refined_prim
        if (max(refined_structure.lattice.abc) != refined_structure.lattice.c):
            if (max(refined_structure.lattice.abc) ==
                    refined_structure.lattice.a):
                refined_rotated = refined_structure.copy()
                refined_rotated.make_supercell(
                    [[0, 0, 1],
                     [0, 1, 0],
                     [1, 0, 0]]
                )
            elif (max(refined_structure.lattice.abc) ==
                    refined_structure.lattice.b):
                refined_rotated = refined_structure.copy()
                refined_rotated.make_supercell(
                    [[0, 1, 0],
                     [1, 0, 0],
                     [0, 0, 1]]
                )
            return 2, refined_rotated
        else:
            return 0, refined_structure


def final_check(structure,
                input_c,
                Li_composition,
                O_composition,
                index,
                tol=0.5,
                symprec=1e-5):
    # slab_check
    """
    Perform final check to see whether the generated slab models are correct.
    1. has the correct geometry
    2. has the inversion symmetry center even with a quite high symmetry
    detection parameter
    3. has the inversion symmetry center shifted to the origin (0, 0, 0)
    Args:
        structure:
        input_c:
        Li_composition:
        O_composition:
        index:
        tol:
        symprec:

    Returns:

    """
    sga = SpacegroupAnalyzer(structure, symprec=symprec)
    ops = sga.get_symmetry_operations()
    inversion = ops[1]
    assert (np.all(inversion.rotation_matrix == -np.identity(3)))
    origin = inversion.translation_vector / 2
    if not sga.is_laue():
        print("{}Li{}O -- structure_{}".format(Li_composition,
                                               O_composition,
                                               index))
        raise NoInversionSymmetryError
    if max(structure.lattice.abc) != structure.lattice.c:
        print("{}Li{}O -- structure_{}".format(Li_composition,
                                               O_composition,
                                               index))
        raise SlabOrientationError
    if any(origin) != 0.:
        print("{}Li{}O -- structure_{}".format(Li_composition,
                                               O_composition,
                                               index))
        raise NonCentralInversionSymmetryError