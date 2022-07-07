import numpy as np

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from surface_pd.core import Slab
from surface_pd.error import TooLargeSlabError


def define_scaling_matrix(a, b, multiple):
    """
    Generate the scaling matrix based on the different values of the
    lattice parameters a and b.

    Args:
        a: Lattice parameter a
        b: Lattice parameter a
        multiple: Maximum target cell size

    Returns:
        Scaling matrix that will be used to create the supercell.
    """
    if multiple > 4:
        raise TooLargeSlabError
    if multiple == 4:
        scaling_matrix = [2, 2, 1]
    else:
        if round(a / b, 0) == 2:
            scaling_matrix = [1, multiple, 1]
        elif round(b / a, 0) == 2:
            scaling_matrix = [multiple, 1, 1]
        else:
            scaling_matrix = [multiple, 1, 1]
    return scaling_matrix


def temp_shift_isc_back(before_refined_structure,
                        after_refined_structure,
                        shift=True):
    """
    Shift / shift back the inversion symmetry center.

    Args:
        before_refined_structure: Before refined structure
        after_refined_structure: After refined structure
        shift: Whether to shift the inversion symmetry center to the origin.

    Returns:
        Slab model with inversion symmetry center shifted to the origin /
        back.

    """
    sga = SpacegroupAnalyzer(before_refined_structure,
                             symprec=1e-1)
    ops = sga.get_symmetry_operations()
    inversion = ops[1]
    assert (np.all(inversion.rotation_matrix == -np.identity(3)))
    origin = inversion.translation_vector / 2
    if shift:
        for site in after_refined_structure:
            site.frac_coords = site.frac_coords + origin
            # wrap_pbc(site.frac_coords, slab_direction=2)
        after_refined_structure = Slab.from_sites(after_refined_structure)
        after_refined_structure.wrap_pbc()
        return after_refined_structure
    else:
        for site in after_refined_structure:
            site.frac_coords = site.frac_coords - origin
            # wrap_pbc(site.frac_coords, slab_direction=2)
        after_refined_structure = Slab.from_sites(after_refined_structure)
        after_refined_structure.wrap_pbc()
        return after_refined_structure
