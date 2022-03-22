"""
Important functions are here.

"""

import copy
import collections

import numpy as np
from pymatgen.core import Structure
from pymatgen.core.surface import get_slab_regions
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.advanced_transformations \
    import EnumerateStructureTransformation
from pymatgen.transformations.standard_transformations \
    import SubstitutionTransformation

__author__ = "Xinhao Li"
__email__ = "xl2778@columbia.edu"
__date__ = "2022-03-15"

# %% Wrap back and equal lattices check
EPS = np.finfo(np.double).eps


class IncompatibleLatticeError(Exception):
    pass


class NoInversionSymmetryError(Exception):
    pass


def wrap_pbc(coo, slab_direction=2, tolerance=0.01):
    """
    Wrap fractional coordinates back into the unit cell.

    Arguments:
      slab_direction:
      coo (ndarray): 3-vector of fractional coordinates
      tolerance (float): tolerance in multiples of the lattice vector applied
        when wrapping in direction of the slab, so that atoms at the bottom
        of the slab are not wrapped to the top of the unit cell
    Returns:
      0.0 <= coo[i] < 1.0  for the directions within the slab planes and
      -tolerance <= coo[i] < (1.0 - tolerance)  for the slab direction

    """
    bounds = np.array([0.0, 0.0, 0.0])
    bounds[slab_direction] = tolerance
    for i in range(3):
        while coo[i] < -bounds[i]:
            coo[i] += 1.0
        while coo[i] >= (1.0 - bounds[i] - EPS):
            coo[i] -= 1.0
    return coo


def equal_lattices(lat1, lat2, dtol=0.001, atol=0.01):
    lencheck = not np.any(
        np.abs(np.array(lat1.lengths) - np.array(lat2.lengths)) > dtol)
    angcheck = not np.any(
        np.abs(np.array(lat1.angles) - np.array(lat2.angles)) > atol)
    return lencheck and angcheck


def group_atoms_by_layer(o_layer, diff=0.01, max_diff=0.03):
    """
    This function is used to group misclassified atoms into the right
    layers. For example, c_atom1 = 0.01, c_atoms2 = 0.02, they should be
    classified in one layer. But actually they are not. So this function
    here will search the difference between two closest atoms, if the
    difference is smaller than diff, they will be regrouped in the same
    layer.

    Args:
        o_layer: dictionary which contains different c fractional
          coordinates as keys and number of atoms as values.
        diff: Accepted c fractional coordinates difference between two
          atoms.
        max_diff: Maximum accepted c fractional coordinates difference
          between two atoms.

    Returns: dictionary which contains different c fractional
      coordinates as keys and number of atoms as values.

    """
    height_sorted = sorted(o_layer.keys(), reverse=True)
    res = dict()
    excluded_heights = set()
    for height in height_sorted:
        if height in excluded_heights:
            continue
        res[height] = o_layer[height]
        for i in range(1, int(max_diff / diff) + 1):
            curr_height = height - i * diff
            if curr_height not in o_layer:
                break
            res[height] += o_layer[curr_height]
            excluded_heights.add(curr_height)
    return res


def layer_classification(input_structure):
    """
    This function is used to classify different layers in the slab model. For
    the (001) surface, it will be classified as Li-plane, TM-plane,
    and O-plane. For the (104) surface, only one Li-TM-O plane will exist.
    Args:
        input_structure: slab models with "selective_dynamics" labeled for
        all sites

    Returns:
        c-fractional coordinates of the upper and lower limit of the
        central fixed region
        c-fractional coordinates of the surface Li/TM and O atoms (for the
        (104) surface, these two values can be the same)
        A dictionary which c-fractional coordinates of Li atoms as key and
        number of Li atoms at this layer as value
        A list which has the c-fractional coordinates of each layer

    """
    # Get the name of the TM species
    structure = input_structure

    # Initialize three dictionaries ot store layer information of Li, TM,
    # and O atoms
    Li_layers = collections.defaultdict(int)
    TM_layers = collections.defaultdict(int)
    O_layers = collections.defaultdict(int)

    for s in structure:
        if 'Li' in s:
            # Since the c fractional coordinates of atoms even in the
            # exactly same layers are not identical, therefore, all c
            # fractional coordinates will be rounded to 2 decimal. Any
            # misclassified atoms will be regrouped by
            # "group_atoms_by_layer" function.
            Li_layers[round(s.frac_coords[2], 2)] += 1
        elif 'O' in s:
            O_layers[round(s.frac_coords[2], 2)] += 1
        else:
            TM_layers[round(s.frac_coords[2], 2)] += 1

    # Further group atoms by layer
    Li_layers = group_atoms_by_layer(Li_layers)
    TM_layers = group_atoms_by_layer(TM_layers)
    O_layers = group_atoms_by_layer(O_layers)
    layers = sorted({**Li_layers, **TM_layers})
    O_layers_c_frac = sorted(O_layers)

    # Sorted out central layers based on the "selective_dynamics"
    center_sites = []
    for s in structure:
        if ('selective_dynamics' in s.properties
                and not any(s.properties['selective_dynamics'])):
            center_sites.append(copy.deepcopy(s))
    # Get central slab boundaries
    center = Structure.from_sites(center_sites)
    center_region = get_slab_regions(center)
    lower_limit, upper_limit = center_region[0]

    # Li-terminated slab model -- (001) surface
    if len(Li_layers) > len(TM_layers):
        surface_Li = sorted(Li_layers.items())[-1][0]
        surface_O = sorted(O_layers.items())[-1][0]
        return [lower_limit, upper_limit, surface_Li,
                surface_O, Li_layers, layers,
                O_layers_c_frac, TM_layers]
    # TM terminated slab model -- (001) surface
    elif len(Li_layers) < len(TM_layers):
        surface_TM = sorted(TM_layers.items())[-1][0]
        surface_O = sorted(O_layers.items())[-1][0]
        return [lower_limit, upper_limit, surface_TM,
                surface_O, Li_layers, layers,
                O_layers_c_frac, TM_layers]
    # Li-TM-O terminated slab model -- (104) surface
    else:
        surface_TM_Li = sorted(Li_layers.items())[-1][0]
        surface_O = sorted(O_layers.items())[-1][0]
        # print('The structure has polar surfaces.')
        return [lower_limit, upper_limit, surface_TM_Li,
                surface_O, Li_layers, layers,
                O_layers_c_frac, TM_layers]


def enum_with_composition(structure_model,
                          subs_li, li_composition,
                          subs_o, o_composition,
                          cell_size):
    """
    This function is used to enumerate the slab model which the target
    surface atoms have been substituted by the dummy species.
    Args:
        structure_model: slab model which included the dummy species on the 
        surface
        subs_o: substituted dummy species for surface O atoms
        o_composition: desired surface O composition after enumeration
        subs_li: substituted dummy species for surface O atoms
        li_composition: desired surface Li composition after enumeration
        cell_size: maximum/minimum cell size that want to be enumerated

    Returns: enumerated slab models

    """
    subs = SubstitutionTransformation(
        {subs_li: {subs_li: li_composition},
         subs_o: {subs_o: o_composition}})
    surface_structure_partial = subs.apply_transformation(
        structure_model)
    enum = EnumerateStructureTransformation(
        min_cell_size=1,
        max_cell_size=cell_size,
        enum_precision_parameter=0.00001
    )
    structures = enum.apply_transformation(
        surface_structure_partial, return_ranked_list=2000)
    return structures


def index_extraction(structure_model, tol=0.01):
    """
    This function is used to get the index of the first and second layers
    which containing Li atoms, surface oxygen atoms as well as the relaxed
    Li atoms in the slab model
    Args:
        structure_model: slab model which has "selective dynamics"
        labeled at the end.
        tol: Maximum spread (in fractional coordinates) that atoms in the
        same plane may exhibit (default: 0.01)

    Returns: indices of the first and second layers which containing Li
    atoms, surface oxygen atoms as well as the relaxed Li atoms in the slab
    model

    """
    [_, _, _, surface_O,
     li_layers, _,
     O_layers_c_frac, _] = layer_classification(structure_model)

    # Initialize two lists to store info of the first and second layers
    # which containing Li atoms (only be used for the (104) surface)
    li_frac_coords = list(li_layers.keys())
    first_layer, second_layer = [li_frac_coords[0], li_frac_coords[-1]], \
                                [li_frac_coords[1], li_frac_coords[-2]]

    # Get indices of the first and second Li layers as well as the surface
    # oxygen layers
    first_layer_index, second_layer_index, oxygen_index = [], [], []
    for index, site in enumerate(structure_model):
        # First layer Li atoms index extraction
        if first_layer[0] - tol <= site.frac_coords[2] <= \
                first_layer[0] + tol and "Li" in site:
            first_layer_index.append(index)
        if first_layer[1] - tol <= site.frac_coords[2] <= \
                first_layer[1] + tol and "Li" in site:
            first_layer_index.append(index)
        # Second layer Li atoms index extraction
        if second_layer[0] - tol <= site.frac_coords[2] <= \
                second_layer[0] + tol and "Li" in site:
            second_layer_index.append(index)
        if second_layer[1] - tol <= site.frac_coords[2] <= \
                second_layer[1] + tol and "Li" in site:
            second_layer_index.append(index)
        # Surface O atoms index extraction (both top and bottom surface O
        # atoms are all included)
        if (site.frac_coords[2] >= surface_O - tol and "O" in site) \
                or (
                site.frac_coords[2] <= O_layers_c_frac[0] + tol and "O" in site
        ):
            oxygen_index.append(index)

    relaxed_li_index = []
    for index, site in enumerate(structure_model):
        if site.properties["selective_dynamics"] == [True, True, True] \
                and "Li" in site:
            relaxed_li_index.append(index)

    return (first_layer_index, second_layer_index,
            oxygen_index, relaxed_li_index)


def remove_sites(structure_model,
                 index,
                 scaling=True,
                 scaling_matrix=None):
    """
    This function is used to remove atoms in the slab model first and
    then make a supercell based on the scaling matrix.
    Args:
        structure_model: parent slab model
        index: indices of the atoms that are going to be removed
        scaling: whether to create a supercell of the slab model with
        vacancies created (default: True)
        scaling_matrix: scaling matrix to be used to create the supercell

    Returns: scaled slab model with vacancies created on the surface

    """
    if scaling_matrix is None:
        scaling_matrix = [2, 1, 1]
    copy_structure = structure_model.copy()
    copy_structure.remove_sites(index)
    if scaling:
        copy_structure.make_supercell(scaling_matrix=scaling_matrix)
    return copy_structure


def symmetrize_top_base(target_slab, symprec=1e-4, direction=2, tol=0.01):
    """
    This function is used to symmetrize the enumerated slab models using
    top surface as base. The inversion symmetry center is determined based
    on the central fixed region only because the input slab of this function
    is not symmetry since one the one hand, the top surface of the slab is
    substituted by other dummy species, on the other hand, it has been
    enumerated with different Li and O compositions on the surface.
    Therefore using the whole slab to detect the inversion symmetry does not
    work.
    Args:
        target_slab: enumerated slab model with Li and oxygen vacancies on
        the top surface
        symprec: tolerance for the symmetry detection (default: 0.1).
        direction: lattice direction perpendicular to the surface (default: 2)
        tol:

    Returns: symmetrized slab models with top and bottom surfaces are
    equivalent

    """

    # Load slab which is going to be symmetrized
    slab_tgt = target_slab

    # Use center region to get the inversion symmetry center and the rotation
    # matrix
    lower_limit, upper_limit, _, _, _, _, _, _ = layer_classification(slab_tgt)

    # Generate the reference slab which is just the central fixed region
    slab_ref_site = []
    for s in slab_tgt:
        if lower_limit - tol < s.frac_coords[2] < upper_limit + tol:
            slab_ref_site.append(copy.deepcopy(s))
    slab_ref = Structure.from_sites(slab_ref_site)

    # Determine symmetry operations of the reference slab and make sure
    # the reference slab has an inversion center
    sga = SpacegroupAnalyzer(slab_ref, symprec=symprec)
    if not sga.is_laue():  # has Laue symmetry (centro-symmetry)
        raise NoInversionSymmetryError(
            "The target slab does not have inversion symmetry.  Try "
            "increasing the tolerance for symmetry detection with the "
            "'--symprec' flag.")

    # get the inversion operation and the inversion center (origin)
    ops = sga.get_symmetry_operations()
    inversion = ops[1]
    # print("trans_vec", inversion.translation_vector)
    assert (np.all(inversion.rotation_matrix == -np.identity(3)))
    origin = inversion.translation_vector / 2

    # wrap target slab models to unit cell
    for s in slab_tgt:
        s.frac_coords = wrap_pbc(s.frac_coords, slab_direction=direction)

    # Symmetrized slab models based on top only
    top_sites = []
    bottom_sites = []
    center_sites = []
    for s in [s for s in slab_tgt
              if s.frac_coords[direction] > origin[direction]]:
        if ('selective_dynamics' in s.properties
                and not any(s.properties['selective_dynamics'])):
            center_sites.append(copy.deepcopy(s))
        else:
            top_sites.append(copy.deepcopy(s))
    for s in [s for s in slab_tgt
              if s.frac_coords[direction] < origin[direction]]:
        if ('selective_dynamics' in s.properties
                and not any(s.properties['selective_dynamics'])):
            center_sites.append(copy.deepcopy(s))
        else:
            bottom_sites.append(copy.deepcopy(s))
    # Extend the center sites list in case of some sites have exactly the
    # same c fractional coordinates of the origin
    center_sites.extend([copy.deepcopy(s) for s in slab_tgt
                         if s.frac_coords[direction] == origin[direction]])

    # Initialize a new slab model which only has the center and top regions of
    # the initial slab model
    sites = center_sites[:] + top_sites[:]
    for s in top_sites:
        s2 = copy.deepcopy(s)
        s2.frac_coords = inversion.operate(s.frac_coords)
        # Add symmetrized top sites into the new slab as the bottom sites
        sites.append(s2)
    symmetrized_slab_top = Structure.from_sites(sites)

    for s in symmetrized_slab_top:
        s.frac_coords = wrap_pbc(s.frac_coords, slab_direction=direction)
    symmetrized_slab_top = symmetrized_slab_top.get_sorted_structure()

    return symmetrized_slab_top


def surface_substitute(target_slab, subs1, subs2,
                       direction=2, tol=0.02):
    """
    This function is used to substitute surface Li and O atoms with "dummy
    species" which will facilitate the enumeration code to detect target atoms.
    Args:
        target_slab: Parent slab model
        subs1: Substitution atom for O atom
        subs2: Substitution atom for Li atom
        direction: Define the direction of substitution (c direction is the
        general case) (default: 2)
        tol: Maximum spread (in fractional coordinates) that atoms in teh
        same plane may exhibit (default: 0.02)

    Returns: Substituted slab model

    """
    # Load slab
    slab_tgt = Structure.from_file(target_slab)

    # Define c-fractional coordinates of upper boundary of fixed region,
    # surface metal, and surface O atoms.
    [_, center_top, surface_metal,
     max_frac_O, _, _, _, _] = layer_classification(slab_tgt)

    # This distance is calculated by using the surface metal (can be Li or
    # TM  atoms) minus the top boundary of central fixed region.
    distance = (surface_metal - center_top)

    # Define criteria to determine surface O and Li atoms
    surface_O = []
    surface_Li = []
    for s in slab_tgt:
        if abs(max_frac_O - s.frac_coords[direction]) < tol and ('O' in s):
            surface_O.append(copy.deepcopy(s))
        if ((abs(surface_metal - s.frac_coords[direction]) < distance - tol)
                and ('Li' in s)):
            surface_Li.append(copy.deepcopy(s))

    # Extract indices for surface O and Li atoms in target slab model
    def indices(tgt_idx, ref_idx):
        idx = []
        for i, s in enumerate(tgt_idx):
            for j, t in enumerate(ref_idx):
                if t == s:
                    idx.append(i)
        return idx

    idx_surface_O = indices(slab_tgt, surface_O)
    idx_surface_Li = indices(slab_tgt, surface_Li)

    # Substitute surface Li and O atoms with dummy species
    slab_surface_substitute = slab_tgt.copy()
    for i in range(len(slab_surface_substitute)):
        if i in idx_surface_O:
            slab_surface_substitute.replace(
                i, subs1,
                properties={'selective_dynamics': [True, True, True]}
            )
        if i in idx_surface_Li:
            slab_surface_substitute.replace(
                i, subs2,
                properties={'selective_dynamics': [True, True, True]}
            )
    return slab_surface_substitute


def get_num_sites(lithiated_structure, slab_substituted,
                  cell_size,
                  Li_composition, O_composition):
    """
    Get number of sites in the slab model after enumeration
    Args:
        lithiated_structure: fully lithiated structure (input structure)
        slab_substituted: the slab model where the surface Li and O atoms
        are substituted
        cell_size: maximum cell size
        Li_composition: enumerated Li composition
        O_composition: enumerated O composition

    Returns: number of sites that should be after the enumeration

    """
    # Get number of lithium, TM, oxygen atoms in the fully lithiated slabs
    # after scaling
    enum_Li, enum_O = [slab_substituted.composition["Na+"] * cell_size * 2,
                       slab_substituted.composition["F2-"] * cell_size * 2]

    num_Li, num_O = [lithiated_structure.composition["Li"] * cell_size,
                     lithiated_structure.composition["O"] * cell_size]
    num_TM = lithiated_structure.num_sites * cell_size - num_Li - num_O

    rest_Li, rest_O = [num_Li - enum_Li, num_O - enum_O]
    curr_Li, curr_O = Li_composition * enum_Li, O_composition * enum_O
    curr_Li, curr_O = curr_Li + rest_Li, curr_O + rest_O
    curr_num_sites = curr_Li + num_TM + curr_O
    return curr_num_sites


def slab_size_check(
        refined_structure,
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
            .get_reduced_structure()
        return -1, refined_prim
    elif enumerated_num_sites < total_num_sites:
        if total_num_sites % enumerated_num_sites != 0:
            raise ValueError("Check primitive structure finder.")
        # else:
            # multiple = total_num_sites / enumerated_num_sites
            # a, b, c = refined_structure.lattice.abc
            # scaling_matrix = define_scaling_matrix(a, b, multiple)
            # refined_super = refined_structure.copy()
            # refined_super.make_supercell(scaling_matrix)
        return 1, refined_structure
    else:
        if max(refined_structure.lattice.abc) > input_c * 2 - 5:
            refined_prim = refined_structure.copy()
            refined_prim = refined_prim.get_primitive_structure() \
                .get_reduced_structure()
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
        print("Check symmetry!")
    if max(structure.lattice.abc) != structure.lattice.c:
        print("{}Li{}O -- structure_{}".format(Li_composition,
                                               O_composition,
                                               index))
        print("Check the orientation of the slab!")
    if any(origin) != 0.:
        print("{}Li{}O -- structure_{}".format(Li_composition,
                                               O_composition,
                                               index))
        print("Check inversion symmetry center shift!")


def define_scaling_matrix(a, b, multiple):
    """
    Generate the scaling matrix based on the differect values of the
    lattice parameters a and b
    Args:
        a: lattice parameter a
        b: lattice parameter a
        multiple: maximum target cell size

    Returns:
        scaling matrix that will be used to create the supercell
    """
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


def get_max_min_c_frac(structure):
    """
    Get the maximum and minimum values of the c fractional coordinates for
    all sites. This is used to make sure whether the slab region is in the
    middle or at the bottom of the unit cell
    Args:
        structure: enumerated, after refined slab model

    Returns:
        maximum and minimum values of the c fractional coordinates
    """
    min, max = 1, 0
    for site in structure:
        if site.frac_coords[2] > max:
            max = site.frac_coords[2]
        if site.frac_coords[2] < min:
            min = site.frac_coords[2]
    return min, max


def Li_TM_layers_finder(structure):
    """
    This function is used to determine the c fractional coordinates of
    central region of the slab models as well as the surface metal and O
    atoms.

    Args:
        input_structure: any slab models

    Returns: upper limit and lower limit of the central region of the
      slab, surface metal and O c fractional coordinates

    """
    Li_layers = collections.defaultdict(int)
    TM_layers = collections.defaultdict(int)
    O_layers = collections.defaultdict(int)
    for s in structure:
        if 'Li' in s:
            # Since the c fractional coordinates of atoms even in the
            # exactly same layers are not identical, therefore, all c
            # fractional coordinates will be rounded to 2 decimal. Any
            # misclassified atoms will be regrouped by
            # "group_atoms_by_layer" function.
            Li_layers[round(s.frac_coords[2], 2)] += 1
        elif 'O' in s:
            O_layers[round(s.frac_coords[2], 2)] += 1
        else:
            TM_layers[round(s.frac_coords[2], 2)] += 1

    Li_layers = group_atoms_by_layer(Li_layers)
    TM_layers = group_atoms_by_layer(TM_layers)
    O_layers = group_atoms_by_layer(O_layers)
    return Li_layers, TM_layers, O_layers


def boundary_define(parent_structure,
                    enumed_structure,
                    num_relaxed):
    """
    Get the boundary of the central fixed slab
    Args:
        parent_structure:
        enumed_structure:
        num_relaxed:

    Returns:

    """
    Li_layers, TM_layers, _ = Li_TM_layers_finder(parent_structure)
    enumed_Li_layers, enumed_TM_layers, _ = Li_TM_layers_finder(
        enumed_structure)
    layers = sorted({**enumed_Li_layers, **enumed_TM_layers})
    if len(Li_layers) == len(TM_layers):
        parent_num_layers = len(TM_layers)
        num_layers = len(enumed_TM_layers)
    else:
        parent_num_layers = len(Li_layers) + len(TM_layers)
        num_layers = len(enumed_Li_layers) + len(enumed_TM_layers)

    # Use number of layers in the parent slab model to define how many
    # layers should be fixed in the middel since this number should be
    # constant for all enumerated slabs.
    num_fixed = parent_num_layers - 2 * num_relaxed

    # Polar surface
    if len(Li_layers) != len(TM_layers):
        if (num_layers % 2) != 0:  # Odd number of layers
            # In the slab model, the number of layers starts from 1, so here
            # we need to add 1.
            center_layer = num_layers // 2 + 1
            # In dict, the number of layers starts from 0, therefore here we
            # need to minus 1.
            if num_fixed // 2 == 0:
                raise ValueError("The polar surface has to have the odd "
                                 "number of fixed layers.")
            else:
                lower_limit = layers[center_layer - 1 - (num_fixed // 2)]
                upper_limit = layers[center_layer - 1 + (num_fixed // 2)]
                return lower_limit, upper_limit
        else:
            raise ValueError("The polar surface has to have the odd number of "
                             "layers.")
    # Non-polar surface
    else:
        if (num_layers % 2) != 0:  # Odd number of layers
            center_layer = num_layers // 2 + 1
            if num_fixed // 2 == 0:
                raise ValueError("The number of fixed layers has to be odd "
                                 "for the non-polar surface with odd number "
                                 "of layers.")
            else:
                lower_limit = layers[center_layer - 1 - (num_fixed // 2)]
                upper_limit = layers[center_layer - 1 + (num_fixed // 2)]
                return lower_limit, upper_limit
        else:  # Even number of layers
            if (num_fixed % 2) != 0:
                raise ValueError("The number of fixed layers has to be even "
                                 "for the non-polar surface with even number "
                                 "of layers.")
            else:
                lower_limit = layers[num_relaxed - 1 + 1]
                upper_limit = layers[-num_relaxed - 1]
                return lower_limit, upper_limit


def add_selective_dynamics(parent_structure,
                           enumed_structure,
                           num_relaxed):
    """
    Add selective dynamics to the after refined slab model based on number
    of layers that will be relaxed on the surface
    Args:
        parent_structure:
        enumed_structure:
        num_relaxed:

    Returns:

    """
    lower_limit, upper_limit = boundary_define(parent_structure,
                                               enumed_structure,
                                               num_relaxed)
    for site in enumed_structure:
        if lower_limit - 0.01 <= site.frac_coords[2] <= upper_limit + 0.01:
            site.properties = {'selective_dynamics': [False, False, False]}
        else:
            site.properties = {'selective_dynamics': [True, True, True]}
    return enumed_structure


def temp_shift_isc_back(before_refine_structure,
                        after_refine_structure,
                        shift=True):
    """
    Shift / shift back the inversion symmetry center
    Args:
        before_refine_structure:
        after_refine_structure:
        shift: whether to turn on the shift

    Returns:

    """
    sga = SpacegroupAnalyzer(before_refine_structure,
                             symprec=1e-1)
    ops = sga.get_symmetry_operations()
    inversion = ops[1]
    assert (np.all(inversion.rotation_matrix == -np.identity(3)))
    origin = inversion.translation_vector / 2
    if shift:
        for site in after_refine_structure:
            site.frac_coords = site.frac_coords + origin
            wrap_pbc(site.frac_coords, slab_direction=2)
        return after_refine_structure
    else:
        for site in after_refine_structure:
            site.frac_coords = site.frac_coords - origin
            wrap_pbc(site.frac_coords, slab_direction=2)
        return after_refine_structure

