"""
Important functions are here.

"""

import copy

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
__date__ = "2021-12-03"

# %% Wrap back and equal lattices check
EPS = np.finfo(np.double).eps


class IncompatibleLatticeError(Exception):
    pass


class NoInversionSymmetryError(Exception):
    pass


def wrap_pbc(coo, slab_direction=2, tolerance=0.1):
    """
    Wrap fractional coordinates back into the unit cell.

    Arguments:
      slab_direction:
      coo (ndarray): 3-vector of fractional coordinates
      tolerance (float): tolerance in multiples of the lattice vector applied
        when wrapping in direction of the slab, so that atoms at the bottom
        of the slab are not wrapped to the top of the unit cell
      ToDo: tolerance could be determined automatically based on the center
            of the vacuum region
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
    symbol_sets = structure.symbol_set
    for species in symbol_sets:
        if species != "Li" and species != "O":
            TM_species = species

    # Initialize three dictionaries ot store layer information of Li, TM,
    # and O atoms
    Li_layers = {}
    TM_layers = {}
    O_layers = {}
    for s in structure:
        if 'Li' in s:
            # Since the c fractional coordinates of atoms even in the
            # exactly same layers are not identical, therefore, all c
            # fractional coordinates will be rounded to 2 decimal. Any
            # misclassified atoms will be regrouped by
            # "group_atoms_by_layer" function.
            if round(s.frac_coords[2], 2) not in Li_layers:
                Li_layers[round(s.frac_coords[2], 2)] = 1
            else:
                Li_layers[round(s.frac_coords[2], 2)] += 1
        elif str(TM_species) in s:
            if round(s.frac_coords[2], 2) not in TM_layers:
                TM_layers[round(s.frac_coords[2], 2)] = 1
            else:
                TM_layers[round(s.frac_coords[2], 2)] += 1
        else:
            if round(s.frac_coords[2], 2) not in O_layers:
                O_layers[round(s.frac_coords[2], 2)] = 1
            else:
                O_layers[round(s.frac_coords[2], 2)] += 1

    # Further group atoms by layer
    Li_layers = group_atoms_by_layer(Li_layers)
    TM_layers = group_atoms_by_layer(TM_layers)
    O_layers = group_atoms_by_layer(O_layers)
    layers = sorted({**Li_layers, **TM_layers})

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
                surface_O, Li_layers, layers]
    # TM terminated slab model -- (001) surface
    elif len(Li_layers) < len(TM_layers):
        surface_TM = sorted(TM_layers.items())[-1][0]
        surface_O = sorted(O_layers.items())[-1][0]
        return [lower_limit, upper_limit, surface_TM,
                surface_O, Li_layers, layers]
    # Li-TM-O terminated slab model -- (104) surface
    else:
        surface_TM_Li = sorted(Li_layers.items())[-1][0]
        surface_O = sorted(O_layers.items())[-1][0]
        # print('The structure has polar surfaces.')
        return [lower_limit, upper_limit, surface_TM_Li,
                surface_O, Li_layers, layers]


def enum_with_composition(structure_model,
                          subs_o, o_composition,
                          subs_li, li_composition,
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
        {subs_o: {subs_o: o_composition},
         subs_li: {subs_li: li_composition}})
    surface_structure_partial = subs.apply_transformation(
        structure_model)
    enum = EnumerateStructureTransformation(
        min_cell_size=cell_size,
        max_cell_size=cell_size)
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
        tol: Maximum spread (in fractional coordinates) that atoms in teh
        same plane may exhibit (default: 0.01)

    Returns: indices of the first and second layers which containing Li
    atoms, surface oxygen atoms as well as the relaxed Li atoms in the slab
    model

    """
    center_bottom, center_top, surface_metal, surface_O, li_layers, _ \
        = layer_classification(structure_model)

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
                or (site.frac_coords[2] <= 3 * tol and "O" in site):
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
    lower_limit, upper_limit, _, _, _, _ = layer_classification(slab_tgt)

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
    _, center_top, surface_metal, max_frac_O, _, _ = layer_classification(
        slab_tgt)

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


def shift_inversion_symmetry_center(file, symprec=1e-1):
    """
    This function is used to shift the inversion symmetry center to the (0,
    0, 0) point of the unit cell. This will facilitate VASP code to detect
    the symmetry of the slab model and can probably accelerate the convergence.
    Args:
        symprec (float): 
        file: enumerated slab model

    Returns: enumerated slab model with inversion symmetry center shifted

    """
    structure = file
    sga = SpacegroupAnalyzer(structure, symprec=symprec)
    if not sga.is_laue():
        raise ValueError("{} is not symmetric".format(str(file)))
    ops = sga.get_symmetry_operations()
    inversion = ops[1]
    assert (np.all(inversion.rotation_matrix == -np.identity(3)))
    origin = inversion.translation_vector / 2
    structure_copy = structure.copy()
    for site in structure_copy:
        site.frac_coords = site.frac_coords - origin
    sga = SpacegroupAnalyzer(structure_copy, symprec=symprec)
    ops = sga.get_symmetry_operations()
    inversion = ops[1]
    assert (np.all(inversion.rotation_matrix == -np.identity(3)))
    origin = inversion.translation_vector / 2
    if any(origin) != 0.:
        print("*****{} didn't shift properly*****!".format(str(file.formula)))
        print("*****The after shifted ISC is {}*****.".format(origin))
    return structure_copy