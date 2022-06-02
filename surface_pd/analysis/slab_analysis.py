"""

"""
from surface_pd.core import Slab
from surface_pd.error import PolarSurfaceError, NonPolarSurfaceError

def get_num_sites(lithiated_structure, slab_substituted,
                  Li_replacement, Li_composition,
                  O_replacement, O_composition,
                  cell_size):
    """
    Get the right number of sites in the slab model after enumeration
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
    enum_Li, enum_O = [
        slab_substituted.composition[Li_replacement] * cell_size * 2,
        slab_substituted.composition[O_replacement] * cell_size * 2
    ]

    num_Li, num_O = [lithiated_structure.composition["Li"] * cell_size,
                     lithiated_structure.composition["O"] * cell_size]
    num_TM = lithiated_structure.num_sites * cell_size - num_Li - num_O

    rest_Li, rest_O = [num_Li - enum_Li, num_O - enum_O]
    curr_Li, curr_O = Li_composition * enum_Li, O_composition * enum_O
    curr_Li, curr_O = curr_Li + rest_Li, curr_O + rest_O
    curr_num_sites = curr_Li + num_TM + curr_O
    return curr_num_sites


def boundary_define(parent_structure,
                    enumed_structure,
                    num_relaxed):
    """
    Get the boundary of the central fixed slab without the "selective
    dynamics" labeled
    Args:
        parent_structure:
        enumed_structure:
        num_relaxed:

    Returns:

    """
    parent_structure = Slab.from_sites(parent_structure)
    enumed_structure = Slab.from_sites(enumed_structure)
    Li_layers, TM_layers, _ = parent_structure.Li_TM_layers_finder()
    # Li_layers, TM_layers, _ = Li_TM_layers_finder(parent_structure)
    [enumed_Li_layers,
     enumed_TM_layers, _] = enumed_structure.Li_TM_layers_finder()
    # enumed_Li_layers, enumed_TM_layers, _ = Li_TM_layers_finder(
    #     enumed_structure)
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
            if num_fixed % 2 == 0:
                raise PolarSurfaceError(num_layers, num_fixed)
            else:
                lower_limit = layers[center_layer - 1 - (num_fixed // 2)]
                upper_limit = layers[center_layer - 1 + (num_fixed // 2)]
                return lower_limit, upper_limit
        else:
            raise PolarSurfaceError(num_layers, num_fixed)
    # Non-polar surface
    else:
        if (num_layers % 2) != 0:  # Odd number of layers
            center_layer = num_layers // 2 + 1
            if num_fixed % 2 == 0:
                raise NonPolarSurfaceError(num_layers, num_fixed)
            else:
                lower_limit = layers[center_layer - 1 - (num_fixed // 2)]
                upper_limit = layers[center_layer - 1 + (num_fixed // 2)]
                return lower_limit, upper_limit
        else:  # Even number of layers
            if (num_fixed % 2) != 0:
                raise NonPolarSurfaceError(num_layers, num_fixed)
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

