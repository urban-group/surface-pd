"""
Slab analysis functions for structure filtering and processing.

This module provides functions to filter enumerated slab structures based on
geometric criteria and complete selective dynamics properties for DFT
calculations.
"""

from surface_pd.core.enumeration_slab import EnumerationSlab


def structure_filter(input_slabs, direction, criteria):
    """
    Filter out the slab models that does not satisfy the criteria.

    Args:
        input_slabs: Input slab model.
        direction: Lattice direction perpendicular to the surface.
        criteria: Lattice parameter perpendicular to the input(parent)
            slab model surface.

    Returns
    -------
        List of structures.

    """
    filtered_structures = []
    for slab in input_slabs:
        lattice = slab["structure"].lattice.abc
        ############################################
        # structure['structure'].remove_site_property(
        #     'selective_dynamics')
        # structure['structure'].to(
        #     fmt='poscar', filename='{}.vasp'.format(k))
        ############################################
        # Strict criteria -- keeps slabs with c lattice
        # as parent slab models
        if any((x - 0.0001) <= criteria <= (x + 0.0001) for x in lattice):
            slab["structure"] = EnumerationSlab.from_sites(
                slab["structure"]
            )._check_rotate(criteria)
            filtered_structures.append(slab["structure"])
    return filtered_structures


def selective_dynamics_completion(
    structure: EnumerationSlab,
    direction: int,
    dummy_species: list,
    fixed_bottom_angstrom: float,
    fixed_top_angstrom: float,
):
    """
    Complete selective dynamics properties after enumeration.

    Args:
        structure: Enumerated slab model.
        direction: Same as before.
        dummy_species: Same as before.
        fixed_bottom_angstrom: Lower Cartesian bound of the fixed region.
        fixed_top_angstrom: Upper Cartesian bound of the fixed region.

    Returns
    -------
        Enumeration slab with selective dynamics on all sites.
    """
    # print(structure)
    # print(dummy_species)
    structure.direction = direction
    coordinates = structure._site_coordinates_angstrom()
    for index, t in enumerate(structure):
        for ds in dummy_species:
            if ds in t:
                t.properties = {"selective_dynamics": [True, True, True]}
        if t.properties["selective_dynamics"] is None:
            if (
                fixed_bottom_angstrom - 1e-8
                <= coordinates[index]
                <= fixed_top_angstrom + 1e-8
            ):
                t.properties = {"selective_dynamics": [False, False, False]}
            else:
                t.properties = {"selective_dynamics": [True, True, True]}
    return structure
