from surface_pd.core import Slab


def structure_filter(input_slabs,
                     direction,
                     criteria):
    """
    Filter out the slab models that does not satisfy the criteria.

    Args:
        input_slabs: Input slab model.
        direction: Lattice direction perpendicular to the surface.
        criteria: Lattice parameter perpendicular to the input(parent)
            slab model surface.

    Returns:
        List of structures.

    """
    filtered_structures = []
    for i, slab in enumerate(input_slabs):
        lattice = slab['structure'].lattice.abc
        ############################################
        # structure['structure'].remove_site_property(
        #     'selective_dynamics')
        # structure['structure'].to(
        #     fmt='poscar', filename='{}.vasp'.format(k))
        ############################################
        # Strict criteria -- keeps slabs with c lattice
        # as parent slab models
        if any((x - 0.0001) <= criteria
               <= (x + 0.0001) for x in lattice):
            slab['structure'] = Slab.from_sites(
                slab['structure']).check_rotate(criteria)
            filtered_structures.append(slab['structure'])
    return filtered_structures


def selective_dynamics_completion(structure: Slab,
                                  direction: int,
                                  dummy_species: list,
                                  center_bottom: float,
                                  center_top: float,
                                  tolerance: float):
    """
    Complete  the selective dynamics properties after the enumeration based
    on the rules of species and locations.

    Args:
        structure: Enumerated slab model.
        direction: Same as before.
        dummy_species: Same as before.
        center_bottom: Same as before.
        center_top: Same as before.
        tolerance: Same as before.

    Returns:
        Slab model with all sites have selective dynamics.
    """
    # print(structure)
    # print(dummy_species)
    for t in structure:
        for ds in dummy_species:
            if ds in t:
                t.properties = {
                    'selective_dynamics': [True, True, True]}
        if t.properties['selective_dynamics'] is None:
            if (center_bottom - tolerance <=
                    t.frac_coords[direction] <=
                    center_top + tolerance):
                t.properties = {
                    'selective_dynamics': [False, False, False]}
            else:
                t.properties = {
                    'selective_dynamics': [True, True, True]}
    return structure
