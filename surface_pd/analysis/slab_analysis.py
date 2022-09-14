from surface_pd.core import Slab


def structure_filter(input_slabs,
                     direction,
                     criteria):
    """
    Filter out the slab models that doese not satisfy the criteria.

    Args:
        input_slabs: Input slab model.
        direction: Lattice direction perpendicular to the surface.
        criteria: Lattice parameter perpendicular to the input(parent)
            slab model surface.

    Returns:
        List of structures.

    """
    filtered_structures = []
    lattice_sets = []
    for i, slab in enumerate(input_slabs):
        lattice = slab['structure'].lattice.abc
        lattice_sets.append(lattice)
        ############################################
        # structure['structure'].remove_site_property(
        #     'selective_dynamics')
        # structure['structure'].to(
        #     fmt='poscar', filename='{}.vasp'.format(k))
        ############################################
        # Strict criteria -- keeps slabs with c lattice
        # as parent slab models
        if (round(lattice[direction], 3) - 0.2) <= criteria <= \
                (round(lattice[direction], 3) + 0.2):
            filtered_structures.append(input_slabs[i]['structure'])

    # In case of that using strict criteria will remove all
    # structures, apply modest criteria to complete the
    # dataset.
    if len(filtered_structures) == 0:
        # print("**The criteria used is too strict! Changing to "
        #       "the modest one.**")
        for i, slab in enumerate(input_slabs):
            a, b, c_new = lattice_sets[i]
            # Keep the c direction of the slab models is
            # perpendicular to x-y plane but the c lattice
            # parameter can be modified for a little.
            if (a and b) < c_new <= criteria * 1.5:
                filtered_structures.append(input_slabs[i]['structure'])

    # In case of the enumerated structures are rotated:
    if len(filtered_structures) == 0:
        for i, slab in enumerate(input_slabs):
            slab['structure'] = Slab.from_sites(
                slab['structure']).check_rotate()
            a, b, c_new = lattice_sets[i]
            # Strict criteria -- keeps slabs with c lattice
            # as parent slab models
            if (round(c_new, 3) - 0.2) <= criteria <= (round(c_new, 3) + 0.2):
                filtered_structures.append(input_slabs[i]['structure'])

        # In case of that using strict criteria will remove all
        # structures, apply modest criteria to complete the
        # dataset.
    if len(filtered_structures) == 0:
        # print("**The criteria used is too strict! Changing to "
        #       "the modest one.**")
        for i, slab in enumerate(input_slabs):
            slab['structure'] = Slab.from_sites(
                slab['structure']).check_rotate()
            a, b, c_new = lattice_sets[i]
            # Keep the c direction of the slab models is
            # perpendicular to x-y plane but the c lattice
            # parameter can be modified for a little.
            if (a and b) < c_new <= criteria * 1.5:
                filtered_structures.append(input_slabs[i]['structure'])
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
