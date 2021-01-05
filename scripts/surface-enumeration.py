#!/usr/bin/env python

"""
TODO: Explain

"""

import argparse
import os

from pymatgen.transformations.standard_transformations \
    import SubstitutionTransformation
from pymatgen.transformations.advanced_transformations \
    import EnumerateStructureTransformation

from surface_pd.surface_enum import (surface_substitute,
                                     layer_classification,
                                     symmetrize_top_base)

__author__ = "Xinhao Li"
__email__ = "xl2778@columbia.edu"
__date__ = "2021-01-05"


def automate_surface(target_slab, target_cell_size=1, to_vasp=False):
    """
    The general purpose of this code is used to generate slab models
    with O and Li vacancies automatically.  The first step we are gonna
    do is to substitute top surface O and Li atoms with other
    elements. This will facilitate the Enumerate code to find the
    surface atoms more easily.  The second step is to enumerate the slab
    model with substituted surface atoms (only top surface atoms).  The
    third step is to do some post-processing to distinguish enumerated
    slab models with desired cell size.  The criteria used here includes
    cell size and length of c lattice parameter.  The final step is to
    symmetrize the enumerated slab models. At the beginning, the initial
    slab model is enumerated based on the top surface only. This is
    because that the slab models have to be symmetry. Enumeration
    process can also be done on both top and bottom surfaces, but
    finally top and bottom surfaces have to be symmetry.  Therefore, in
    order to reduce the redundant enumeration process on bottom surface,
    we will just enumerate the top surface first and then symmetrize it
    to the bottom, making a symmetry slab model.

    Args:
        target_slab: inpur slab model with no vacancies on the surface,
          desired cell size = 1x2x1 or 2x1x1 for now
        to_vasp: Whether generate the output structures

    Returns:

    """

    print("target_cell_size = {}".format(target_cell_size))

    # Load initial slab model with no vacancies on the surface
    slab_tgt = surface_substitute(target_slab, subs1='F', subs2='Na')

    # Extract c parameter and volume of the primitive unit cell (will be
    # used as criteria next)
    c = slab_tgt.lattice.abc[2]
    volume = slab_tgt.volume

    # Enumerate with maximum unit cell of 2
    # For 3 layers relaxed only (generates ~1000 structures)
    # composition_Li = [1, 0.833, 0.667, 0.5, 0.333, 0.167, 0]
    # For 2 layers
    composition_Li = [1, 0.75, 0.5, 0.25]
    composition_O = [1, 0.75, 0.5, 0.25]
    # composition = [1, 0.5, 0]
    # composition1 = [1]
    # composition2 = [0.5]
    num = 0

    for i in composition_O:
        for j in composition_Li:
            if (i == 1 and j == 1) or (i == 0 and j == 0):
                continue

            subs = SubstitutionTransformation({"F": {"F": i}, "Na": {"Na": j}})
            surface_structure_partial = subs.apply_transformation(slab_tgt)
            # FIX: min_cell_size cannot be set to anything greater than 1
            enum = EnumerateStructureTransformation(
                min_cell_size=target_cell_size, max_cell_size=target_cell_size)
            structures = enum.apply_transformation(
                surface_structure_partial, return_ranked_list=2000)

            # C-parameter and cell size check
            new_structures = []
            for k, s in enumerate(structures):
                # Keep volume constant and c lattice parameter unchanged
                if ((-1 < s['structure'].lattice.abc[2] - c < 1) and
                    (-1 < s['structure'].volume
                     - target_cell_size*volume < 1)):
                    new_structures.append(structures[k]['structure'])
            num += len(new_structures)

            print(f'The enumeration found {len(structures)} '
                  f'and {len(new_structures)} '
                  f'distinct structures for {j * 100}% Li and {i * 100}% O '
                  f'before and after filter respectively.')

            # Add selective dynamics for enumerated sites
            for structure in new_structures:
                for t in structure:
                    if 'Na' in t:
                        t.properties = {
                            'selective_dynamics': [True, True, True]}
                    if 'F' in t:
                        t.properties = {
                            'selective_dynamics': [True, True, True]}
                    if t.properties['selective_dynamics'] is None:
                        if layer_classification(slab_tgt)[0]-0.01 \
                           <= t.frac_coords[2] \
                           <= layer_classification(slab_tgt)[1]+0.01:
                            t.properties = {
                                'selective_dynamics': [False, False, False]}
                        else:
                            t.properties = {
                                'selective_dynamics': [True, True, True]}

            # Symmetrize structure models
            symmetrized_structure = []
            for s in new_structures:
                s.replace_species({'Na': 'Li', 'F': 'O'})
                symmetrized_structure.append(symmetrize_top_base(s))

            if to_vasp:
                # Generate structure models
                dirname = str(j) + 'Li' + str(i) + 'O'
                if not os.path.exists(dirname):
                    os.makedirs(dirname)
                for l, s in enumerate(symmetrized_structure):
                    s.to(fmt='poscar', filename=os.path.join(
                        dirname, "structure-{}.vasp".format(l)))

    print(f'{num} distinct structures are found totally.')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__+"\n{} {}".format(__date__, __author__),
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(
        "input_file",
        help="Path to an input file in VASP's POSCAR format.")

    parser.add_argument(
        "--generate-poscars", "-g",
        help="Generate POSCAR files of enumerated structures.",
        action="store_true")

    parser.add_argument(
        "--target-cell-size", "-s",
        help="Target cell size relative to input cell (default: 1).",
        type=int,
        default=1)

    args = parser.parse_args()

    automate_surface(args.input_file,
                     target_cell_size=args.target_cell_size,
                     to_vasp=args.generate_poscars)
