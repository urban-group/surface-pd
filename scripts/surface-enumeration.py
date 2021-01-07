#!/usr/bin/env python

"""
This code will generate slab models with O and Li vacancies on the surface
    automatically. The input slab structure should be as small as
    possible because it will increase the possibilities it can have.
    automatically.
"""

__author__ = "Xinhao Li"
__email__ = "xl2778@columbia.edu"
__date__ = "2021-01-05"

import argparse
import os

from pymatgen.transformations.standard_transformations \
    import SubstitutionTransformation
from pymatgen.transformations.advanced_transformations \
    import EnumerateStructureTransformation
from pymatgen.core.periodic_table import DummySpecies, Species

from surface_pd.surface_enum import (surface_substitute,
                                     layer_classification,
                                     symmetrize_top_base)

"""
Finished:
1. More arguments are added, such as the desired compositions.

Problems: 
1. When 0% Li and 75.0% O, 0% Li and 50.0% O and 0% Li and 25.0% 
    O, it can enumerate. But for 0% Li and 100% O, 100% Li and 0% O, 
    75.0% Li and 0% O, 50.0% Li and 0% O and 25.0% Li and 0% O, it can 
    not enumerate. 
2. Tried to replace S and Na with dummy species, but it does not work.
    The periodicity of atoms is not represented. For example, whatever
    the composition of oxygen is set on the surface, all the oxygen
    atoms on the surface will be deleted.
"""


def automate_surface(target_slab,
                     surface_oxygen_composition,
                     surface_lithium_composition,
                     target_cell_size=1,
                     to_vasp=False):
    """
    The first step here is to substitute top surface O and Li atoms with
    other elements. This will facilitate the Enumerate code to find the
    surface atoms more easily.  The second step is to enumerate the
    slab model with substituted surface atoms (only top surface atoms).
    The third step is to do some post-processing to distinguish enumerated
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
        target_slab:
        surface_oxygen_composition:
        surface_lithium_composition:
        target_cell_size:
        to_vasp:

    Returns: All enumerate structures in different folders.

    """

    print("target_cell_size = {}".format(target_cell_size))
    print("Composition of oxygen on the surface will be {}.".format(
        surface_oxygen_composition))
    print("Composition of lithium on the surface will be {}.".format(
        surface_lithium_composition))

    # Load initial slab model with no vacancies on the surface
    # O_replacement = DummySpecies(symbol='X', oxidation_state=-2)
    # Li_replacement = DummySpecies(symbol='Z', oxidation_state=+1)
    O_replacement = Species(symbol='S', oxidation_state=-2)
    Li_replacement = Species(symbol='Na', oxidation_state=+1)
    slab_tgt = surface_substitute(target_slab,
                                  subs1=O_replacement,
                                  subs2=Li_replacement)

    # Extract c parameter and volume of the primitive unit cell (will be
    # used as criteria next)
    c = slab_tgt.lattice.abc[2]
    volume = slab_tgt.volume

    # Enumerate with maximum unit cell of 2
    # For 3 layers relaxed only (generates ~1000 structures)
    # composition_Li = [1, 0.833, 0.667, 0.5, 0.333, 0.167, 0]
    # For 2 layers
    composition_Li = surface_lithium_composition
    composition_O = surface_oxygen_composition

    num = 0

    for i in composition_O:
        for j in composition_Li:
            if (i == 1 and j == 1) or (i == 0 and j == 0):
                continue
            if (i == 0) or (j == 0):
                try:
                    subs = SubstitutionTransformation(
                        {O_replacement: {O_replacement: i},
                         Li_replacement: {Li_replacement: j}})
                    surface_structure_partial = subs.apply_transformation(
                        slab_tgt)
                    enum = EnumerateStructureTransformation(
                        min_cell_size=target_cell_size,
                        max_cell_size=target_cell_size)
                    structures = enum.apply_transformation(
                        surface_structure_partial, return_ranked_list=2000)
                except AttributeError:
                    print(f'******************** Unable to enumerate the'
                          f' {j * 100}% lithium and {i * 100}% oxygen directly'
                          f'.********************')
                    continue
                else:
                    pass

            subs = SubstitutionTransformation(
                {O_replacement: {O_replacement: i},
                 Li_replacement: {Li_replacement: j}})
            surface_structure_partial = subs.apply_transformation(
                slab_tgt)
            enum = EnumerateStructureTransformation(
                min_cell_size=target_cell_size,
                max_cell_size=target_cell_size)
            structures = enum.apply_transformation(
                surface_structure_partial, return_ranked_list=2000)
            new_structures = []
            for k, s in enumerate(structures):
                # Keep volume constant and c lattice parameter unchanged
                if ((-1 < s['structure'].lattice.abc[2] - c < 1) and
                        (-1 < s['structure'].volume
                         - target_cell_size * volume < 1)):
                    new_structures.append(structures[k]['structure'])
            num += len(new_structures)

            print(f'The enumeration found {len(structures)} '
                  f'and {len(new_structures)} '
                  f'distinct structures for {j * 100}% Li and {i * 100}% O '
                  f'before and after filter respectively.')

            # Add selective dynamics for enumerated sites
            for structure in new_structures:
                for t in structure:
                    if Li_replacement in t:
                        t.properties = {
                            'selective_dynamics': [True, True, True]}
                    if O_replacement in t:
                        t.properties = {
                            'selective_dynamics': [True, True, True]}
                    if t.properties['selective_dynamics'] is None:
                        if layer_classification(slab_tgt)[0] - 0.01 \
                                <= t.frac_coords[2] \
                                <= layer_classification(slab_tgt)[1] + 0.01:
                            t.properties = {
                                'selective_dynamics': [False, False, False]}
                        else:
                            t.properties = {
                                'selective_dynamics': [True, True, True]}

            # Symmetrize structure models
            symmetrized_structures = []
            for s in new_structures:
                if j == 0:
                    s.replace_species({O_replacement: 'O'})
                if i == 0:
                    s.replace_species({Li_replacement: 'Li'})
                if (i and j) != 0:
                    s.replace_species({Li_replacement: 'Li',
                                       O_replacement: 'O'})

                symmetrized_structures.append(symmetrize_top_base(s))

            if to_vasp:
                # Generate structure models
                dirname = str(j) + 'Li' + str(i) + 'O'
                if not os.path.exists(dirname):
                    os.makedirs(dirname)
                for k, s in enumerate(symmetrized_structures):
                    s.to(fmt='poscar', filename=os.path.join(
                        dirname, "structure-{}.vasp".format(k)))

    print(f'{num} distinct structures are found totally.')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__ + "\n{} {}".format(__date__, __author__),
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(
        "input_file",
        help="Path to an input file in VASP's POSCAR format.")

    parser.add_argument(
        "--target-cell-size", "-s",
        help="Target cell size relative to input cell (default: 1).",
        type=int,
        default=1)

    parser.add_argument(
        "--oxygen-composition", "-o",
        nargs="+",
        type=float,
        help="All desired surface oxygen composition.",
        default=[1, 0.75, 0.5, 0.25, 0])

    parser.add_argument(
        "--lithium-composition", "-L",
        help="All desired surface lithium composition.",
        nargs="+",
        type=float,
        default=[1, 0.75, 0.5, 0.25, 0])

    parser.add_argument(
        "--generate-poscars", "-g",
        help="Generate POSCAR files of enumerated structures.",
        action="store_true")

    args = parser.parse_args()

    automate_surface(args.input_file,
                     target_cell_size=args.target_cell_size,
                     surface_oxygen_composition=args.oxygen_composition,
                     surface_lithium_composition=args.lithium_composition,
                     to_vasp=args.generate_poscars)
