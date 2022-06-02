#!/usr/bin/env python

"""
This code will generate slab models with Li and O vacancies on the surface
    automatically. The input slab structure should be as small as
    possible because it will increase the number of possibilities that the
    enumeration can have.
"""

__author__ = "Xinhao Li"
__email__ = "xl2778@columbia.edu"
__date__ = "2022-03-22"

import argparse
import os

import numpy as np

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.periodic_table import Species
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from surface_pd.surface_enum import (surface_substitute,
                                     layer_classification,
                                     enum_with_composition,
                                     index_extraction,
                                     remove_sites,
                                     symmetrize_top_base,
                                     get_num_sites,
                                     slab_size_check,
                                     final_check,
                                     define_scaling_matrix,
                                     get_max_min_c_frac,
                                     add_selective_dynamics,
                                     temp_shift_isc_back)


def automate_surface(target_slab,
                     surface_lithium_composition,
                     surface_oxygen_composition,
                     target_cell_size,
                     num_layers_relaxed,
                     to_vasp=False):
    """
    ToDo: generalize the composition of surface atoms (not just Li and O)
    This function contains the general framework to enumerate the parent
    slab model with different composition of lithium and oxygen vacancies.
    Args:
        target_slab: Fully lithiated slab with no oxygen vacancies on the
        surface.
        surface_lithium_composition: Desired composition of lithium on the
        surface. This composition is determined by number of lithium
        remaining on the surface divided by the total number of lithium atoms
        that will be relaxed on the surface.
        surface_oxygen_composition: Desired composition of oxygen on the
        surface.
        target_cell_size: Maximum number of supercells of the input slab.
        num_layers_relaxed: Define how many layers on the surface will be
        relaxed during the geometry optimization.
        to_vasp: Whether to generate all the output slab models in a
        well-organized format.

    Returns: All enumerated slabs with different lithium and oxygen
    vacancy compositions
    """

    surface_oxygen_composition.sort(reverse=True)
    surface_lithium_composition.sort(reverse=True)
    print("target_cell_size = {}".format(target_cell_size))
    print("Composition of lithium on the surface will be {}.".format(
        surface_lithium_composition))
    print("Composition of oxygen on the surface will be {}.".format(
        surface_oxygen_composition))

    # Load initial slab model with no lithium and oxygen vacancies on the
    # surface
    input_structure = Structure.from_file(target_slab)

    # Replace the top surface Li and O atoms with other species which will
    # be easier to enumerate using the enumlib
    Li_replacement = Species(symbol='Na', oxidation_state=+1)
    O_replacement = Species(symbol='F', oxidation_state=-2)
    slab_substituted = surface_substitute(target_slab,
                                          subs1=O_replacement,
                                          subs2=Li_replacement)

    # Extract the c parameter of the parent slab model (will be used as
    # criteria next)
    a, b, c = slab_substituted.lattice.abc

    # Enumerate with maximum unit cell of 4, but the cell size can also be
    # like 2x1 or 1x2 or 2x2
    composition_Li = surface_lithium_composition
    composition_O = surface_oxygen_composition
    scaling_matrix = define_scaling_matrix(a, b, target_cell_size)
    print("Scaling matrix used here is: {}".format(scaling_matrix))

    # Initialize a number store total number of enumerated slab models
    num = 0

    # 0, 1 - c_fractional coordinates of the Lower and upper boundaies of the
    # fixed region in the center
    [center_bottom, center_top,
     _, _, _, _, _, _] = layer_classification(input_structure)

    # Get indices of the surface oxygen atoms as well as the relaxed Li atoms
    # in the slab model
    [_, _,
     oxygen_index, relaxed_li_index] = index_extraction(input_structure)

    # Initialize a dictionary to store the successfully enumerated slab
    # models by composition
    enumerated_slabs_by_composition = {}

    # Begin to enumerate
    for i in composition_Li:
        for j in composition_O:
            prev = 0
            # Create a list to store supplemental slabs, this mainly
            # contains edge cases
            supplemental_structures = []
            if i == 1 and j == 1:  # Fully lithiated
                structure1 = input_structure.copy()
                structure1.make_supercell(scaling_matrix=scaling_matrix)
                supplemental_structures = [structure1]
                symmetrized_structures = [structure1]
            elif i == 0 and j == 1:  # remove all surface Li atoms
                structure1 = remove_sites(input_structure,
                                          index=relaxed_li_index,
                                          scaling_matrix=scaling_matrix)
                supplemental_structures = [structure1]
                symmetrized_structures = [structure1]
            elif i == 1 and j == 0:  # remove all surface O atoms
                structure1 = remove_sites(input_structure,
                                          index=oxygen_index,
                                          scaling_matrix=scaling_matrix)
                supplemental_structures = [structure1]
                symmetrized_structures = [structure1]
            elif i == 0 and j == 0:  # remove all surface atoms
                structure1 = remove_sites(
                    input_structure,
                    index=relaxed_li_index + oxygen_index,
                    scaling_matrix=scaling_matrix)
                supplemental_structures = [structure1]
                symmetrized_structures = [structure1]
            else:
                # When all surface oxygen atoms are removed, the slab models
                # can not be enumerated, this might be the bug in the enumlib
                # library or pymatgen does no implement it properly. But
                # luckly, all the not-working cases are actually the
                # sub-cases(remove all surface O atoms) of previous
                # successfully enumerated composition
                if str([i, j + 1]) in enumerated_slabs_by_composition.keys():
                    direct_structures = []
                    for symmetrized_structure in \
                            enumerated_slabs_by_composition[str([i, j + 1])]:
                        # Remove surface oxygen atoms by extracting
                        # their indices first
                        oxygen_index_enumed = index_extraction(
                            symmetrized_structure)[2]
                        direct_structures.append(
                            remove_sites(structure_model=symmetrized_structure,
                                         index=oxygen_index_enumed,
                                         scaling_matrix=[1, 1, 1]))
                    supplemental_structures = direct_structures
                    symmetrized_structures = direct_structures

                else:
                    structures = enum_with_composition(
                        slab_substituted,
                        subs_li=Li_replacement,
                        li_composition=i,
                        subs_o=O_replacement,
                        o_composition=j,
                        cell_size=target_cell_size)

                    # Filtered out the structures which has c lattice as the
                    # largest lattice (a tall cuboid)
                    filtered_structures = []
                    for k, structure in enumerate(structures):
                        lattice = structure['structure'].lattice.abc
                        # Strict criteria -- keeps slabs with c lattice
                        # as parent slab models
                        if (round(lattice[2], 3) - 0.2) <= c <= \
                                (round(lattice[2], 3) + 0.2):
                            filtered_structures.append(
                                structures[k]['structure'])
                    # In case of that using strict criteria will remove all
                    # structures, apply modest criteria to complete the
                    # dataset.
                    if len(filtered_structures) == 0:
                        # print("**The criteria used is too strict! Changing to "
                        #       "the modest one.**")
                        for k, structure in enumerate(structures):
                            lattice = structure['structure'].lattice.abc
                            # Keep the c direction of the slab models is
                            # perpendicular to x-y plane but the c lattice
                            # parameter can be modified for a little bit.
                            if (lattice[0] and lattice[1]) < lattice[2] \
                                    <= c * 1.5:
                                filtered_structures.append(
                                    structures[k]['structure'])

                    # Add selective dynamics for enumerated sites
                    for filtered_structure in filtered_structures:
                        for t in filtered_structure:
                            if Li_replacement in t:
                                t.properties = {
                                    'selective_dynamics': [True, True, True]
                                }
                            if O_replacement in t:
                                t.properties = {
                                    'selective_dynamics': [True, True, True]
                                }
                            if t.properties['selective_dynamics'] is None:
                                if (center_bottom - 0.01 <= t.frac_coords[2]
                                        <= center_top + 0.01):
                                    t.properties = {
                                        'selective_dynamics':
                                            [False, False, False]}
                                else:
                                    t.properties = {
                                        'selective_dynamics':
                                            [True, True, True]}

                    # Symmetrize slab models based on the top enumerated
                    # surface
                    symmetrized_structures = []
                    for filtered_structure in filtered_structures:
                        if i == 0 and j == 0:
                            pass
                        elif i == 0 and j != 0:
                            filtered_structure.replace_species(
                                {O_replacement: 'O'})
                        elif i != 0 and j == 0:
                            filtered_structure.replace_species(
                                {Li_replacement: 'Li'})
                        elif (i and j) != 0:
                            filtered_structure.replace_species(
                                {Li_replacement: 'Li', O_replacement: 'O'})
                        symmetrized_structures.append(
                            symmetrize_top_base(filtered_structure))
                    prev = len(symmetrized_structures)

            # Add composition and enumerated structures to the
            # dictionary as the keys and values
            enumerated_slabs_by_composition[str([i, j])] \
                = symmetrized_structures

            for k, symmetrized_structure in enumerate(
                    symmetrized_structures):
                sga = SpacegroupAnalyzer(symmetrized_structure,
                                         symprec=0.1)

                # Create the refined structures
                refined_structure = sga.get_refined_structure()

                # Define number of sites that should be after enumeration
                total_num_sites = get_num_sites(input_structure,
                                                slab_substituted,
                                                target_cell_size,
                                                i,
                                                j)
                num_sites = refined_structure.num_sites

                indicator = 0
                while ((num_sites > total_num_sites) or
                       (max(refined_structure.lattice.abc) > c * 2 - 5) or
                       (max(refined_structure.lattice.abc) !=
                        refined_structure.lattice.c)):
                    indicator, refined_structure = slab_size_check(
                        refined_structure,
                        total_num_sites,
                        num_sites,
                        c)
                    num_sites = refined_structure.num_sites

                # if num_sites < total_num_sites:
                #     if total_num_sites % num_sites != 0:
                #         raise ValueError("Check primitive structure finder.")
                #     indicator = 1

                # Add selective dynamics
                min_c, max_c = get_max_min_c_frac(refined_structure)
                if max_c - min_c > 0.8:
                    refined_structure = temp_shift_isc_back(
                        before_refine_structure=symmetrized_structure,
                        after_refine_structure=refined_structure,
                        shift=True
                    )
                    refined_structure = add_selective_dynamics(
                        input_structure,
                        refined_structure,
                        num_relaxed=num_layers_relaxed
                    )
                    refined_structure = temp_shift_isc_back(
                        before_refine_structure=symmetrized_structure,
                        after_refine_structure=refined_structure,
                        shift=False
                    )
                else:
                    refined_structure = add_selective_dynamics(
                        input_structure,
                        refined_structure,
                        num_relaxed=num_layers_relaxed
                    )

                # Perform final check
                final_check(refined_structure, c, i, j, k)

                # Generate slab models
                if to_vasp:
                    # Make directories
                    dirname = str(i) + 'Li' + str(j) + 'O'
                    if not os.path.exists(dirname):
                        os.makedirs(dirname)

                    if indicator == -1:
                        refined_structure.to(
                            fmt='poscar',
                            filename=os.path.join(
                                dirname,
                                "refined-prim-{}.vasp".format(k)))
                    elif indicator == 1:
                        refined_structure.to(
                            fmt='poscar',
                            filename=os.path.join(
                                dirname,
                                "refined-reduced-{}.vasp".format(k)))
                    elif indicator == 2:
                        refined_structure.to(
                            fmt='poscar',
                            filename=os.path.join(
                                dirname,
                                "refined-rotated-{}.vasp".format(k)))
                    else:
                        refined_structure.to(
                            fmt='poscar',
                            filename=os.path.join(
                                dirname,
                                "refined-structure-{}.vasp".format(k)))

            # Calculate total number of enumerated structures
            num += len(symmetrized_structures)

            # Print the results
            print(f'The enumeration found {len(symmetrized_structures)} '
                  f'({prev}+{len(supplemental_structures)}) '
                  f'distinct structures for {i * 100}% Li and {j * 100}% O.')

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
        "--lithium-composition", "-L",
        help="All desired surface lithium composition.",
        nargs="+",
        type=float,
        default=[1.0, 0.75, 0.5, 0.25, 0.0])

    parser.add_argument(
        "--oxygen-composition", "-o",
        nargs="+",
        type=float,
        help="All desired surface oxygen composition.",
        default=[1.0, 0.75, 0.5, 0.25, 0.0])

    parser.add_argument(
        "--num-of-relaxed-layers", "-nr",
        help="Number of layers that will be relaxed on the surface.",
        type=int,
        default=2
    )

    parser.add_argument(
        "--generate-poscar", "-g",
        help="Generate POSCAR files of enumerated structures.",
        action="store_true")

    args = parser.parse_args()

    automate_surface(args.input_file,
                     surface_lithium_composition=args.lithium_composition,
                     surface_oxygen_composition=args.oxygen_composition,
                     target_cell_size=args.target_cell_size,
                     num_layers_relaxed=args.num_of_relaxed_layers,
                     to_vasp=args.generate_poscar)
