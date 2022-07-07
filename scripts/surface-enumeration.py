#!/usr/bin/env python

"""
This code will generate slab models with Li and O vacancies on the surface
    automatically. The input slab structure should be as small as
    possible because it will increase the number of possibilities that the
    enumeration can have.
"""

__author__ = "Xinhao Li"
__email__ = "xl2778@columbia.edu"
__date__ = "2022-06-01"

import os
import string
import argparse

import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Species, DummySpecies
from pymatgen.analysis.structure_matcher import StructureMatcher

from surface_pd.core import Slab, EnumWithComposition, PreCheck, PostCheck
from surface_pd.analysis.slab_analysis import (get_num_sites,
                                               add_selective_dynamics)
from surface_pd.util import (define_scaling_matrix, temp_shift_isc_back)
from surface_pd.error import *


def automate_surface(target_slab_path,
                     surface_lithium_composition,
                     surface_oxygen_composition,
                     target_cell_size,
                     num_layers_relaxed,
                     to_vasp=False):
    """
    This function contains the general framework to enumerate the parent
    slab model with different composition of lithium and oxygen vacancies.
    Args:
        target_slab_path: Fully lithiated slab with no oxygen vacancies on the
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

    Returns: All enumerated slabs with different composition of lithium and
    oxygen vacancies
    """
    surface_lithium_composition.sort(reverse=True)
    surface_oxygen_composition.sort(reverse=True)

    print("target_cell_size = {}".format(target_cell_size))
    print("Composition of lithium on the surface will be {}.".format(
        surface_lithium_composition))
    print("Composition of oxygen on the surface will be {}.".format(
        surface_oxygen_composition))

    # Load initial slab model with no lithium and oxygen vacancies on the
    # surface
    input_structure = Slab.from_file(target_slab_path)

    # PreCheck
    precheck = PreCheck(input_structure)

    if precheck.is_not_cuboid():
        raise SlabOrientationError

    if precheck.is_not_slab():
        raise NonSlabError

    if precheck.not_all_has_selective_dynamics():
        raise NonDefinedSelectiveDynamicsError

    if precheck.has_no_inversion_symmetry():
        raise NoInversionSymmetryError


    # Define the dummy species that will be used to substitute the target
    # species
    Li_replacement = DummySpecies(symbol='A', oxidation_state=+1)
    O_replacement = DummySpecies(symbol='D', oxidation_state=-2)

    # Replace the top surface Li and O atoms with dummy species which will
    # be easier to enumerate using the enumlib
    slab_substituted = input_structure.surface_substitute(subs1=Li_replacement,
                                                          subs2=O_replacement)

    # Extract the c parameter of the parent slab model (will be used as
    # criteria next)
    a, b, c = slab_substituted.lattice.abc

    # Enumerate with maximum unit cell of 4, but the cell size can also be
    # like 2x1 or 1x2 or 2x2
    composition_Li = surface_lithium_composition
    composition_O = surface_oxygen_composition
    scaling_matrix = define_scaling_matrix(a, b, target_cell_size)
    print("Scaling matrix used here is: {}".format(scaling_matrix))

    # Initialize a number to store total number of enumerated slab models
    num = 0

    # c_fractional coordinates of the lower and upper boundaries of the
    # fixed region in the centeral slab
    [center_bottom, center_top,
     _, _, _, _, _, _] = input_structure.layer_distinguisher()

    # Get indices of the relaxed Li atoms as well as the surface oxygen atoms
    # in the slab model
    relaxed_Li_index, oxygen_index = input_structure.index_extraction()

    # Initialize a dictionary to store the successfully enumerated slab
    # models by composition
    enumerated_slabs_by_composition = {}

    # Begin to enumerate
    for i in composition_Li:
        for j in composition_O:
            prev = 0
            # Create a list to store supplemental slabs, this mainly
            # contains slabs generated from the edge cases
            supplemental_structures = []
            if i == 1 and j == 1:  # Fully lithiated
                structure1 = input_structure.copy()
                structure1.make_supercell(scaling_matrix=scaling_matrix)
                supplemental_structures = [structure1]
                symmetrized_structures = [structure1]
            elif i == 0 and j == 1:  # remove all surface Li atoms
                structure1 = input_structure.remove_sites_with_scaling(
                    index=relaxed_Li_index,
                    scaling_matrix=scaling_matrix)
                # structure1 = remove_sites(input_structure,
                #                           index=relaxed_li_index,
                #                           scaling_matrix=scaling_matrix)
                supplemental_structures = [structure1]
                symmetrized_structures = [structure1]
            elif i == 1 and j == 0:  # remove all surface O atoms
                structure1 = input_structure.remove_sites_with_scaling(
                    index=oxygen_index,
                    scaling_matrix=scaling_matrix
                )
                # structure1 = remove_sites(input_structure,
                #                           index=oxygen_index,
                #                           scaling_matrix=scaling_matrix)
                supplemental_structures = [structure1]
                symmetrized_structures = [structure1]
            elif i == 0 and j == 0:  # remove all surface atoms
                structure1 = input_structure.remove_sites_with_scaling(
                    index=relaxed_Li_index + oxygen_index,
                    scaling_matrix=scaling_matrix
                )
                # structure1 = remove_sites(
                #     input_structure,
                #     index=relaxed_li_index + oxygen_index,
                #     scaling_matrix=scaling_matrix)
                supplemental_structures = [structure1]
                symmetrized_structures = [structure1]
            else:
                # When all surface oxygen atoms are removed, the slab models
                # can not be enumerated. This might be the bug in the enumlib
                # library source code or pymatgen does no implement it
                # properly. But luckly, all the not-working cases are
                # actually the sub-cases(remove all surface O atoms) of
                # previous successfully enumerated composition
                if str([i, j + 1]) in enumerated_slabs_by_composition.keys():
                    direct_structures = []
                    for enumed_structure in \
                            enumerated_slabs_by_composition[str([i, j + 1])]:
                        # Remove surface oxygen atoms by extracting
                        # their indices first
                        oxygen_index_enumed = \
                            enumed_structure.index_extraction()[1]
                        # oxygen_index_enumed = index_extraction(
                        #     symmetrized_structure)[2]
                        direct_structures.append(
                            enumed_structure.remove_sites_with_scaling(
                                index=oxygen_index_enumed,
                                scaling_matrix=[1, 1, 1]
                            )
                        )
                        # direct_structures.append(
                        #     remove_sites(structure_model=symmetrized_structure,
                        #                  index=oxygen_index_enumed,
                        #                  scaling_matrix=[1, 1, 1]))
                    supplemental_structures = direct_structures
                    symmetrized_structures = direct_structures

                else:
                    # New
                    ewc = EnumWithComposition(
                        subs_species=[Li_replacement, O_replacement],
                        subs_composition=[i, j],
                        max_cell_size=target_cell_size
                    )
                    structures = ewc.apply_enumeration(slab_substituted)

                    # structures = enum_with_composition(
                    #     slab_substituted,
                    #     subs_li=Li_replacement,
                    #     li_composition=i,
                    #     subs_o=O_replacement,
                    #     o_composition=j,
                    #     cell_size=target_cell_size)

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
                        elif i != 0 and j == 0:
                            filtered_structure.replace_species(
                                {Li_replacement: 'Li'})
                        elif i == 0 and j != 0:
                            filtered_structure.replace_species(
                                {O_replacement: 'O'})
                        elif (i and j) != 0:
                            filtered_structure.replace_species(
                                {Li_replacement: 'Li', O_replacement: 'O'})
                        filtered_structure = \
                            Slab.from_sites(filtered_structure)
                        symmetrized_structures.append(
                            filtered_structure.symmetrize_top_base()
                        )
                        # symmetrized_structures.append(
                        #     symmetrize_top_base(filtered_structure))
                    prev = len(symmetrized_structures)

            # Add composition and enumerated structures to the
            # dictionary as the keys and values
            enumerated_slabs_by_composition[str([i, j])] \
                = symmetrized_structures

            for k, symmetrized_structure in enumerate(symmetrized_structures):
                sga = SpacegroupAnalyzer(symmetrized_structure,
                                         symprec=0.1)

                # Create the refined structures
                refined_structure = sga.get_refined_structure()

                # Define number of sites that should be after enumeration
                total_num_sites = get_num_sites(input_structure,
                                                slab_substituted,
                                                Li_replacement,
                                                i,
                                                O_replacement,
                                                j,
                                                target_cell_size)
                num_sites = refined_structure.num_sites

                indicator = 0
                while ((num_sites > total_num_sites) or
                       (max(refined_structure.lattice.abc) > c * 2 - 5) or
                       (max(refined_structure.lattice.abc) !=
                        refined_structure.lattice.c)):
                    pc = PostCheck(refined_structure)
                    indicator, refined_structure = pc.slab_size_check(
                        total_num_sites=total_num_sites,
                        enumerated_num_sites=num_sites,
                        input_c=c
                    )
                    num_sites = refined_structure.num_sites

                # Add selective dynamics
                refined_structure = Slab.from_sites(refined_structure)
                min_c, max_c = refined_structure.get_max_min_c_frac()
                # min_c, max_c = get_max_min_c_frac(refined_structure)
                if max_c - min_c > 0.8:
                    refined_structure = temp_shift_isc_back(
                        before_refined_structure=symmetrized_structure,
                        after_refined_structure=refined_structure,
                        shift=True
                    )
                    refined_structure = add_selective_dynamics(
                        input_structure,
                        refined_structure,
                        num_relaxed=num_layers_relaxed
                    )
                    refined_structure = temp_shift_isc_back(
                        before_refined_structure=symmetrized_structure,
                        after_refined_structure=refined_structure,
                        shift=False
                    )
                else:
                    refined_structure = add_selective_dynamics(
                        input_structure,
                        refined_structure,
                        num_relaxed=num_layers_relaxed
                    )

                # Perform final check
                pc = PostCheck(refined_structure)
                pc.final_check(
                    Li_composition=i,
                    O_composition=j,
                    index=k
                )

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
        help="All desired surface lithium composition "
             "(default: [[1.0, 0.75, 0.5, 0.25, 0.0]]).",
        nargs="+",
        type=float,
        default=[1.0, 0.75, 0.5, 0.25, 0.0])

    parser.add_argument(
        "--oxygen-composition", "-O",
        nargs="+",
        type=float,
        help="All desired surface oxygen composition "
             "(default: [[1.0, 0.75, 0.5, 0.25, 0.0]]).",
        default=[1.0, 0.75, 0.5, 0.25, 0.0])

    parser.add_argument(
        "--num-of-relaxed-layers", "-nr",
        help="Number of layers that will be relaxed on the surface "
             "(default: 2).",
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
