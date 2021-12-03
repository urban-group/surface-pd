#!/usr/bin/env python

"""
This code will generate slab models with O and Li vacancies on the surface
    automatically. The input slab structure should be as small as
    possible because it will increase the possibilities it can have.
    automatically.
"""

__author__ = "Xinhao Li"
__email__ = "xl2778@columbia.edu"
__date__ = "2021-12-03"

import argparse
import os

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
                                     shift_inversion_symmetry_center)


def automate_surface(target_slab,
                     surface_oxygen_composition,
                     surface_lithium_composition,
                     target_cell_size,
                     num_layers_relaxed,
                     shift_isc=False,
                     to_vasp=False):
    """
    This function contains the general framework to  enumerate the parent
    slab model with different composition of oxygen and lithium vacancies.
    Args:
        target_slab: Fully lithiated slab with no oxygen vacancies on the
        surface
        surface_oxygen_composition: Desired composition of oxygen on the
        surface. This composition is determined by number of oxygen remaining
        on the surface divided by the total number of oxygen atoms can be
        removed.
        surface_lithium_composition: Desired composition of lithium on the
        surface. This composition is determined by number of lithium remaining
        on the surface divided by the total number of lithium atoms can be
        removed.
        target_cell_size: Maximum and minimum number of supercells of the
        input slab.
        num_layers_relaxed: Define how many layers on the surface will be
        relaxed during the geometry optimization. This is mainly useful for the
        (104) surface.
        shift_isc: Whether to shift the inversion symmetry center to
        (0, 0, 0). This argument is used to facilitate VASP to detect the
        symmetry of the slab model and probably accelerate the time to
        converge.
        to_vasp: Whether to generate all the output slab models in a
        well-organized format.

    Returns: All enumerated slabs with different lithium and oxygen
    vacancy compositions
    """

    surface_oxygen_composition.sort(reverse=True)
    surface_lithium_composition.sort(reverse=True)
    print("target_cell_size = {}".format(target_cell_size))
    print("Composition of oxygen on the surface will be {}.".format(
        surface_oxygen_composition))
    print("Composition of lithium on the surface will be {}.".format(
        surface_lithium_composition))

    # Load initial slab model with no oxygen and lithium vacancies on the
    # surface
    input_structure = Structure.from_file(target_slab)

    # Replace the top surface Li and O atoms with other species which will
    # be easier to enumerate using the enumlib
    # O_replacement = DummySpecies(symbol='X', oxidation_state=-2)
    # Li_replacement = DummySpecies(symbol='Z', oxidation_state=+1)
    O_replacement = Species(symbol='F', oxidation_state=-2)
    Li_replacement = Species(symbol='Na', oxidation_state=+1)
    slab_tgt = surface_substitute(target_slab,
                                  subs1=O_replacement,
                                  subs2=Li_replacement)

    # Extract c parameter and volume of the primitive unit cell (will be
    # used as criteria next)
    a, b, c = slab_tgt.lattice.abc

    # Enumerate with maximum unit cell of 4, but the cell size can also be
    # like 2x1 or 1x2 or 2x2
    composition_Li = surface_lithium_composition
    composition_O = surface_oxygen_composition
    if target_cell_size == 4:
        scaling_matrix = [2, 2, 1]
    else:
        if round(a / b, 0) == 2:
            scaling_matrix = [1, target_cell_size, 1]
        elif round(b / a, 0) == 2:
            scaling_matrix = [target_cell_size, 1, 1]
        else:
            scaling_matrix = [target_cell_size, 1, 1]
    print("Scaling matrix used here is: {}".format(scaling_matrix))

    # Initialize a number store total number of enumerated slab models
    num = 0

    # c fractional coordinates of the boundaries of fixed region (top and
    # bottom), surface metal, surface O
    center_bottom, center_top, surface_metal, surface_O, li_layers, layers \
        = layer_classification(input_structure)

    # Get indices of first and second layers which containing Li atoms,
    # surface oxygen atoms as well as the relaxed Li atoms in the slab model
    first_layer_index, second_layer_index, oxygen_index, relaxed_li_index = \
        index_extraction(input_structure)

    # Initialize a dictionary to store the successfully enumerated slab
    # models by composition
    enumerated_slabs_by_composition = {}

    # Begin to enumerate
    for i in composition_O:
        for j in composition_Li:
            prev = 0
            # create a list to store supplemental slabs
            supplemental_structures = []
            if i == 1 and j == 1:
                structure1 = input_structure.copy()
                structure1.make_supercell(scaling_matrix=scaling_matrix)
                supplemental_structures = [structure1]
                symmetrized_structures = [structure1]
            elif i == 0 and j == 1:
                structure1 = remove_sites(input_structure,
                                          index=oxygen_index,
                                          scaling_matrix=scaling_matrix)
                supplemental_structures = [structure1]
                symmetrized_structures = [structure1]
            elif i == 1 and j == 0:
                structure1 = remove_sites(input_structure,
                                          index=relaxed_li_index,
                                          scaling_matrix=scaling_matrix)
                supplemental_structures = [structure1]
                symmetrized_structures = [structure1]
            elif i == 0 and j == 0:
                structure1 = remove_sites(
                    input_structure,
                    index=oxygen_index + relaxed_li_index,
                    scaling_matrix=scaling_matrix)
                supplemental_structures = [structure1]
                symmetrized_structures = [structure1]
            else:
                # Check whether the not-working case is actually a
                # sub-case (remove all surface O atoms) of previous
                # successfully enumerated composition
                if str([i + 1, j]) \
                        in enumerated_slabs_by_composition.keys():
                    direct_structures = []
                    for struct in enumerated_slabs_by_composition[str(
                            [i + 1, j])]:
                        # Remove surface oxygen atoms by extracting
                        # their indices first
                        oxygen_index_enumed = index_extraction(struct)[2]
                        direct_structures.append(remove_sites(
                            structure_model=struct,
                            index=oxygen_index_enumed,
                            scaling_matrix=[1, 1, 1]))
                    prev = 0
                    if to_vasp:
                        dirname = str(j) + 'Li' + str(i) + 'O'
                        if not os.path.exists(dirname):
                            os.makedirs(dirname)
                        for k, s in enumerate(direct_structures):
                            s.to(fmt='poscar', filename=os.path.join(
                                dirname, "structure-{}.vasp".format(k)))
                    supplemental_structures = direct_structures
                    symmetrized_structures = direct_structures

                else:
                    structures = enum_with_composition(
                        slab_tgt,
                        subs_o=O_replacement,
                        o_composition=i,
                        subs_li=Li_replacement,
                        li_composition=j,
                        cell_size=target_cell_size)

                    if j == 0.5:
                        # Used to check whether it is a polar or non-polar
                        # surface.
                        if 0 <= abs(surface_metal - surface_O) <= 0.01:
                            #  these supplemental slabs are only needed for
                            #  104 surface
                            if i == 1:
                                structure1 = remove_sites(
                                    input_structure,
                                    index=first_layer_index,
                                    scaling_matrix=scaling_matrix)
                                structure2 = remove_sites(
                                    input_structure,
                                    index=second_layer_index,
                                    scaling_matrix=scaling_matrix)
                                supplemental_structures = [structure1,
                                                           structure2]
                            if i == 0:
                                structure1 = remove_sites(
                                    input_structure,
                                    index=first_layer_index + oxygen_index,
                                    scaling_matrix=scaling_matrix)
                                structure2 = remove_sites(
                                    input_structure,
                                    index=second_layer_index + oxygen_index,
                                    scaling_matrix=scaling_matrix)
                                supplemental_structures = [structure1,
                                                           structure2]
                    # Check whether supplemental structures are symmetry or not
                    for struct in supplemental_structures:
                        sga = SpacegroupAnalyzer(struct, symprec=0.1)
                        if not sga.is_laue():
                            print("The supplemental structures are not all "
                                  "symmetric!")

                    # Filtered out the structures which has c lattice as the
                    # largest lattice (a tall cuboid)
                    new_structures = []
                    for k, s in enumerate(structures):
                        lattice = s['structure'].lattice.abc
                        # Strict criteria -- keeps slabs with exact c lattice
                        # as parent slab models
                        if lattice[2] == c:
                            new_structures.append(structures[k]['structure'])
                    # In case of that using strict criteria will remove all
                    # structures, apply modest criteria to complete the
                    # dataset.
                    if len(new_structures) == 0:
                        print("**The criteria used is too strict! Changing to "
                              "the modest one.**")
                        for k, s in enumerate(structures):
                            lattice = s['structure'].lattice.abc
                            # Keep the c direction of the slab models is
                            # perpendicular to x-y plane but the c lattice
                            # parameter can be modified for a little bit.
                            if (lattice[0] and lattice[1]) < lattice[2] \
                                    <= c * 1.5:
                                new_structures.append(
                                    structures[k]['structure'])

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
                    for s in new_structures:
                        if i == 0 and j == 0:
                            pass
                        elif j == 0 and i != 0:
                            s.replace_species({O_replacement: 'O'})
                        elif i == 0 and j != 0:
                            s.replace_species({Li_replacement: 'Li'})
                        elif (i and j) != 0:
                            s.replace_species({Li_replacement: 'Li',
                                               O_replacement: 'O'})
                        symmetrized_structures.append(symmetrize_top_base(s))
                    prev = len(symmetrized_structures)

                    # Checking whether any supplemental_structure is redundant
                    # structure of the symmetrized_structures dataset. Also
                    # combined supplemental_structures dataset into the
                    # symmetrized_structures dataset.
                    if len(supplemental_structures) != 0:
                        new_added = []
                        for s1 in supplemental_structures:
                            same = []
                            for s2 in symmetrized_structures:
                                sm = StructureMatcher(primitive_cell=False)
                                fit = sm.fit(s1, s2)
                                same.append(fit)
                                if (s2.lattice.b - s2.lattice.a) < 2:
                                    reference = s2
                            if not any(same):
                                sm = StructureMatcher(primitive_cell=False,
                                                      ignored_species=["Li"])
                                transformed_structure = sm.get_s2_like_s1(
                                    struct1=reference,
                                    struct2=s1,
                                    include_ignored_species=True)
                                symmetrized_structures.append(
                                    transformed_structure)
                                new_added.append(0)
                            else:
                                supplemental_structures = new_added

                    # Adjust number of relaxed layers (if no need to adjust,
                    # just ignore this argument), mainly used for the (104)
                    # surface
                    fixed_layers_c = layers
                    relaxed_layers_c = []
                    if num_layers_relaxed != 0:
                        for index in range(0, num_layers_relaxed):
                            relaxed_layers_c.append(layers[index])
                            relaxed_layers_c.append(layers[-index - 1])
                        relaxed_layers_c.sort()
                        for item in relaxed_layers_c:
                            fixed_layers_c.remove(item)

                        # Add selective dynamics
                        for symmetrized_structure in symmetrized_structures:
                            for s in symmetrized_structure:
                                if (s.frac_coords[2] >= (
                                        fixed_layers_c[-1] + 0.02)
                                        or (s.frac_coords[2] <=
                                            fixed_layers_c[0] - 0.02)):
                                    s.properties = {
                                        'selective_dynamics': [True, True,
                                                               True]}
                                else:
                                    s.properties = {
                                        'selective_dynamics': [False, False,
                                                               False]}
                    # Shift the inversion symmetry center
                    if shift_isc:
                        shifted_structure = []
                        for struct in symmetrized_structures:
                            shifted_structure.append(
                                shift_inversion_symmetry_center(struct))
                    else:
                        shifted_structure = symmetrized_structures

                    # Add composition and enumerated structures to the
                    # dictionary as the keys and values
                    enumerated_slabs_by_composition[str([i, j])] \
                        = shifted_structure

                    # Generate slab models
                    if to_vasp:
                        dirname = str(j) + 'Li' + str(i) + 'O'
                        if not os.path.exists(dirname):
                            os.makedirs(dirname)
                        for k, s in enumerate(shifted_structure):
                            s.to(fmt='poscar', filename=os.path.join(
                                dirname, "structure-{}.vasp".format(k)))

            # Calculate total number of enumerated structures
            num += len(symmetrized_structures)

    # Print the results
            print(f'The enumeration found {len(symmetrized_structures)} '
                  f'({prev}+{len(supplemental_structures)}) '
                  f'distinct structures for {j * 100}% Li and {i * 100}% O.')
    print(f'{num} distinct structures are found totally.')
    if shift_isc:
        print(f'The isc of {num} structures are all shifted.')


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
        default=[1.0, 0.75, 0.5, 0.25, 0.0])

    parser.add_argument(
        "--lithium-composition", "-L",
        help="All desired surface lithium composition.",
        nargs="+",
        type=float,
        default=[1.0, 0.75, 0.5, 0.25, 0.0])

    parser.add_argument(
        "--shift-isc", "-m",
        help="Whether to shift isc (inversion symmetry center) to the "
             "origin.",
        action="store_true")

    parser.add_argument(
        "--num-of-relaxed-layers", "-r",
        help="Number of layers that will be relaxed on the surface.",
        type=int,
        default=0
    )

    parser.add_argument(
        "--generate-poscars", "-g",
        help="Generate POSCAR files of enumerated structures.",
        action="store_true")

    args = parser.parse_args()

    automate_surface(args.input_file,
                     target_cell_size=args.target_cell_size,
                     surface_oxygen_composition=args.oxygen_composition,
                     surface_lithium_composition=args.lithium_composition,
                     shift_isc=args.shift_isc,
                     num_layers_relaxed=args.num_of_relaxed_layers,
                     to_vasp=args.generate_poscars)