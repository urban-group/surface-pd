#!/usr/bin/env python

"""
This code will enumerate the input slab model with user defined target
species and composition. The input slab structure should be as small as
possible because the unit cell of the slab will also be enumerated, and it
will increase the number of possibilities that the enumeration can have.
For the detailed algorithm behind the enumeration, please see the
following references below.
(1) Morgan, W. S.; Hart, G. L. W.; Forcade, R. W.
Computational Materials Science 2017, 136, 144–149.
https://doi.org/10.1016/j.commatsci.2017.04.015.
(2) Hart, G. L. W.; Nelson, L. J.; Forcade, R. W.
Computational Materials Science 2012, 59, 101–107.
https://doi.org/10.1016/j.commatsci.2012.02.015.
(3) Hart, G. L. W.; Forcade, R. W.
Phys. Rev. B 2008, 77 (22), 224115.
https://doi.org/10.1103/PhysRevB.77.224115.
(4) Hart, G. L. W.; Forcade, R. W.
Phys. Rev. B 2009, 80 (1), 014120.
https://doi.org/10.1103/PhysRevB.80.014120.

"""

__author__ = "Xinhao Li"
__email__ = "xinhao.li@columbia.edu"
__date__ = "2022-09-06"

import os
import json
import copy
import argparse
import warnings

from itertools import product

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import DummySpecies
from pymatgen.core.surface import get_slab_regions

from surface_pd.core import Slab, EnumWithComposition, PreCheck, PostCheck
from surface_pd.analysis.slab_analysis import (structure_filter,
                                               selective_dynamics_completion)
from surface_pd.util import (get_values_nested_dict,
                             all_int, replace_dummy, have_zero)
from surface_pd.error import *

warnings.filterwarnings("ignore")


def automate_surface(target_slab_path: str,
                     species: list,
                     replace: list,
                     num_layers_relaxed: dict,
                     max_cell_size: int,
                     symmetric: bool,
                     max_structures: int,
                     to_vasp: bool = False):
    """
    ToDo:
        1. Li_100 slab works, but any atom with c frac < 0 has issue to
        assign the selective dynamics.
        2. Re-test whether the new framework works for all cases.
        3. Check the documentation of PostCheck class.

    This function contains the general framework to enumerate the parent
    slab model with different target species and compositions.

    Args:
        target_slab_path: Input slab model
        species: Target species that will be enumerated.
        replace: Species and occupancy dictionaries containing the species
            mapping in string-string pairs. E.g. {'Li': {'Li': 0.5}}
            (This dict represents that only half of the original Li
            atoms in the slab model will be kept), stored in the list.
        num_layers_relaxed: Dictionaries where species and number of relaxed
            layers as keys and values, respectively. E.g. {'Li': 2, 'O': 1},
            which represents 2 and 1 layers of Li and O will be relaxed
            during the DFT geometry optimization.
        max_cell_size: Maximum number of supercells of the input slab.
        symmetric: Whether the symmetry of the slab model will be kept
            after the enumeration.
        max_structures: Number of structures to be returned at most for
            each composition.
        to_vasp: Whether to generate all the output slab models in a
            well-organized format. Defaults to False.

    Returns:
        All enumerated slab models with different composition of target
        species.

    """

    # Load initial slab model
    input_structure = Slab.from_file(target_slab_path)
    input_structure.wrap_pbc()

    # Define important parameters which will be used later, including to be
    # enumerated species, number of layers relaxed, whether to keep the
    # symmetry, direction of the slab, and the proximity tolerance for
    # adjacent atoms
    input_structure.to_be_enumerated_species = species
    input_structure.num_layers_relaxed = num_layers_relaxed
    input_structure.symmetric = symmetric
    direction = input_structure.direction
    tolerance = input_structure.tolerance

    # PreCheck
    pre_check = PreCheck(input_structure)

    if not pre_check.is_cuboid():
        raise SlabOrientationError

    if not pre_check.is_slab():
        raise NonSlabError

    if not pre_check.all_has_selective_dynamics():
        raise NonDefinedSelectiveDynamicsError

    if pre_check.relax_both_surfaces() != input_structure.symmetric:
        raise IncompatibleSymmError

    if symmetric:
        if not pre_check.has_inversion_symmetry():
            raise NoInversionSymmetryError

    # Check the user defined to-be-enumerated composition is valid
    if not pre_check.has_validate_composition(replace, max_cell_size):
        raise InvalidCompositionError

    # Extract the lattice parameter of the parent slab model along the slab
    # direction defined by user (will be used as criteria next)
    criteria = input_structure.lattice.abc[direction]

    # Get the indices of the to-be-enumerated atoms on the surface.
    # Get the c_frac coordinates of the lower and upper boundaries of the
    # fixed region in the central slab.
    center_bottom, center_top, relaxed_index = \
        input_structure.index_extraction(only_top=False)

    # print(relaxed_index)
    # Define the "dummy" species that will be used to substitute the target
    # species
    dummy_species = [DummySpecies(symbol='X' + str(i + 1))
                     for i in range(len(species))]

    # Replace the to-be-enumerated species from the slab top surface with
    # "dummy" species. This will facilitate the enumeration code.
    input_substituted = input_structure.surface_substitute(dummy_species)

    # Initialize a number to store total number of enumerated slab models
    num = 0

    # Begin to enumerate
    for subs_dict in replace:
        slab_substituted = copy.deepcopy(input_substituted)

        prev = 0
        supplemental_structures = []

        composition_list = list(get_values_nested_dict(subs_dict))
        if all_int(composition_list):
            structure = input_structure.supplemental_structures_gene(
                subs_dict, relaxed_index)
            supplemental_structures = [structure]
            enumerated_structures = [structure]
        else:
            if have_zero(composition_list):
                slab_substituted = slab_substituted \
                    .supplemental_structures_gene(subs_dict, relaxed_index)

            subs_dict = replace_dummy(subs_dict, dummy_species)

            ewc = EnumWithComposition(
                subs_dict=subs_dict,
                max_cell_size=max_cell_size
            )

            structures = ewc.apply_enumeration(slab_substituted,
                                               max_structures)

            # Filtered out the structures which has c lattice as the
            # largest lattice (a tall cuboid)
            filtered_structures = structure_filter(structures,
                                                   direction=direction,
                                                   criteria=criteria)
            # print(filtered_structures)
            # Add selective dynamics for enumerated sites
            for filtered_structure in filtered_structures:
                # print(filtered_structure)
                selective_dynamics_completion(
                    structure=filtered_structure,
                    direction=direction,
                    dummy_species=dummy_species,
                    center_bottom=center_bottom,
                    center_top=center_top,
                    tolerance=tolerance
                )

            # Symmetrize slab models based on the top enumerated
            # surface if needed
            enumerated_structures = []
            for filtered_structure in filtered_structures:
                for i in range(len(species)):
                    filtered_structure.replace_species(
                        {dummy_species[i]: species[i]})
                filtered_structure = Slab.from_sites(filtered_structure)
                if symmetric:
                    enumerated_structures.append(
                        filtered_structure.symmetrize_top_base())
                else:
                    enumerated_structures.append(filtered_structure)

            prev = len(enumerated_structures)

        unique_index = "-{:0%dd}" % (len(str(len(enumerated_structures))))
        for k, symmetrized_structure in enumerate(enumerated_structures):
            indicator = 0
            if symmetric:
                _, origin, _ = symmetrized_structure.is_symmetry(
                    return_isc=True)
                min_c, _ = symmetrized_structure.get_max_min_c_frac()
                sga = SpacegroupAnalyzer(symmetrized_structure,
                                         symprec=0.1)
                # Create the refined structures
                refined_structure = sga.get_refined_structure()

                # Define number of sites that should be after enumeration
                total_num_sites = input_structure.calculate_num_sites(
                    composition_list=composition_list,
                    relaxed_index=relaxed_index,
                    max_cell_size=max_cell_size)

                num_sites = refined_structure.num_sites

                while ((num_sites > total_num_sites) or
                       (max(refined_structure.lattice.abc) >
                        criteria * 2 - 5) or
                       (max(refined_structure.lattice.abc) !=
                        refined_structure.lattice.abc[direction])):
                    pc = PostCheck(refined_structure)
                    indicator, refined_structure = pc.slab_size_check(
                        total_num_sites=total_num_sites,
                        enumerated_num_sites=num_sites,
                        criteria=criteria)
                    num_sites = refined_structure.num_sites

                refined_structure = Slab.from_sites(refined_structure)
                noncontiguous = (len(get_slab_regions(refined_structure)) == 2)
                # print(input_structure.layers_finder())

                if noncontiguous:
                    refined_structure = refined_structure.tune_isc(
                        origin=origin, shift_isc_back=True)
                    refined_structure = refined_structure.tune_c(
                        target_min_c=min_c)
                    # print(refined_structure.layers_finder())
                    refined_structure = \
                        refined_structure.add_selective_dynamics(
                            lower_limit=center_bottom,
                            upper_limit=center_top
                        )
                else:
                    refined_structure = refined_structure.tune_c(
                        target_min_c=min_c)
                    # print(refined_structure.layers_finder())
                    refined_structure = \
                        refined_structure.add_selective_dynamics(
                            lower_limit=center_bottom,
                            upper_limit=center_top
                        )

                refined_structure = refined_structure.tune_isc(
                    origin=origin, shift_isc_back=False)
            else:
                refined_structure = symmetrized_structure

            # Perform final check
            pc = PostCheck(refined_structure)
            pc.final_check(
                species=species,
                composition_list=composition_list,
                keep_symmetric=symmetric,
                criteria=criteria,
                index=k)

            # Generate slab models
            if to_vasp:
                # Make directories
                dirname = ''
                for n in range(len(species)):
                    dirname += str(composition_list[n]) + str(species[n])
                if not os.path.exists(dirname):
                    os.makedirs(dirname)

                if indicator == -1:
                    refined_structure.to(
                        fmt='poscar',
                        filename=os.path.join(
                            dirname,
                            "refined-prim{}.vasp".format(
                                unique_index.format(k))))
                elif indicator == 1:
                    refined_structure.to(
                        fmt='poscar',
                        filename=os.path.join(
                            dirname,
                            "refined-reduced{}.vasp".format(
                                unique_index.format(k))))
                elif indicator == 2:
                    refined_structure.to(
                        fmt='poscar',
                        filename=os.path.join(
                            dirname,
                            "refined-rotated{}.vasp".format(
                                unique_index.format(k))))
                else:
                    refined_structure.to(
                        fmt='poscar',
                        filename=os.path.join(
                            dirname,
                            "refined-structure{}.vasp".format(
                                unique_index.format(k))))

        num += len(enumerated_structures)
        print(
            'The enumeration found {}({}+{}) distinct structures for {} '
            'with {} composition.'.format(
                len(enumerated_structures), prev,
                len(supplemental_structures), species,
                composition_list))
    print(f'{num} distinct structures are found totally.')


def run(json_file_path, max_structures, to_vasp):
    """
    Function to read the json file and load the parameters.

    Args:
        json_file_path: JSON file path direct to the file which includes
            all the important parameters.
        max_structures: Number of structures to be returned at most.
        to_vasp: Whether to generate all output slab models.

    """

    # Load json input file
    with open(json_file_path) as fp:
        data = json.load(fp)

    # Generate all composition combo
    to_be_enumerated_species = list(data['replacements'])
    num_layers_relaxed = data['num_layers_relaxed']
    composition_list = list(
        get_values_nested_dict(data['replacements']))
    composition_list = [sorted(x, reverse=True) for x in
                        composition_list]

    print("target_cell_size = {}".format(data['max_cell_size']))
    print("Composition of {} on the surface will be {}, "
          "respectively.".format(to_be_enumerated_species, composition_list))

    composition_combo = []
    for combo in product(*composition_list):
        replace_dict = {}
        for key, value in zip(to_be_enumerated_species, combo):
            replace_dict[key] = {key: float(value)}
        composition_combo.append(replace_dict)

    automate_surface(target_slab_path=data['target_slab_path'],
                     species=to_be_enumerated_species,
                     replace=composition_combo,
                     num_layers_relaxed=num_layers_relaxed,
                     max_cell_size=data['max_cell_size'],
                     symmetric=data['symmetric'],
                     max_structures=max_structures,
                     to_vasp=to_vasp)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__ + "\n{} {}".format(__date__, __author__),
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(
        "json_file_path",
        help="Path to an input file in VASP's POSCAR format.")

    parser.add_argument(
        "--max-structures", "-max",
        help="Number of structures to be returned at most for each "
             "composition.",
        type=int,
        default=2000
    )

    parser.add_argument(
        "--generate-poscar", "-g",
        help="Generate POSCAR files of enumerated structures.",
        action="store_true")

    args = parser.parse_args()

    run(args.json_file_path,
        args.max_structures,
        args.generate_poscar)
