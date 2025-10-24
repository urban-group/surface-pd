#!/usr/bin/env python

"""
This code will enumerate the surface of the input slab model with user defined
target species and composition.
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
__date__ = "2023-08-23"
__version__ = "0.1.0"

import argparse
import copy
import json
import os
import warnings
from itertools import product

from pymatgen.core.periodic_table import DummySpecies
from pymatgen.core.surface import get_slab_regions
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from surface_pd.analysis.slab_analysis import (
    selective_dynamics_completion,
    structure_filter,
)
from surface_pd.core import EnumWithComposition, PostCheck, PreCheck, Slab
from surface_pd.error import (
    InvalidInputFormatError,
    IncompatibleSymmError,
    NonDefinedSelectiveDynamicsError,
    NonSlabError,
    InvalidCompositionError,
    SlabOrientationError,
    NoInversionSymmetryError
)
from surface_pd.util import (
    all_int,
    get_values_nested_dict,
    have_zero,
    replace_dummy,
)

warnings.filterwarnings("ignore")


def automate_surface(
    target_slab_path: str,
    species: list,
    replace: list,
    num_layers_enumed: dict,
    max_cell_size: int,
    symmetric: bool,
    to_vasp: bool = False,
):
    """
    This function contains the general framework to enumerate the parent
    slab model surface with different target species and compositions.
    TODO:
        1. Check the meaning of num_layers_enumed. Is this really the number of
        layers enumed or number of layers relaxed? Can the number of layers
        enumed be different than the number of layers relaxed?
        2. symmetric parameter: If turned off, is the only top surface
        enumerated?
        4. Check enumlib code surface enumeration and see how difficult it
        is to implement.
    Args:
        target_slab_path: Input slab model
        species: Target species on the surface that will be enumerated.
        replace: Species and occupancy dictionaries containing the species
            mapping in string-string pairs. E.g. {'Li': {'Li': 0.5}}
            (This dict means that only half of the original Li
            atoms in the slab model will be kept), stored in the list.
        num_layers_enumed: Dictionaries where species and number of relaxed
            layers as keys and values, respectively. E.g. {'Li': 2, 'O': 1},
            which represents 2 and 1 layers of Li and O will be relaxed
            during the DFT geometry optimization.
        max_cell_size: Maximum number of supercells of the input slab.
        symmetric: Whether the symmetry of the slab model will be kept
            after the enumeration.
        to_vasp: Whether to generate all the output slab models in VASP
        format. Defaults to False.

    Returns
    -------
        All enumerated slab models with different composition of target
        species.

    """
    # Load parent slab model
    input_structure = Slab.from_file(target_slab_path)

    # Wrap out of the boundary fractional coordinates back into the unit cell.
    input_structure.wrap_pbc()

    # Define important parameters which will be used later, including to be
    # enumerated species, number of layers relaxed, whether to keep the
    # symmetry, direction of the slab, and the proximity tolerance for
    # adjacent atoms
    input_structure.to_be_enumerated_species = species
    input_structure.num_layers_enumed = num_layers_enumed
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
    center_bottom, center_top, relaxed_index = (
        input_structure.index_extraction(only_top=False)
    )

    # Define the "dummy" species that will be used to substitute the target
    # species
    dummy_species = [
        DummySpecies(symbol="X" + str(i + 1)) for i in range(len(species))
    ]

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
                subs_dict, relaxed_index
            )
            supplemental_structures = [structure]
            enumerated_structures = [structure]
        else:
            if have_zero(composition_list):
                slab_substituted = (
                    slab_substituted.supplemental_structures_gene(
                        subs_dict, relaxed_index
                    )
                )

            subs_dict = replace_dummy(subs_dict, dummy_species)

            ewc = EnumWithComposition(
                subs_dict=subs_dict, max_cell_size=max_cell_size
            )

            structures = ewc.apply_enumeration(slab_substituted)

            # Filtered out the structures which has c lattice as the
            # largest lattice (a tall cuboid)
            filtered_structures = structure_filter(
                structures, direction=direction, criteria=criteria
            )

            # Add selective dynamics for enumerated sites
            for filtered_structure in filtered_structures:
                selective_dynamics_completion(
                    structure=filtered_structure,
                    direction=direction,
                    dummy_species=dummy_species,
                    center_bottom=center_bottom,
                    center_top=center_top,
                    tolerance=tolerance,
                )

            # Symmetrize slab models based on the top enumerated
            # surface if needed
            enumerated_structures = []
            for filtered_structure in filtered_structures:
                for i in range(len(species)):
                    filtered_structure.replace_species(
                        {dummy_species[i]: species[i]}
                    )
                filtered_structure = Slab.from_sites(filtered_structure)
                if symmetric:
                    enumerated_structures.append(
                        filtered_structure.symmetrize_top_base()
                    )
                else:
                    enumerated_structures.append(filtered_structure)

            prev = len(enumerated_structures)

        unique_index = "-{:0%dd}" % (len(str(len(enumerated_structures))))
        for k, symmetrized_structure in enumerate(enumerated_structures):
            indicator = 0
            if symmetric:
                _, origin, _ = symmetrized_structure.is_symmetry(
                    return_isc=True
                )
                min_c, _ = symmetrized_structure.get_max_min_c_frac()
                sga = SpacegroupAnalyzer(symmetrized_structure, symprec=0.1)
                # Create the refined structures
                refined_structure = sga.get_refined_structure()

                # Define number of sites that should be after enumeration
                total_num_sites = input_structure.calculate_num_sites(
                    composition_list=composition_list,
                    relaxed_index=relaxed_index,
                    max_cell_size=max_cell_size,
                )

                num_sites = refined_structure.num_sites

                while (
                    (num_sites > total_num_sites)
                    or (max(refined_structure.lattice.abc) > criteria * 2 - 5)
                    or (
                        max(refined_structure.lattice.abc)
                        != refined_structure.lattice.abc[direction]
                    )
                ):
                    pc = PostCheck(refined_structure)
                    indicator, refined_structure = pc.slab_size_check(
                        total_num_sites=total_num_sites,
                        enumerated_num_sites=num_sites,
                        criteria=criteria,
                    )
                    num_sites = refined_structure.num_sites

                refined_structure = Slab.from_sites(refined_structure)
                noncontiguous = len(get_slab_regions(refined_structure)) == 2
                # print(input_structure.layers_finder())

                if noncontiguous:
                    refined_structure = refined_structure.tune_isc(
                        origin=origin, shift_isc_back=True
                    )
                    refined_structure = refined_structure.tune_c(
                        target_min_c=min_c
                    )
                    # print(refined_structure.layers_finder())
                    refined_structure = (
                        refined_structure.add_selective_dynamics(
                            lower_limit=center_bottom, upper_limit=center_top
                        )
                    )
                else:
                    refined_structure = refined_structure.tune_c(
                        target_min_c=min_c
                    )
                    # print(refined_structure.layers_finder())
                    refined_structure = (
                        refined_structure.add_selective_dynamics(
                            lower_limit=center_bottom, upper_limit=center_top
                        )
                    )

                refined_structure = refined_structure.tune_isc(
                    origin=origin, shift_isc_back=False
                )
            else:
                refined_structure = symmetrized_structure

            # Perform final check
            pc = PostCheck(refined_structure)
            pc.post_check(
                species=species,
                composition_list=composition_list,
                keep_symmetric=symmetric,
                criteria=criteria,
                index=k,
            )

            # Generate slab models
            if to_vasp:
                # Make directories
                dirname = ""
                for n in range(len(species)):
                    dirname += str(composition_list[n]) + str(species[n])
                if not os.path.exists(dirname):
                    os.makedirs(dirname)

                if indicator == -1:
                    refined_structure.to(
                        fmt="poscar",
                        filename=os.path.join(
                            dirname,
                            f"refined-prim{unique_index.format(k)}.vasp",
                        ),
                    )
                elif indicator == 1:
                    refined_structure.to(
                        fmt="poscar",
                        filename=os.path.join(
                            dirname,
                            f"refined-reduced{unique_index.format(k)}.vasp",
                        ),
                    )
                elif indicator == 2:
                    refined_structure.to(
                        fmt="poscar",
                        filename=os.path.join(
                            dirname,
                            f"refined-rotated{unique_index.format(k)}.vasp",
                        ),
                    )
                else:
                    refined_structure.to(
                        fmt="poscar",
                        filename=os.path.join(
                            dirname,
                            f"refined-structure{unique_index.format(k)}.vasp",
                        ),
                    )

        num += len(enumerated_structures)
        print(
            f"The enumeration found {len(enumerated_structures)}"
            f"({prev}+{len(supplemental_structures)}) "
            f"distinct structures for {species} "
            f"with {composition_list} composition."
        )
    print(f"{num} distinct structures are found totally.")


def run(json_file_path, to_vasp):
    """
    Function to read the json file and load the parameters.

    Args:
        json_file_path: JSON file path direct to the file which includes
            all the important parameters.
        to_vasp: Whether to generate all output slab models.

    """
    # Load json input file
    try:
        with open(json_file_path) as fp:
            data = json.load(fp)
    except (json.decoder.JSONDecodeError, IsADirectoryError):
        raise InvalidInputFormatError

    # Generate all composition combo
    to_be_enumerated_species = list(data["replacements"])
    num_layers_enumed = data["num_layers_enumed"]
    composition_list = list(get_values_nested_dict(data["replacements"]))
    composition_list = [sorted(x, reverse=True) for x in composition_list]

    print("target_cell_size = {}".format(data["max_cell_size"]))
    print(
        f"Composition of {to_be_enumerated_species} on the surface "
        f"will be {composition_list}, respectively."
    )

    composition_combo = []
    for combo in product(*composition_list):
        replace_dict = {}
        for key, value in zip(to_be_enumerated_species, combo):
            replace_dict[key] = {key: float(value)}
        composition_combo.append(replace_dict)

    automate_surface(
        target_slab_path=data["target_slab_path"],
        species=to_be_enumerated_species,
        replace=composition_combo,
        num_layers_enumed=num_layers_enumed,
        max_cell_size=data["max_cell_size"],
        symmetric=data["symmetric"],
        to_vasp=to_vasp,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Enumerate surface structures with different compositions.\n\n"
            "This tool systematically generates all symmetrically unique "
            "surface structures by creating vacancies or substitutions on "
            "specified surface layers.\n\n"
            "Example usage:\n"
            "  %(prog)s input.json\n"
            "  %(prog)s input.json --generate-poscar\n"
            "  %(prog)s examples/enumeration-examples/input/input-Li.json -g"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            f"Version: {__version__}\n"
            f"Author: {__author__} ({__email__})\n"
            f"Date: {__date__}\n\n"
            "For detailed algorithm, see:\n"
            "  Morgan et al., Comput. Mater. Sci. 2017, 136, 144-149.\n"
            "  Hart et al., Comput. Mater. Sci. 2012, 59, 101-107."
        ),
    )

    parser.add_argument(
        "json_file_path",
        metavar="INPUT_JSON",
        help=(
            "Path to JSON configuration file containing:\n"
            "  - target_slab_path: VASP structure file\n"
            "  - replacements: species compositions to enumerate\n"
            "  - num_layers_enumed: layers to modify per species\n"
            "  - max_cell_size: maximum supercell expansion\n"
            "  - symmetric: maintain top-bottom symmetry"
        ),
    )

    parser.add_argument(
        "--generate-poscar",
        "-g",
        dest="generate_poscar",
        action="store_true",
        help=(
            "Generate VASP POSCAR files for all enumerated structures. "
            "Output organized in directories by composition (e.g., 1.0Li/)."
        ),
    )

    args = parser.parse_args()

    run(args.json_file_path, args.generate_poscar)
