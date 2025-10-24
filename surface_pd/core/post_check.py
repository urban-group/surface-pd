"""
Post-processing validation module for enumerated structures.

This module provides the PostCheck class for validating enumerated slab
structures after enumeration, ensuring they maintain required symmetry
properties and geometric consistency.
"""

import numpy as np
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from surface_pd.core.slab import Slab
from surface_pd.error import (
    IncompatibleSymmError,
    NoInversionSymmetryError,
    NonCentralInversionSymmetryError,
    PrimitiveStructureFinderError,
)


class PostCheck:
    def __init__(
        self,
        structure: Slab,
        tolerance: float = 0.03,
        direction: int = 2,
        symprec: float = 1e-2,
    ):
        self.structure = structure
        self.tolerance = tolerance
        self.direction = direction
        self.symprec = symprec

    def has_correct_num_sites(self, curr_num_sites: int):
        """
        Check if the enumerated structure has the correct number of sites.

        Args:
            curr_num_sites: Expected number of sites.

        Returns
        -------
            True if number of sites matches expectation.
        """
        return self.structure.num_sites == curr_num_sites

    def post_check(
        self,
        species,
        composition_list,
        index,
        keep_symmetric,
        criteria,
        symprec=1e-5,
    ):
        r"""
        Perform final check to see whether the enumerated slab models are
        correct.

        1. Has the correct geometry.
        2. Has the inversion symmetry center even with a quite high symmetry
           detection parameter.
        3. Has the inversion symmetry center shifted to the origin (0, 0, 0).

        Args:
            species: Target species that will be enumerated.
            composition_list: Compositions of enumerated species.
            index: Unique index of this slab model.
            keep_symmetric: Whether the symmetry of the slab model will be
                kept after the enumeration.
            criteria:
            symprec: Tolerance for symmetry finding. Defaults to 1e-5.

        Returns
        -------
            If everything is fine, it will just pass. If something
            unexpected happens, the specific error report will generate.

        """
        if keep_symmetric:
            try:
                symmetric, origin, _ = self.structure.is_symmetry(
                    symprec=symprec, return_isc=True
                )
            except TypeError:
                symmetric = self.structure.is_symmetry(
                    symprec, return_isc=False
                )

            if not symmetric:
                species_str = ""
                if isinstance(species, list):
                    for spec in species:
                        species_str += f"{spec}"
                else:
                    species_str = f"{species}"
                composition_str = (
                    f"{composition_list}".format()
                    .replace("[", "")
                    .replace("]", "")
                    .replace(", ", "-")
                )

                # Save structure for debugging
                self.structure.to(
                    fmt="poscar",
                    filename=(f"debug-enumed-{species_str}_"
                              + f"{composition_str}-{index}.vasp"),
                )
                print(
                    f'Please see the saved "debug-enumed-{species_str}'
                    f'_{composition_str}_{index}.vasp" '
                    "structure and see if it makes sense."
                )
                raise NoInversionSymmetryError

            # Make sure that the inversion symmetry center is in the c
            # direction origin
            if (
                abs(origin[self.direction])
                > self.tolerance + self.symprec * 100
            ):
                print(origin)
                raise NonCentralInversionSymmetryError
        return self.structure

    def get_refined_structure(
        self,
        species,
        composition_list,
        index,
        symprec=1e-1,
    ):
        """
        Get refined primitive structure.

        Args:
            species: Target species that will be enumerated.
            composition_list: Compositions of enumerated species.
            index: Unique index of this slab model.
            symprec: Tolerance for symmetry finding. Defaults to 1e-1.

        Returns
        -------
            Refined primitive structure.
        """
        refined_structure = Slab.from_sites(
            SpacegroupAnalyzer(
                self.structure, symprec=symprec
            ).get_refined_structure()
        )
        refined_structure.direction = self.direction
        refined_structure.tolerance = self.tolerance
        refined_structure.to_be_enumerated_species = (
            self.structure.to_be_enumerated_species
        )
        refined_structure.num_layers_enumed = self.structure.num_layers_enumed
        refined_structure.symmetric = self.structure.symmetric

        if (
            refined_structure.num_sites % self.structure.num_sites == 0
            and refined_structure.num_sites != self.structure.num_sites
        ):
            refined_structure.make_supercell(
                np.array(
                    [
                        [
                            self.structure.num_sites
                            / refined_structure.num_sites,
                            0,
                            0,
                        ],
                        [
                            0,
                            self.structure.num_sites
                            / refined_structure.num_sites,
                            0,
                        ],
                        [0, 0, 1],
                    ]
                )
            )

        if self.structure.num_sites != refined_structure.num_sites:
            species_str = ""
            if isinstance(species, list):
                for spec in species:
                    species_str += f"{spec}"
            else:
                species_str = f"{species}"
            composition_str = (
                f"{composition_list}".format()
                .replace("[", "")
                .replace("]", "")
                .replace(", ", "-")
            )

            # Save structure for debugging
            refined_structure.to(
                fmt="poscar",
                filename=(f"debug-refined-{species_str}"
                          + f"_{composition_str}-{index}.vasp"),
            )
            print(
                f'Please see the saved "debug-refined-{species_str}'
                f'_{composition_str}_{index}.vasp" '
                "structure and see if it makes sense."
            )
            raise PrimitiveStructureFinderError

        return refined_structure

    def user_symmetric_compatible(
        self, user_symmetric: bool, code_determined_symmetric: bool
    ):
        """
        Check if the user-defined symmetric parameter is compatible.

        Args:
            user_symmetric: User-defined symmetric parameter.
            code_determined_symmetric: Code-determined symmetric parameter.

        Returns
        -------
            True if compatible.

        Raises
        ------
            IncompatibleSymmError: If parameters are incompatible.
        """
        if user_symmetric != code_determined_symmetric:
            raise IncompatibleSymmError
        else:
            return True
