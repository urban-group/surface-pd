"""
Post-processing validation module for enumerated structures.

This module provides the PostCheck class for validating enumerated slab
structures after enumeration, ensuring they maintain required symmetry
properties and geometric consistency.
"""

import copy
import logging

import numpy as np
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from surface_pd.core.enumeration_slab import EnumerationSlab
from surface_pd.error import (
    IncompatibleSymmError,
    NoInversionSymmetryError,
    NonCentralInversionSymmetryError,
    PrimitiveStructureFinderError,
)

logger = logging.getLogger(__name__)


class PostCheck:
    """Validate enumerated slab structures after generation."""

    def __init__(
        self,
        structure: EnumerationSlab,
        inversion_center_tolerance_fractional: float = 0.03,
        direction: int = 2,
        symprec: float = 1e-2,
    ):
        self.structure = structure
        self.inversion_center_tolerance_fractional = (
            inversion_center_tolerance_fractional
        )
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

    def repair_refined_slab_geometry(
        self,
        total_num_sites: int,
        enumerated_num_sites: int,
        criteria: float,
    ):
        """
        Repair refined slab geometry after symmetry refinement.

        Symmetric enumeration refines enumerated structures with pymatgen's
        space-group analyzer. The refined cell can be primitive, expanded, or
        oriented with the longest lattice vector away from the slab direction.
        This method applies the historical post-refinement adjustments used by
        the enumeration workflow so the returned slab has the expected size and
        orientation.

        Args:
            total_num_sites: Expected number of sites after enumeration.
            enumerated_num_sites: Actual number of sites in the refined slab.
            criteria: Lattice parameter perpendicular to the parent slab
                surface.

        Returns
        -------
            Tuple of ``(indicator, structure)``. The indicator preserves the
            historical filename behavior: ``-1`` means primitive/reduced,
            ``1`` means reduced, ``2`` means rotated, and ``0`` means no
            repair was needed.

        Raises
        ------
            PrimitiveStructureFinderError: If a refined slab has fewer sites
                than expected and cannot tile the expected site count.
        """
        if enumerated_num_sites > total_num_sites:
            refined_prim = self.structure.copy()
            refined_prim = (
                refined_prim.get_primitive_structure()
                .get_reduced_structure()
                .get_sorted_structure()
            )
            return -1, refined_prim

        if enumerated_num_sites < total_num_sites:
            if total_num_sites % enumerated_num_sites != 0:
                raise PrimitiveStructureFinderError
            return 1, self.structure

        if max(self.structure.lattice.abc) > criteria * 2 - 5:
            refined_prim = copy.deepcopy(self.structure)
            refined_prim = (
                refined_prim.get_primitive_structure()
                .get_reduced_structure()
                .get_sorted_structure()
            )
            return -1, refined_prim

        if max(self.structure.lattice.abc) != self.structure.lattice.c:
            refined_rotated = copy.deepcopy(self.structure)
            if max(self.structure.lattice.abc) == self.structure.lattice.a:
                refined_rotated.make_supercell(
                    [[0, 0, 1], [0, 1, 0], [1, 0, 0]]
                )
            elif max(self.structure.lattice.abc) == self.structure.lattice.b:
                refined_rotated.make_supercell(
                    [[1, 0, 0], [0, 0, 1], [0, 1, 0]]
                )
            return 2, refined_rotated

        return 0, self.structure

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
        Perform final checks on enumerated slab models.

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
                debug_filename = (
                    f"debug-enumed-{species_str}_"
                    f"{composition_str}-{index}.vasp"
                )
                self.structure.to(
                    fmt="poscar",
                    filename=debug_filename,
                )
                logger.error(
                    "Saved debug structure %s for enumerated slab without "
                    "inversion symmetry.",
                    debug_filename,
                )
                raise NoInversionSymmetryError

            # Make sure that the inversion symmetry center is in the c
            # direction origin
            if (
                abs(origin[self.direction])
                > self.inversion_center_tolerance_fractional
                + self.symprec * 100
            ):
                logger.error(
                    "Detected non-central inversion symmetry origin: %s",
                    origin,
                )
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
        refined_structure = EnumerationSlab.from_sites(
            SpacegroupAnalyzer(
                self.structure, symprec=symprec
            ).get_refined_structure()
        )
        refined_structure.direction = self.direction
        refined_structure.layer_tolerance_angstrom = (
            self.structure.layer_tolerance_angstrom
        )
        refined_structure.enumerated_species = (
            self.structure.enumerated_species
        )
        refined_structure.num_enumerated_layers = (
            self.structure.num_enumerated_layers
        )
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
            debug_filename = (
                f"debug-refined-{species_str}_"
                f"{composition_str}-{index}.vasp"
            )
            refined_structure.to(
                fmt="poscar",
                filename=debug_filename,
            )
            logger.error(
                "Saved debug structure %s for refined slab with unexpected "
                "site count.",
                debug_filename,
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
