"""Pre-enumeration validation helpers for slab models."""

import logging
from shutil import which

from monty.dev import requires
from pymatgen.core.surface import get_slab_regions

from surface_pd.core.enumeration_slab import EnumerationSlab
from surface_pd.util.util import check_int

logger = logging.getLogger(__name__)


class PreCheck:
    """
    PreCheck class to check whether the input slab model is valid.

    Args:
        structure: Input slab model.

    """

    def __init__(self, structure: EnumerationSlab, analysis=None):
        self.structure = structure
        self.analysis = analysis

    enum_cmd = which('enum.x')
    makestr_cmd = which('makestr.x') or which("makeStr.x") or which(
        "makeStr.py")

    @requires(
        enum_cmd and makestr_cmd,
        "In order to run this script, you need to have the executable "
        "'enum.x' "
        "or 'makestr.x' files and add them in the path. For the detailed "
        "compilation, please refer to the documentation of this package ("
        "Installation)."
    )
    def is_slab(self):
        """Check whether the input structure is a slab."""
        try:
            ranges = get_slab_regions(self.structure)
        except ValueError:
            return False
        else:
            if len(ranges) == 2:
                slab_width = ranges[0][1] - ranges[1][0]
            else:
                slab_width = ranges[0][1] - ranges[0][0]
            if slab_width < 0.9:
                return True
            else:
                return False

    def is_cuboid(self):
        """Check whether the input structure is a cuboid."""
        if max(self.structure.lattice.abc) != self.structure.lattice.c:
            return False
        else:
            return True

    def all_has_selective_dynamics(self):
        """Check if all sites have selective dynamics properties."""
        for site in self.structure:
            try:
                sd = site.properties["selective_dynamics"]
                # Handle both list and numpy array
                if (all(not x for x in sd) or all(x for x in sd)):
                    return True
            except KeyError:
                return False
            else:
                sd = site.properties['selective_dynamics']
                if len(sd) == 1 and not sd[0]:
                    return False

    def has_inversion_symmetry(self, symprec=0.1):
        """
        Check if the input structure has inversion symmetry.

        Args:
            symprec: Symmetry detection tolerance. Defaults to 0.1. (
                fairly strict).

        """
        return EnumerationSlab.from_sites(self.structure).is_symmetry(
            symprec=symprec, return_isc=False
        )

    def relax_both_surfaces(self):
        """Check whether only one side of the surface is relaxed."""
        analysis = self.analysis or self.structure.analyze()
        fixed = set(analysis.fixed_site_indices)
        bottom_surface = set(analysis.layers[0].site_indices)
        if fixed & bottom_surface:
            logger.info(
                "The slab model provided has the whole bottom surface fixed; "
                "no symmetrization is needed."
            )
            return False
        else:
            return True

    def has_validate_composition(self,
                                 replace,
                                 max_cell_size):
        """
        Check whether the user defined composition can be enumerated.

        Args:
            replace: Species and occupancy dictionaries containing the species
                mapping in string-string pairs. E.g. {'Li': {'Li': 0.5}},
                stored in the list.
            max_cell_size: Maximum number of supercells of the input slab.

        """
        is_int = []

        analysis = self.analysis or self.structure.analyze()
        selected_indices = analysis.enumerated_site_indices
        num_atoms = {
            key: len(value) for key, value in selected_indices.items()
        }
        for comp in replace:
            for s, n_s in num_atoms.items():
                if self.structure.symmetric:
                    f = comp[s][s] * n_s * max_cell_size / 2
                    f = check_int(f)
                else:
                    f = comp[s][s] * n_s * max_cell_size
                    f = check_int(f)
                if f.is_integer():
                    is_int.append(True)
                else:
                    is_int.append(False)
                    break
        if all(is_int):
            return True
        else:
            return False
