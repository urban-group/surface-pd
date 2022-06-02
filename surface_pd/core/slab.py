import copy
import numpy as np
import collections
from typing import Sequence, Union

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.surface import get_slab_regions
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.typing import ArrayLike, CompositionLike, SpeciesLike

from surface_pd.error import NoInversionSymmetryError


def group_atoms_by_layer(o_layer, diff=0.01, max_diff=0.03):
    """
    This function is used to group misclassified atoms into the right
    layers. For example, c_atom1 = 0.01, c_atoms2 = 0.02, they should be
    classified in one layer. But actually they are not. So this function
    here will search the difference between two closest atoms, if the
    difference is smaller than diff, they will be regrouped in the same
    layer.

    Args:
        o_layer: dictionary which contains different c fractional
          coordinates as keys and number of atoms as values.
        diff: Accepted c fractional coordinates difference between two
          atoms.
        max_diff: Maximum accepted c fractional coordinates difference
          between two atoms.

    Returns: dictionary which contains different c fractional
      coordinates as keys and number of atoms as values.

    """
    height_sorted = sorted(o_layer.keys(), reverse=True)
    res = dict()
    excluded_heights = set()
    for height in height_sorted:
        if height in excluded_heights:
            continue
        res[height] = o_layer[height]
        for i in range(1, int(max_diff / diff) + 1):
            curr_height = height - i * diff
            if curr_height not in o_layer:
                break
            res[height] += o_layer[curr_height]
            excluded_heights.add(curr_height)
    return res


class Slab(Structure):
    """

    """
    def __init__(self, lattice: Union[ArrayLike, Lattice],
                 species: Sequence[CompositionLike],
                 coords: Sequence[ArrayLike],
                 charge: float = None,
                 validate_proximity: bool = False,
                 to_unit_cell: bool = False,
                 coords_are_cartesian: bool = False,
                 site_properties: dict = None,
                 ):
        """
        Constructor of Slab calss.
        This child class is inherited from the pymatgen.core.structure.
        The args defined have the same meaning as in the pymatgen
        documentation.
        """
        super().__init__(lattice, species, coords,
                                   charge, validate_proximity,
                                   to_unit_cell, coords_are_cartesian,
                                   site_properties)

    def wrap_pbc(self,
                 slab_direction=2,
                 tolerance=0.01):
        """
        Wrap fractional coordinates back into the unit cell.

        Arguments:
          slab_direction:
          coo (ndarray): 3-vector of fractional coordinates
          tolerance (float): tolerance in multiples of the lattice vector applied
            when wrapping in direction of the slab, so that atoms at the bottom
            of the slab are not wrapped to the top of the unit cell
        Returns:
          0.0 <= coo[i] < 1.0  for the directions within the slab planes and
          -tolerance <= coo[i] < (1.0 - tolerance)  for the slab direction

        """
        EPS = np.finfo(np.double).eps
        for site in self:
            bounds = np.array([0.0, 0.0, 0.0])
            bounds[slab_direction] = tolerance
            for i in range(3):
                while site.frac_coords[i] < -bounds[i]:
                    site.frac_coords[i] += 1.0
                while site.frac_coords[i] >= (1.0 - bounds[i] - EPS):
                    site.frac_coords[i] -= 1.0
        return self

    def layer_classification(self):
        """
        This function is used to classify different layers in the slab model. For
        the (001) surface, it will be classified as Li-plane, TM-plane,
        and O-plane. For the (104) surface, only one Li-TM-O plane will exist.

        Returns:
            c-fractional coordinates of the upper and lower limit of the
                central fixed region
                c-fractional coordinates of the surface Li/TM and O atoms (for
                the (104) surface, these two values can be the same)
                A dictionary which c-fractional coordinates of Li atoms as key
                and number of Li atoms at this layer as value
                A list which has the c-fractional coordinates of each layer

        """

        # Initialize three dictionaries to store layer information of Li,
        # TM, and O atoms
        Li_layers = collections.defaultdict(int)
        TM_layers = collections.defaultdict(int)
        O_layers = collections.defaultdict(int)

        for s in self:
            if 'Li' in s:
                # Since the c fractional coordinates of atoms even in the
                # exactly same layers are not identical, therefore, all c
                # fractional coordinates will be rounded to 2 decimal. Any
                # misclassified atoms will be regrouped by
                # "group_atoms_by_layer" function.
                Li_layers[round(s.frac_coords[2], 2)] += 1
            elif 'O' in s:
                O_layers[round(s.frac_coords[2], 2)] += 1
            else:
                TM_layers[round(s.frac_coords[2], 2)] += 1

        # Further group atoms by layer
        Li_layers = group_atoms_by_layer(Li_layers)
        TM_layers = group_atoms_by_layer(TM_layers)
        O_layers = group_atoms_by_layer(O_layers)

        layers = sorted({**Li_layers, **TM_layers})
        O_layers_c_frac = sorted(O_layers)

        # Collect the central layers based on the "selective_dynamics"
        center_sites = []
        for s in self:
            if ('selective_dynamics' in s.properties
                    and not any(s.properties['selective_dynamics'])):
                center_sites.append(copy.deepcopy(s))

        # Get the central slab boundaries
        center = Structure.from_sites(center_sites)
        center_region = get_slab_regions(center)
        lower_limit, upper_limit = center_region[0]

        # Li-terminated slab model -- (001) surface
        if len(Li_layers) > len(TM_layers):
            c_max_surface_Li = sorted(Li_layers.items())[-1][0]
            c_max_surface_O = sorted(O_layers.items())[-1][0]
            return [lower_limit, upper_limit, c_max_surface_Li,
                    c_max_surface_O, Li_layers, layers,
                    O_layers_c_frac, TM_layers]
        # TM terminated slab model -- (001) surface
        elif len(Li_layers) < len(TM_layers):
            c_surface_TM = sorted(TM_layers.items())[-1][0]
            c_max_surface_O = sorted(O_layers.items())[-1][0]
            return [lower_limit, upper_limit, c_surface_TM,
                    c_max_surface_O, Li_layers, layers,
                    O_layers_c_frac, TM_layers]
        # Li-TM-O terminated slab model -- (104) surface
        else:
            c_surface_TM_Li = sorted(Li_layers.items())[-1][0]
            c_max_surface_O = sorted(O_layers.items())[-1][0]
            return [lower_limit, upper_limit, c_surface_TM_Li,
                    c_max_surface_O, Li_layers, layers,
                    O_layers_c_frac, TM_layers]

    def index_extraction(self,
                         tol=0.01):
        """
        This function is used to get the index of the first and second layers
        which containing Li atoms, surface oxygen atoms as well as the relaxed
        Li atoms in the slab model
        Args:
            tol: Maximum spread (in fractional coordinates) that atoms in the
            same plane may exhibit (default: 0.01)

        Returns: indices of the first and second layers which containing Li
        atoms, surface oxygen atoms as well as the relaxed Li atoms in the slab
        model

        """
        [_, _, _, c_max_surface_O,
         _, _, O_layers_c_frac, _] = self.layer_classification()

        # Initialize two lists to store info of the first and second layers
        # which containing Li atoms (only be used for the (104) surface)
        # Li_layers_frac_coords = list(Li_layers.keys())
        # first_layer, second_layer = [Li_layers_frac_coords[0],
        #                              Li_layers_frac_coords[-1]], \
        #                             [Li_layers_frac_coords[1],
        #                              Li_layers_frac_coords[-2]]

        # Get indices of the first and second Li layers as well as the surface
        # oxygen layers
        relaxed_li_index, oxygen_index = [], []
        for index, site in enumerate(self):
            if 'Li' in site:
                if site.properties["selective_dynamics"] == [True, True, True]:
                    relaxed_li_index.append(index)
                # # First layer Li atoms index extraction
                # if first_layer[0] - tol <= site.frac_coords[2] <= \
                #         first_layer[0] + tol:
                #     first_layer_index.append(index)
                # if first_layer[1] - tol <= site.frac_coords[2] <= \
                #         first_layer[1] + tol:
                #     first_layer_index.append(index)
                # # Second layer Li atoms index extraction
                # if second_layer[0] - tol <= site.frac_coords[2] <= \
                #         second_layer[0] + tol:
                #     second_layer_index.append(index)
                # if second_layer[1] - tol <= site.frac_coords[2] <= \
                #         second_layer[1] + tol:
                #     second_layer_index.append(index)
            # Surface O atoms index extraction (both top and bottom surface O
            # atoms are all included)
            elif 'O' in site:
                if (site.frac_coords[2] >= c_max_surface_O - tol) or\
                        (site.frac_coords[2] <= O_layers_c_frac[0] + tol):
                    oxygen_index.append(index)

        return relaxed_li_index, oxygen_index

    def remove_sites_with_scaling(self,
                                  index,
                                  scaling=True,
                                  scaling_matrix=None):
        """
        This function is used to remove atoms in the slab model first and
        then make a supercell based on the scaling matrix.
        Args:
            index: indices of the atoms that are going to be removed
            scaling: whether to create a supercell of the slab model with
            vacancies created (default: True)
            scaling_matrix: scaling matrix to be used to create the supercell

        Returns: scaled slab model with vacancies created on the surface

        """
        if scaling_matrix is None:
            scaling_matrix = [2, 1, 1]
        copy_structure = self.copy()
        copy_structure.remove_sites(index)
        if scaling:
            copy_structure.make_supercell(scaling_matrix=scaling_matrix)
        return copy_structure

    def symmetrize_top_base(self,
                            symprec=1e-4,
                            direction=2,
                            tol=0.01):
        """
        This function is used to symmetrize the enumerated slab models using
        top surface as base. The inversion symmetry center is determined based
        on the central fixed region only because the input slab of this function
        is not symmetry since one the one hand, the top surface of the slab is
        substituted by other dummy species, on the other hand, it has been
        enumerated with different Li and O compositions on the surface.
        Therefore using the whole slab to detect the inversion symmetry does not
        work.
        Args:
            target_slab: enumerated slab model with Li and oxygen vacancies on
            the top surface
            symprec: tolerance for the symmetry detection (default: 0.1).
            direction: lattice direction perpendicular to the surface (default: 2)
            tol:

        Returns: symmetrized slab models with top and bottom surfaces are
        equivalent

        """

        # Load slab which is going to be symmetrized
        slab_tgt = self

        # Use center region to get the inversion symmetry center and the rotation
        # matrix
        lower_limit, upper_limit, _, _, _, _, _, _ = self.layer_classification()

        # Generate the reference slab which is just the central fixed region
        slab_ref_site = []
        for s in slab_tgt:
            if lower_limit - tol < s.frac_coords[2] < upper_limit + tol:
                slab_ref_site.append(copy.deepcopy(s))
        slab_ref = Structure.from_sites(slab_ref_site)

        # Determine symmetry operations of the reference slab and make sure
        # the reference slab has an inversion center
        sga = SpacegroupAnalyzer(slab_ref, symprec=symprec)
        if not sga.is_laue():  # has Laue symmetry (centro-symmetry)
            raise NoInversionSymmetryError

        # get the inversion operation and the inversion center (origin)
        ops = sga.get_symmetry_operations()
        inversion = ops[1]
        # print("trans_vec", inversion.translation_vector)
        assert (np.all(inversion.rotation_matrix == -np.identity(3)))
        origin = inversion.translation_vector / 2

        # wrap target slab models to unit cell
        # for s in slab_tgt:
        #     s.frac_coords = wrap_pbc(s.frac_coords, slab_direction=direction)
        # New way to implememt wrap_pbc
        self = self.wrap_pbc()

        # Symmetrized slab models based on top only
        top_sites = []
        bottom_sites = []
        center_sites = []
        for s in [s for s in self
                  if s.frac_coords[direction] > origin[direction]]:
            if ('selective_dynamics' in s.properties
                    and not any(s.properties['selective_dynamics'])):
                center_sites.append(copy.deepcopy(s))
            else:
                top_sites.append(copy.deepcopy(s))
        for s in [s for s in self
                  if s.frac_coords[direction] < origin[direction]]:
            if ('selective_dynamics' in s.properties
                    and not any(s.properties['selective_dynamics'])):
                center_sites.append(copy.deepcopy(s))
            else:
                bottom_sites.append(copy.deepcopy(s))
        # Extend the center sites list in case of some sites have exactly the
        # same c fractional coordinates of the origin
        center_sites.extend(
            [copy.deepcopy(s) for s in self if s.frac_coords[direction] ==
             origin[direction]])

        # Initialize a new slab model which only has the center and top regions of
        # the initial slab model
        sites = center_sites[:] + top_sites[:]
        for s in top_sites:
            s2 = copy.deepcopy(s)
            s2.frac_coords = inversion.operate(s.frac_coords)
            # Add symmetrized top sites into the new slab as the bottom sites
            sites.append(s2)

        symmetrized_slab_top = Slab.from_sites(sites)

        # New way to implememt wrap_pbc
        # for s in symmetrized_slab_top:
        #     s.frac_coords = wrap_pbc(s.frac_coords, slab_direction=direction)
        symmetrized_slab_top = symmetrized_slab_top.wrap_pbc()
        symmetrized_slab_top = symmetrized_slab_top.get_sorted_structure()

        return symmetrized_slab_top

    def surface_substitute(self, subs1: str, subs2: str,
                           direction=2, tol=0.02):
        """
        This function is used to substitute surface Li and O atoms with "dummy
        species" which will facilitate the enumeration code to detect target atoms.
        Args:
            subs1: Substitution atom for Li atom
            subs2: Substitution atom for O atom
            direction: Define the direction of substitution (c direction is the
            general case) (default: 2)
            tol: Maximum spread (in fractional coordinates) that atoms in the
            same plane may exhibit (default: 0.02)

        Returns: Substituted slab model

        """
        # Load slab
        slab_tgt = self

        # Define c-fractional coordinates of upper boundary of the fixed
        # region, surface metal, and surface O atoms.
        [_, center_top, c_max_surface_metal,
         c_max_surface_O, _, _, _, _] = self.layer_classification()

        # This distance is calculated by using the surface metal (can be Li or
        # TM  atoms) minus the top boundary of central fixed region,
        # which is the thickness of the top relaxed region.
        distance = (c_max_surface_metal - center_top)

        # Determine the surface Li and O atoms
        surface_Li = []
        surface_O = []
        for s in slab_tgt:
            if (abs(c_max_surface_metal - s.frac_coords[direction]) <
                    distance - tol and ('Li' in s)):
                surface_Li.append(copy.deepcopy(s))
            if (abs(c_max_surface_O - s.frac_coords[direction]) <
                    tol and ('O' in s)):
                surface_O.append(copy.deepcopy(s))

        # Extract indices for surface Li and O atoms in target slab model
        def extract_indices(tgt_idx, ref_idx):
            idx = []
            for i, s in enumerate(tgt_idx):
                for j, t in enumerate(ref_idx):
                    if t == s:
                        idx.append(i)
            return idx

        idx_surface_Li = extract_indices(slab_tgt, surface_Li)
        idx_surface_O = extract_indices(slab_tgt, surface_O)

        # Substitute surface Li and O atoms with dummy species
        slab_surface_substitute = slab_tgt.copy()
        for i in range(len(slab_surface_substitute)):
            if i in idx_surface_Li:
                slab_surface_substitute.replace(
                    i, subs1,
                    properties={'selective_dynamics': [True, True, True]}
                )
            if i in idx_surface_O:
                slab_surface_substitute.replace(
                    i, subs2,
                    properties={'selective_dynamics': [True, True, True]}
                )
        return slab_surface_substitute

    def get_max_min_c_frac(self):
        """
        Get the maximum and minimum values of the c fractional coordinates for
        all sites. This is used to make sure whether the slab region is in the
        middle or at the bottom of the unit cell
        Args:
            structure: enumerated, after refined slab model

        Returns:
            maximum and minimum values of the c fractional coordinates
        """
        min, max = 1, 0
        for site in self:
            if site.frac_coords[2] > max:
                max = site.frac_coords[2]
            if site.frac_coords[2] < min:
                min = site.frac_coords[2]
        return min, max

    def Li_TM_layers_finder(self):
        """
        This function is used to determine the c fractional coordinates of
        central region of the slab models as well as the surface metal and O
        atoms.

        Args:
            input_structure: any slab models

        Returns: upper limit and lower limit of the central region of the
          slab, surface metal and O c fractional coordinates

        """
        Li_layers = collections.defaultdict(int)
        TM_layers = collections.defaultdict(int)
        O_layers = collections.defaultdict(int)
        for s in self:
            if 'Li' in s:
                # Since the c fractional coordinates of atoms even in the
                # exactly same layers are not identical, therefore, all c
                # fractional coordinates will be rounded to 2 decimal. Any
                # misclassified atoms will be regrouped by
                # "group_atoms_by_layer" function.
                Li_layers[round(s.frac_coords[2], 2)] += 1
            elif 'O' in s:
                O_layers[round(s.frac_coords[2], 2)] += 1
            else:
                TM_layers[round(s.frac_coords[2], 2)] += 1

        Li_layers = group_atoms_by_layer(Li_layers)
        TM_layers = group_atoms_by_layer(TM_layers)
        O_layers = group_atoms_by_layer(O_layers)
        return Li_layers, TM_layers, O_layers