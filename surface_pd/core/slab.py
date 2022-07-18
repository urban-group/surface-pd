import copy
import numpy as np
import collections
from typing import Union, List, Sequence

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.surface import get_slab_regions
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core import Element, Species, DummySpecies, Composition

from surface_pd.error import NoInversionSymmetryError


class Slab(Structure):
    """
    Constructor of periodic Slab class.
    This child class is inherited from the pymatgen.core.structure.Structure.
    The args defined have the same meaning as in the pymatgen
    documentation.

    Args:
        lattice (Union): Either a pymatgen.core.lattice.Lattice or any 2D
            array.
        species (Sequence): List of species on each site.
        coords (Sequence): List of fractional or cartesian coordinates of each
            species.
        charge (float): Overall charge of the structure.
        validate_proximity (bool): Whether to check if the distance between
            two sites is two close, i.e. less than 0.01 Ang. Defaults to False.
        to_unit_cell (bool): Whether to wrap fractional coordinates back
            into the unit cell. Defaults to False.
        coords_are_cartesian (bool): Whether the coordinates of sites are
            provided in cartesian coordinates. Defaults to False.
        site_properties (dict): Properties of each site in the slab model as a
            dict of sequences.

    """
    def __init__(self, lattice: Union[List, np.ndarray, Lattice],
                 species: Sequence[Union[str, Element, Species,
                                         DummySpecies, Composition]],
                 coords: Sequence[Sequence[float]],
                 charge: float = None,
                 validate_proximity: bool = False,
                 to_unit_cell: bool = False,
                 coords_are_cartesian: bool = False,
                 site_properties: dict = None,
                 ):
        super().__init__(lattice, species, coords,
                         charge, validate_proximity,
                         to_unit_cell, coords_are_cartesian,
                         site_properties)
        
        self.slab_direction = slab_direction  # Questions: How about defining this instance variable? because the slab_direction is used in many instance functions

    @staticmethod
    def group_atoms_by_layer(o_layer,
                             diff=0.01,
                             max_diff=0.03):
        """
        Group misclassified atoms into the right layers.
        For example, c_atom1 = 0.01, c_atoms2 = 0.02, they should be
        classified in one layer. But actually they are not. So this function
        here will search the difference between two closest atoms, if the
        difference is smaller than diff, they will be regrouped in the same
        layer.

        Args:
            o_layer: Dict which contains different c fractional
              coordinates as keys and number of atoms as values.
            diff: Accepted c fractional coordinates difference between two
              atoms. Defaults to 0.01.
            max_diff: Maximum accepted c fractional coordinates difference
              between two atoms. Defaults to 0.03.

        Returns:
            Dict which contains different c fractional
            coordinates as keys and number of atoms as values.

        """
#         height_sorted = sorted(o_layer.keys(), reverse=True)
#         res = dict()
#         excluded_heights = set()
#         for height in height_sorted:
#             if height in excluded_heights:
#                 continue
#             res[height] = o_layer[height]
#             for i in range(1, int(max_diff / diff) + 1):
#                 curr_height = height - i * diff  # Questions: Why do you only check it in units of diff (0.01)? The diff value may have to be consistent with the value used in proximity check. 
#                 if curr_height not in o_layer:
#                     break
#                 res[height] += o_layer[curr_height]
#                 excluded_heights.add(curr_height)

        # Suggestions (check and reflect it if this suggestion looks better):
        res = dict()
        o_layer_sorted = dict(sorted(o_layer.items(), key=lambda x: x[0], reverse=True))
        for i, (height_, num_) in enumerate(o_layer_sorted.items()):
            if i == 0:
                res[height_] = num_
                previous_height = height_
                continue

            if abs(previous_height - height_) <= max_diff:
                res[previous_height] += num_
            else:
                res[height_] = num_
                previous_height = height_
            
        return res

    def wrap_pbc(self,
                 slab_direction=2,
                 tolerance=0.01):
        """
        Wrap fractional coordinates back into the unit cell.

        Args:
          slab_direction: Lattice direction perpendicular to the surface,
            i.e. parallel to the c lattice parameter. Defaults to 2.
          tolerance: Tolerance in multiples of the lattice vector
            applied when wrapping in direction of the slab, so that atoms at
            the bottom of the slab are not wrapped to the top of the unit
            cell. Defaults to 0.01.

        Returns:
            Slab model with all sites in the unit cell

        """
        EPS = np.finfo(np.double).eps
        for site in self:
            bounds = np.array([0.0, 0.0, 0.0])
            bounds[slab_direction] = tolerance
            for i in range(3):
                while site.frac_coords[i] < -bounds[i]:  # Questions: I cannot understand what this part is for. I think it is enough to take the positive decimal fraction making use of the "np.modf" function 
                    site.frac_coords[i] += 1.0
                while site.frac_coords[i] >= (1.0 - bounds[i] - EPS):
                    site.frac_coords[i] -= 1.0
        return self

    def layer_distinguisher(self, precision=2):
        """
        This function is used to distinguish different layers in the slab
        model.
        For the (001) surface, it will be distinguished as Li-plane,
        TM-plane, and O-plane.
        For the (104) surface, only the Li-TM-O plane will exist.

        Args:
            precision: Round c fraction coordinates to a given precision in
                decimal digits. Defaults to 2.

        Returns:
            c fractional coordinates of the upper and lower boundaries
            of the central fixed region. (1, 2)\n
            c fractional coordinates of the surface Li/TM and O atoms
            (for the (104) surface, these two values can be the same).
            (3, 4) \n
            Dict where c fractional coordinates as key
            and number of atoms at this layer as values (Li layers,
            O layers, and TM layers, respectively). (5, 7, and 8) \n
            List which has the c-fractional coordinates of each layer. (6)

        """

        # Initialize three dictionaries to store layer information of Li,
        # TM, and O atoms
        Li_layers = collections.defaultdict(int)
        TM_layers = collections.defaultdict(int)
        O_layers = collections.defaultdict(int)

        # Questions 1: The following lines are only applicable to Li-TM-O assuming that the Li and TM are mixed. Do you think it is enough to only consider the ternary Li-TM-O?
        # Questions 2: The following lines are assuming that the surface normal direction is parallel to the c lattice. However, the possibility that the surface normal direction is along a or b lattice should also be considered using the "slab_direction" argument as in self.wrap_pbc
        for s in self:
            if 'Li' in s:
                # Since the c fractional coordinates of atoms even in the
                # exactly same layers are not identical, therefore, all c
                # fractional coordinates will be rounded to 2 decimal. Any
                # misclassified atoms will be regrouped by
                # "group_atoms_by_layer" function.
                Li_layers[round(s.frac_coords[2], precision)] += 1
            elif 'O' in s:
                O_layers[round(s.frac_coords[2], precision)] += 1
            else:
                TM_layers[round(s.frac_coords[2], precision)] += 1

        # Further group atoms by layer
        Li_layers = self.group_atoms_by_layer(Li_layers)
        TM_layers = self.group_atoms_by_layer(TM_layers)
        O_layers = self.group_atoms_by_layer(O_layers)

        # Get the fractional coordinates of distinct vertical layers for Metals and Oxygen
        M_layers_c_frac = sorted({**Li_layers, **TM_layers})
        O_layers_c_frac = sorted(O_layers)

        # Collect the central layers based on the "selective_dynamics" tag: atoms with the "F" property are assigned into the central layers 
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
        # Questions: Is it okay that the terminated metal species are determined by the number of Li and TM atoms as follows? I think this algorithm depends on the assumption that the given structure is "Layered" and the number of metal atoms on each layer is the same.
        if len(Li_layers) > len(TM_layers):
            c_max_surface_Li = sorted(Li_layers.items())[-1][0]
            c_max_surface_O = sorted(O_layers.items())[-1][0]
            return [lower_limit, upper_limit, c_max_surface_Li,
                    c_max_surface_O, Li_layers, M_layers_c_frac,
                    O_layers_c_frac, TM_layers]
        # TM terminated slab model -- (001) surface
        elif len(Li_layers) < len(TM_layers):
            c_surface_TM = sorted(TM_layers.items())[-1][0]
            c_max_surface_O = sorted(O_layers.items())[-1][0]
            return [lower_limit, upper_limit, c_surface_TM,
                    c_max_surface_O, Li_layers, M_layers_c_frac,
                    O_layers_c_frac, TM_layers]
        # Li-TM-O terminated slab model -- (104) surface
        else:
            c_surface_TM_Li = sorted(Li_layers.items())[-1][0]
            c_max_surface_O = sorted(O_layers.items())[-1][0]
            return [lower_limit, upper_limit, c_surface_TM_Li,
                    c_max_surface_O, Li_layers, M_layers_c_frac,
                    O_layers_c_frac, TM_layers]

    def index_extraction(self,
                         tol=0.01,):
        """
        This function is used to get the index of atoms in the relaxed
        surface layers.

        Args:
            tol: Maximum spread (in fractional coordinates) that atoms in the
                same plane may exhibit. Defaults to 0.01.

        Returns:
            Index of Li and O atoms in the relaxed surface layers.

        """
        [_, _, _, c_max_surface_O,
         _, _, O_layers_c_frac, _] = self.layer_distinguisher()

        # Get indices of the first and second Li layers as well as the surface
        # oxygen layers
        relaxed_li_index, oxygen_index = [], []
        # Questions 1: The following lines are only applicable to Li-TM-O assuming that the Li and TM are mixed. Do you think it is enough to only consider the ternary Li-TM-O?
        # Questions 2: The following lines are assuming that the surface normal direction is parallel to the c lattice. However, the possibility that the surface normal direction is along a or b lattice should also be considered using the "slab_direction" argument as in self.wrap_pbc
        for s in self:
        for index, site in enumerate(self):
            # Surface Li atoms index extraction
            if 'Li' in site:
                if site.properties["selective_dynamics"] == [True, True, True]:
                    relaxed_li_index.append(index)

            # Surface O atoms index extraction (both top and bottom surface O
            # atoms are all included)
            elif 'O' in site:
                if (site.frac_coords[2] >= c_max_surface_O - tol) or \
                        (site.frac_coords[2] <= O_layers_c_frac[0] + tol):  # Questions 3: Although the "O_layers" and "O_layers_c_frac" are given after calling the self.group_atoms_by_layer in self.layer_distinguisher, why this "tol" is needed? Or I think the use of the same value as "max_diff" in self.group_atoms_by_layer is more appropriate than this independent "tol" value.
                    oxygen_index.append(index)

        return relaxed_li_index, oxygen_index

    def remove_sites_with_scaling(self,
                                  index,
                                  scaling=True,
                                  scaling_matrix=None):
        """
        This function is used to remove sites in the slab model and
        then create a supercell based on the scaling matrix.

        Args:
            index: Indices of the atoms that are going to be removed.
            scaling: Whether to create a supercell of the slab model with
                existing vacancies. Defaults to True.
            scaling_matrix: Scaling matrix to be used to create the
                supercell. Defaults to None.

        Returns:
            Supercell slab model with sites removed on the surface.

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
                            direction=2,  # Questions 1: How about defining the slab_direction as instance varialbe of Slab Class instead of separatly assigning in many instance functions
                            tol=0.01):
        """
        This function is used to symmetrize the top surface enumerated slab
        model on the basis of the inversion symmetry center to create a
        fully enumerated slab model. The inversion symmetry center is
        determined based on the central fixed region only because the input
        slab model is not symmetric since on the one hand, the top surface
        of the slab is substituted by other dummy species,
        on the other hand, it has been enumerated with different Li and O
        compositions on the surface. Therefore, using the whole slab to
        detect the inversion symmetry does not work.

        Args:
            symprec: Tolerance for the symmetry detection. Defaults to 1e-4.
            direction: Lattice direction perpendicular to the surface,
                i.e. parallel to the c lattice parameter. Defaults to 2.
            tol: Tolerance to determine the central fixed region. Defaults
                to 0.01.

        Returns:
            Symmetrized slab model with equivalent top and bottom
            surfaces.
        """

        # Load slab which is going to be symmetrized
        slab_tgt = self  # Questions 2: copy.deepcopy is not needed in this part?

        # Use center region to get the inversion symmetry center and the
        # rotation matrix
        lower_limit, upper_limit, _, _, _, _, _, _ = self.layer_distinguisher()

        # Generate the reference slab which is just the central fixed region
        slab_ref_site = []
        for s in slab_tgt:
            if lower_limit - tol < s.frac_coords[2] < upper_limit + tol:
                slab_ref_site.append(copy.deepcopy(s))
        slab_ref = Structure.from_sites(slab_ref_site) # Questions 3: "center_sites" is already defined in self.layer_distinguisher. If same variables are defined more than two times, it would be better to define the variable as instance variable and share it in multiple instance functions.

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

        # Initialize a new slab model which only has the center and top
        # regions of the initial slab model
        sites = center_sites[:] + top_sites[:]
        for s in top_sites:
            s2 = copy.deepcopy(s)
            s2.frac_coords = inversion.operate(s.frac_coords)
            # Add symmetrized top sites into the new slab as the bottom sites
            sites.append(s2)

        symmetrized_slab_top = Slab.from_sites(sites)
        symmetrized_slab_top = symmetrized_slab_top.wrap_pbc()
        symmetrized_slab_top = symmetrized_slab_top.get_sorted_structure()

        return symmetrized_slab_top

    def surface_substitute(self,
                           subs1: str,
                           subs2: str,
                           direction=2,
                           tol=0.02):
        """
        This function is used to substitute surface Li and O atoms with "dummy
        species". This will facilitate the EnumWithComposition class to detect
        to-be-enumerated atoms.

        Args:
            subs1: Substitution atom for Li atom.
            subs2: Substitution atom for O atom.
            direction: Define the direction which is parallel to the c
                lattice parameter. Defaults to 2.
            tol: Maximum spread (in fractional coordinates) that atoms in the
                same plane may exhibit. Defaults to 0.02.

        Returns:
            Substituted slab model

        """
        # Load slab
        slab_tgt = self

        # Define c-fractional coordinates of upper boundary of the fixed
        # region, surface metal, and surface O atoms.
        [_, center_top, c_max_surface_metal,
         c_max_surface_O, _, _, _, _] = self.layer_distinguisher()

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
        all sites. This is used to find the location of the slab region in
        the model.

        Returns:
            maximum and minimum c fractional coordinates of a slab model

        """
        minimum, maximum = 1, 0
        for site in self:
            if site.frac_coords[2] > maximum:
                maximum = site.frac_coords[2]
            if site.frac_coords[2] < minimum:
                minimum = site.frac_coords[2]
        return minimum, maximum

    def Li_TM_layers_finder(self):
        """
        This function is used to find the c fractional coordinates of Li, TM,
        and O layers as well as number of atoms in each layer.

        Returns:
            Dict which c fractional coordinates as the key and number of
            atoms as the value

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

        Li_layers = self.group_atoms_by_layer(Li_layers)
        TM_layers = self.group_atoms_by_layer(TM_layers)
        O_layers = self.group_atoms_by_layer(O_layers)
        return Li_layers, TM_layers, O_layers
