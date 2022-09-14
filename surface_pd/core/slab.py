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
from surface_pd.util import check_int


class Slab(Structure):
    """
    Constructor of periodic Slab class.
    This child class is inherited from the pymatgen.core.structure.Structure.
    The args defined have the same meanings as in the pymatgen
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
            dict of sequences. Defaults to None.
        _direction (int): Lattice direction perpendicular to the surface,
            i.e. parallel to the c lattice parameter. Defaults to 2.
        _tolerance (float): If the distance between two sites in c direction is
            larger than this tolerance value, these two sites will be
            treated as located on two difference layers.
        _to_be_enumerated_species (list):
        _num_layers_relaxed (dict):
        _symmetric (bool):

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
                 _direction: int = 2,
                 _tolerance: float = 0.03,
                 _to_be_enumerated_species: list = None,
                 _num_layers_relaxed: dict = None,
                 _symmetric: bool = None
                 ):
        super().__init__(lattice, species, coords,
                         charge, validate_proximity,
                         to_unit_cell, coords_are_cartesian,
                         site_properties)
        self._direction = _direction
        self._tolerance = _tolerance
        self._to_be_enumerated_species = _to_be_enumerated_species
        self._num_layers_relaxed = _num_layers_relaxed
        self._symmetric = _symmetric

    @property
    def direction(self):
        return self._direction

    @direction.setter
    def direction(self, a: int):
        self._direction = a

    @property
    def tolerance(self):
        return self._tolerance

    @tolerance.setter
    def tolerance(self, a: float):
        self._tolerance = a

    @property
    def to_be_enumerated_species(self):
        return self._to_be_enumerated_species

    @to_be_enumerated_species.setter
    def to_be_enumerated_species(self, a: list):
        self._to_be_enumerated_species = a

    @property
    def num_layers_relaxed(self):
        return self._num_layers_relaxed

    @num_layers_relaxed.setter
    def num_layers_relaxed(self, a: dict):
        self._num_layers_relaxed = a

    @property
    def symmetric(self):
        return self._symmetric

    @symmetric.setter
    def symmetric(self, a: bool):
        self._symmetric = a

    def group_atoms_by_layer(self,
                             layers: dict):
        """
        Group misclassified atoms into the right layers.
        For example, c_atom1 = 0.01, c_atoms2 = 0.02, they should be
        classified in one layer. But actually they are not. So this function
        here will search the difference between two closest atoms, if the
        difference is smaller than diff, they will be regrouped in the same
        layer.

        Args:
            layers: Dict which contains different c fractional
              coordinates as keys and number of atoms as values.

        Returns:
            Dict which contains different c fractional
            coordinates as keys and number of atoms as values.

        """

        res = dict()
        o_layer_sorted = dict(
            sorted(layers.items(), key=lambda x: x[0], reverse=True))
        for i, (height_, num_) in enumerate(o_layer_sorted.items()):
            if i == 0:
                res[height_] = num_
                previous_height = height_
                continue
            if abs(previous_height - height_) <= self.tolerance:
                res[previous_height] += num_
            else:
                res[height_] = num_
                previous_height = height_
        return res

    def wrap_pbc(self):
        """
        Wrap fractional coordinates back into the unit cell.

        Returns:
            Slab model with all sites in the unit cell.

        """
        EPS = np.finfo(np.double).eps
        # min_c, _ = self.get_max_min_c_frac()
        # if min_c < -0.01:
        #     shift = abs(min_c) + 0.01
        # else:
        shift = 0

        for site in self:
            bounds = np.array([0.0, 0.0, 0.0])
            bounds[self.direction] = self.tolerance
            for i in range(3):
                while site.frac_coords[i] < -bounds[i]:
                    site.frac_coords[i] += 1.0
                while site.frac_coords[i] >= (1.0 - bounds[i] - EPS):
                    site.frac_coords[i] -= 1.0
            site.frac_coords[self.direction] += shift
        return self

    def get_center_sites(self):
        """
        Get the center sites in the slab model on the basis of the
        "selective_dynamics" defined.

        Returns:
            lower and upper limits of the central fixed region, and the
            central fixed structure
        """

        center_sites = []
        for s in self:
            if ('selective_dynamics' in s.properties
                    and not any(s.properties['selective_dynamics'])):
                center_sites.append(copy.deepcopy(s))

        # Get the central slab boundaries
        center = Slab.from_sites(center_sites)
        center_region = get_slab_regions(center)
        lower_limit, upper_limit = center_region[0]
        return lower_limit, upper_limit, center

    def layers_finder(self,
                      precision: int = 2):
        """
        This function is used to find the c fractional coordinates of all
        species in the input slab model and the number of atoms in each
        layer.

        Args:
            precision: Round c fraction coordinates to a given precision in
                decimal digits. Defaults to 2.

        Returns:
            {"species": {"c fractional coordinates": number of atoms at here}}

        """
        layers = {}
        for species in self.symbol_set:
            partial_layers = collections.defaultdict(int)
            for s in self:
                if species in s:
                    # Since the c fractional coordinates of atoms even in the
                    # exactly same layers are not identical, therefore, all c
                    # fractional coordinates will be rounded to 2 decimal. Any
                    # misclassified atoms will be regrouped by
                    # "group_atoms_by_layer" function.
                    partial_layers[
                        round(s.frac_coords[self.direction], precision)] += 1
            partial_layers = self.group_atoms_by_layer(partial_layers)
            layers[species] = partial_layers
        return layers

    def layer_distinguisher(self):
        """
        This function is used to distinguish the layers that will be relaxed.

        Returns:
            c fractional coordinates of the upper and lower boundaries
            of the central fixed region. (1, 2)\n
            {"species": [top and bottom c fractional coordinates]}

        """

        # Initialize three dictionaries to store layer information of Li,
        # TM, and O atoms
        layers = self.layers_finder()

        # Get the fractional coordinates of distinct vertical layers for
        # metals and oxygen
        target_layers = {}
        for species, c_frac in zip(self.to_be_enumerated_species,
                                   self.num_layers_relaxed):
            all_c_frac = list(layers[species])
            if self.symmetric:

                target_layers[species] = \
                    [all_c_frac[self.num_layers_relaxed[species] - 1],
                     all_c_frac[-self.num_layers_relaxed[species]]]
            else:
                target_layers[species] = \
                    [all_c_frac[self.num_layers_relaxed[species] - 1]]

        # Collect the central layers based on the "selective_dynamics"
        lower_limit, upper_limit, _ = self.get_center_sites()

        return lower_limit, upper_limit, target_layers

    def index_extraction(self,
                         only_top: bool = False):
        """
        This function is used to get the indices of the target species in the
        relaxed surface layers.

        Args:
            only_top: Whether to only return indices of the target species in
                the top relaxed surface layer.

        Returns:
            {"species": [indices of the target species]}

        """

        [lower_limit, upper_limit, target_layers] = self.layer_distinguisher()

        relaxed_index_by_species = {}
        for species in self.to_be_enumerated_species:
            relaxed_index = []
            for index, site in enumerate(self):
                # Surface Li atoms index extraction (both top and bottom)
                if species in site:
                    if (site.frac_coords[self.direction] >=
                            target_layers[species][0] - self.tolerance):
                        relaxed_index.append(index)
                    if self.symmetric:
                        if only_top:
                            pass
                        else:
                            if (site.frac_coords[self.direction] <=
                                    target_layers[species][1] +
                                    self.tolerance):
                                relaxed_index.append(index)
                    else:
                        pass
            relaxed_index_by_species[species] = relaxed_index

        return lower_limit, upper_limit, relaxed_index_by_species

    def surface_substitute(self,
                           dummy_species: list):
        """
        This function is used to substitute the to-be-enumerated species from
        the slab top surface with dummy species. This will facilitate the
        EnumWithComposition class to find the targeted species.

        Args:
            dummy_species: A special specie for representing non-traditional
                elements or species.

        Returns:
            Substituted slab model

        """
        slab_surface_substitute = copy.deepcopy(self)
        _, _, relaxed_index = slab_surface_substitute.index_extraction(
            only_top=True)

        for i, species in enumerate(self.to_be_enumerated_species):
            for j in relaxed_index[species]:
                slab_surface_substitute.replace(
                    j, dummy_species[i],
                    properties={'selective_dynamics': [True, True, True]}
                )

        return slab_surface_substitute

    def supplemental_structures_gene(self,
                                     subs_dict: dict,
                                     relaxed_index: dict):
        """
        This function is used to generate supplemental structures.
        Supplemental structures include the cases that the to-be-enumerated
        structure is straight forward to generate and does not require the
        EnumWithComposition class to involve.

        Args:
            subs_dict: Species and occupancy dictionaries containing the
                species mapping in string-string pairs.
            relaxed_index: The indices of the to-be-enumerated atoms on the
                surface.

        Returns:
            New slab model with some sites removed.

        """
        structure = copy.deepcopy(self)
        removed_index = []
        for species in self.to_be_enumerated_species:
            if subs_dict[species][species] == 0:
                removed_index += relaxed_index[species]
        structure.remove_sites(removed_index)
        return structure

    def symmetrize_top_base(self,
                            symprec: float = 1e-4):
        """
        Symmetrize the top surface enumerated slab
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

        Returns:
            Symmetrized slab model with equivalent top and bottom
            surfaces.
        """

        # Generate the reference slab which is just the central fixed region
        _, _, slab_ref = self.get_center_sites()

        # Determine symmetry operations of the reference slab and make sure
        # the reference slab has an inversion center
        try:
            symmetric, origin, inversion = slab_ref.is_symmetry(
                symprec=symprec, return_isc=True)
        except TypeError:
            symmetric = slab_ref.is_symmetry(symprec, return_isc=False)

        if not symmetric:
            slab_ref.to(fmt='poscar', filename='debug-center.vasp')
            print("Please see the saved \"debug-center.vasp\" structure and "
                  "see if it mases sense.")
            raise NoInversionSymmetryError

        # New way to implement wrap_pbc
        self.wrap_pbc()

        # Symmetrized slab models based on top only
        top_sites = []
        bottom_sites = []
        center_sites = slab_ref.sites

        for s in [s for s in self
                  if s.frac_coords[self.direction] > origin[self.direction]]:
            if ('selective_dynamics' in s.properties and
                    any(s.properties['selective_dynamics'])):
                top_sites.append(copy.deepcopy(s))
        for s in [s for s in self
                  if s.frac_coords[self.direction] < origin[self.direction]]:
            if ('selective_dynamics' in s.properties and
                    any(s.properties['selective_dynamics'])):
                bottom_sites.append(copy.deepcopy(s))

        # # Extend the center sites list in case of some sites have exactly the
        # # same c fractional coordinates of the origin
        # center_sites.extend(
        #     [copy.deepcopy(s) for s in self
        #      if s.frac_coords[self.direction] ==
        #      origin[self.direction]]
        # )

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
            if site.frac_coords[self.direction] > maximum:
                maximum = site.frac_coords[self.direction]
            if site.frac_coords[self.direction] < minimum:
                minimum = site.frac_coords[self.direction]
        return minimum, maximum

    def is_symmetry(self,
                    symprec: float = 1e-1,
                    return_isc: bool = False):
        """
        Check whether the slab model is symmetric. If it is, the inversion
        symmetry center and the inversion operation are able to return.

        Args:
            symprec: Tolerance for symmetry finding. Defaults to 1e-1.
            return_isc: Whether to return the inversion symmetry center and
                the inversion operation.

        """
        sga = SpacegroupAnalyzer(self, symprec=symprec)
        if sga.is_laue():  # has Laue symmetry (centro-symmetry)
            if return_isc:
                ops = sga.get_symmetry_operations()
                inversion = ops[1]
                assert (np.all(inversion.rotation_matrix == -np.identity(3)))
                origin = inversion.translation_vector / 2
                return True, origin, inversion
            else:
                return True
        else:
            return False

    def check_rotate(self):
        """
        Check if the enumerated slab models need to rotate to satisfy the
        shape (cuboid) requirement.

        Returns:
            Cuboid slab model

        """
        if max(self.lattice.abc) != self.lattice.c:
            if max(self.lattice.abc) == self.lattice.a:
                slab_rotated = copy.deepcopy(self)
                slab_rotated.make_supercell(
                    [[0, 0, 1],
                     [0, 1, 0],
                     [1, 0, 0]]
                )
            elif max(self.lattice.abc) == self.lattice.b:
                slab_rotated = copy.deepcopy(self)
                slab_rotated.make_supercell(
                    [[1, 0, 0],
                     [0, 0, 1],
                     [0, 1, 0]]
                )
            return slab_rotated
        else:
            return self

    def calculate_num_sites(self,
                            composition_list: list,
                            relaxed_index: dict,
                            max_cell_size: int):
        """
        Calculate the num of sites should be in the slab model on the basis
        of the user defined composition.

        Args:
            composition_list: Compositions of enumerated species.
            relaxed_index: {"species": [indices of the target species]}
            max_cell_size: Maximum number of supercells of the input slab.

        Returns:
            Number of sites in the slab model.

        """
        num_relaxed_atoms = {}
        total_target_atoms = {}
        num_rest_atoms = {}
        num_curr_atoms = {}
        for i, species in enumerate(self.to_be_enumerated_species):
            num_relaxed_atoms[species] = len(
                relaxed_index[species]) * max_cell_size
            total_target_atoms[species] = \
                self.composition[species] * max_cell_size
            num_rest_atoms[species] = \
                total_target_atoms[species] - num_relaxed_atoms[species]
            num_curr_atoms[species] = \
                num_relaxed_atoms[species] * composition_list[i] + \
                num_rest_atoms[species]

        total_rest_atoms = self.num_sites * max_cell_size - sum(
            list(total_target_atoms.values()))

        curr_num_sites = check_int(sum(list(num_curr_atoms.values())) +
                                   total_rest_atoms)
        return curr_num_sites

    def tune_isc(self,
                 origin,
                 shift_isc_back: bool = True):
        """
        Shift or shift back the inversion symmetry center.

        Args:
            origin: Inversion symmetry center.
            shift_isc_back: Whether to shift the inversion symmetry center to
                the origin.

        Returns:
            Slab model with inversion symmetry center shifted to the origin /
                back.

        """
        if shift_isc_back:
            for site in self:
                site.frac_coords = site.frac_coords + origin
            self.wrap_pbc()
        else:
            symmetric, origin, _ = self.is_symmetry(
                symprec=0.1,
                return_isc=True)
            for site in self:
                site .frac_coords -= origin
            self.wrap_pbc()
        return self

    def tune_c(self,
               target_min_c):
        """
        Slightly adjust the slab along the slab direction defined before.
        This is to make sure that the central region is still located at
        around the same place.

        Args:
            target_min_c: Desired minimum fractional coordinate along the
                slab direction defined before.

        Returns:
            Slight adjusted slab model

        """
        _min_c, _ = self.get_max_min_c_frac()
        displace = target_min_c - _min_c

        if abs(displace) < self.tolerance:
            displace = 0

        for site in self:
            site.frac_coords[self.direction] += displace
        return self

    def add_selective_dynamics(self,
                               lower_limit,
                               upper_limit):
        """
        Add selective dynamics to the after refined slab model based on the
        number of layers that will be relaxed on the surface.

        Args:
            lower_limit: Lower limit of the central fixed region
            upper_limit: Upper limit of the central fixed region

        Returns:
            Enumerated slab model with "selective dynamics" labeled.

        """
        direction = self.direction
        tolerance = self.tolerance
        structure = copy.deepcopy(self)
        for site in structure:
            if lower_limit - tolerance <= site.frac_coords[direction] <= \
                    upper_limit + tolerance:
                site.properties = {'selective_dynamics': [False, False, False]}
            else:
                site.properties = {'selective_dynamics': [True, True, True]}
        return structure

    # def selective_dynamics_completion(self,
    #                                   dummy_species: list,
    #                                   center_bottom: float,
    #                                   center_top: float,
    #                                   ):
    #     for t in self:
    #         for ds in dummy_species:
    #             if ds in t:
    #                 t.properties = {
    #                     'selective_dynamics': [True, True, True]}
    #         if t.properties['selective_dynamics'] is None:
    #             if (center_bottom - self.tolerance <=
    #                     t.frac_coords[self.direction] <=
    #                     center_top + self.tolerance):
    #                 t.properties = {
    #                     'selective_dynamics': [False, False, False]}
    #             else:
    #                 t.properties = {
    #                     'selective_dynamics': [True, True, True]}
    #     return self
