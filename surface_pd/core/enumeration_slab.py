"""Enumeration slab representation and manipulation.

This module extends pymatgen's ``Structure`` with layer identification,
symmetry operations, and surface-enumeration configuration.
"""

import collections
import copy
import logging
import math
from collections.abc import Mapping, Sequence
from numbers import Integral, Real

import numpy as np
from pymatgen.core import Composition, DummySpecies, Element, Species
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.surface import get_slab_regions
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from surface_pd.error import NoInversionSymmetryError
from surface_pd.util.util import check_int

logger = logging.getLogger(__name__)


class EnumerationSlab(Structure):
    """Represent and manipulate a periodic surface slab.

    Parameters
    ----------
    lattice : Lattice or array-like of shape (3, 3)
        Periodic lattice accepted by :class:`pymatgen.core.Structure`.
    species : sequence
        Species for the sites, in the same order as ``coords``.
    coords : sequence of sequence of float
        Fractional coordinates unless ``coords_are_cartesian`` is true.
    charge : float, optional
        Overall structure charge.
    validate_proximity : bool, default=False
        Reject sites separated by less than pymatgen's proximity threshold.
    to_unit_cell : bool, default=False
        Wrap fractional coordinates into the unit cell during construction.
    coords_are_cartesian : bool, default=False
        Interpret ``coords`` as Cartesian coordinates when true.
    site_properties : dict, optional
        Per-site property sequences accepted by pymatgen.
    labels : sequence of str or None, optional
        Optional label for each site, accepted by pymatgen.
    properties : dict, optional
        Structure-level properties accepted by pymatgen.
    direction : int, default=2
        Lattice-axis index normal to the surface. Must be 0, 1, or 2.
    tolerance : float, default=0.03
        Positive finite fractional-coordinate tolerance used to group layers.
    enumerated_species : sequence of str, optional
        Unique, nonempty species names whose surface sites will be enumerated.
    num_enumerated_layers : mapping of str to int, optional
        Positive number of outer layers to enumerate for each target species.
    symmetric : bool, optional
        Whether enumeration operates on both slab surfaces. ``None`` denotes
        an unconfigured slab.

    Notes
    -----
    Surface-pd configuration parameters are keyword-only. They initialize
    validated public properties while their values are stored privately.

    """

    def __init__(
        self,
        lattice: list | np.ndarray | Lattice,
        species: Sequence[
            str | Element | Species | DummySpecies | Composition
        ],
        coords: Sequence[Sequence[float]],
        charge: float = None,
        validate_proximity: bool = False,
        to_unit_cell: bool = False,
        coords_are_cartesian: bool = False,
        site_properties: dict = None,
        labels: Sequence[str | None] | None = None,
        properties: dict | None = None,
        *,
        direction: int = 2,
        tolerance: float = 0.03,
        enumerated_species: Sequence[str] | None = None,
        num_enumerated_layers: Mapping[str, int] | None = None,
        symmetric: bool | None = None,
    ):
        super().__init__(
            lattice,
            species,
            coords,
            charge,
            validate_proximity,
            to_unit_cell,
            coords_are_cartesian,
            site_properties,
            labels,
            properties,
        )
        self.direction = direction
        self.tolerance = tolerance
        self.enumerated_species = enumerated_species
        self.num_enumerated_layers = num_enumerated_layers
        self._validate_layer_species_match()
        self.symmetric = symmetric

    @classmethod
    def from_structure(
        cls,
        structure: Structure,
        *,
        direction: int = 2,
        tolerance: float = 0.03,
        enumerated_species: Sequence[str] | None = None,
        num_enumerated_layers: Mapping[str, int] | None = None,
        symmetric: bool | None = None,
    ):
        """Construct an enumeration slab from a pymatgen structure.

        Parameters
        ----------
        structure : Structure
            Source structure whose lattice, sites, charge, labels, site
            properties, and structure properties are copied.
        direction : int, default=2
            Lattice-axis index normal to the surface.
        tolerance : float, default=0.03
            Fractional-coordinate tolerance used to group layers.
        enumerated_species : sequence of str, optional
            Species whose outer surface layers may be modified.
        num_enumerated_layers : mapping of str to int, optional
            Independent outer-layer count for each enumerated species.
        symmetric : bool, optional
            Whether corresponding layers on both surfaces are selected.

        Returns
        -------
        EnumerationSlab
            Independent copy configured for surface enumeration.
        """
        if not isinstance(structure, Structure):
            raise TypeError("structure must be a pymatgen Structure")
        return cls(
            structure.lattice,
            structure.species,
            structure.frac_coords,
            charge=structure.charge,
            site_properties=structure.site_properties,
            labels=structure.labels,
            properties=structure.properties,
            direction=direction,
            tolerance=tolerance,
            enumerated_species=enumerated_species,
            num_enumerated_layers=num_enumerated_layers,
            symmetric=symmetric,
        )

    @property
    def direction(self):
        """int: Lattice-axis index perpendicular to the surface.

        Valid axis indices are 0, 1, and 2.
        """
        return self._direction

    @direction.setter
    def direction(self, value: int):
        if isinstance(value, bool) or not isinstance(value, Integral):
            raise TypeError("direction must be an integer")
        if value not in (0, 1, 2):
            raise ValueError("direction must be 0, 1, or 2")
        self._direction = int(value)

    @property
    def tolerance(self):
        """float: Fractional-coordinate tolerance used to group layers.

        The value must be positive and finite.
        """
        return self._tolerance

    @tolerance.setter
    def tolerance(self, value: float):
        if (
            isinstance(value, bool)
            or not isinstance(value, Real)
            or not math.isfinite(value)
        ):
            raise TypeError("tolerance must be a finite real number")
        if value <= 0:
            raise ValueError("tolerance must be positive")
        self._tolerance = float(value)

    @property
    def enumerated_species(self):
        """List of str or None: Species selected for enumeration.

        Values must be unique, nonempty strings. Assigned sequences are copied.
        """
        return self._enumerated_species

    @enumerated_species.setter
    def enumerated_species(self, value: Sequence[str] | None):
        if value is None:
            self._enumerated_species = None
            return
        if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
            raise TypeError("enumerated_species must be a sequence")
        species = list(value)
        if not species:
            raise ValueError("enumerated_species must not be empty")
        if any(
            not isinstance(item, str) or not item.strip() for item in species
        ):
            raise ValueError("enumerated species must be nonempty strings")
        if len(set(species)) != len(species):
            raise ValueError("enumerated species must be unique")
        self._enumerated_species = species
        self._validate_layer_species_match()

    @property
    def num_enumerated_layers(self):
        """Dict or None: Number of enumerated layers for each target species.

        Keys must be nonempty strings and values positive integers. Assigned
        mappings are copied.
        """
        return self._num_enumerated_layers

    @num_enumerated_layers.setter
    def num_enumerated_layers(self, value: Mapping[str, int] | None):
        if value is None:
            self._num_enumerated_layers = None
            return
        if not isinstance(value, Mapping):
            raise TypeError("num_enumerated_layers must be a mapping")
        if not value:
            raise ValueError("num_enumerated_layers must not be empty")
        layers = {}
        for species, count in value.items():
            if not isinstance(species, str) or not species.strip():
                raise ValueError("layer species must be nonempty strings")
            if isinstance(count, bool) or not isinstance(count, Integral):
                raise TypeError("enumerated layer counts must be integers")
            if count <= 0:
                raise ValueError("enumerated layer counts must be positive")
            layers[species] = int(count)
        self._num_enumerated_layers = layers
        self._validate_layer_species_match()

    def _validate_layer_species_match(self):
        species = getattr(self, "_enumerated_species", None)
        layers = getattr(self, "_num_enumerated_layers", None)
        if (
            species is not None
            and layers is not None
            and set(species) != set(layers)
        ):
            raise ValueError(
                "num_enumerated_layers must contain the same species as "
                "enumerated_species"
            )

    @property
    def symmetric(self):
        """Bool or None: Whether both slab surfaces are enumerated.

        The value must be a boolean or ``None``.
        """
        return self._symmetric

    @symmetric.setter
    def symmetric(self, value: bool | None):
        if value is not None and not isinstance(value, bool):
            raise TypeError("symmetric must be a boolean or None")
        self._symmetric = value

    def group_atoms_by_layer(self, layers: dict):
        """
        Group misclassified atoms into the right layers.

        For example, c_atom1 = 0.01, c_atoms2 = 0.02, they should be
        classified in one layer. But actually they are not. So this function
        here will search the difference between two closest atoms, if the
        difference is smaller than diff, they will be regrouped in the same
        layer.

        Parameters
        ----------
        layers : dict[float, int]
            Fractional heights mapped to atom populations.

        Returns
        -------
        dict[float, int]
            Layer populations after adjacent heights within ``tolerance`` are
            merged. The input mapping is not mutated.

        """
        res = {}
        o_layer_sorted = dict(
            sorted(layers.items(), key=lambda x: x[0], reverse=True)
        )
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

        Returns
        -------
        EnumerationSlab
            This slab after its site coordinates have been wrapped.

        Notes
        -----
        This method mutates the receiver and returns it for fluent use.

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
        Get center sites from fixed selective dynamics flags.

        Returns
        -------
        tuple[float, float, EnumerationSlab]
            Lower and upper fractional-coordinate bounds of the fixed region,
            followed by a new slab containing its fixed sites.
        """
        center_sites = []
        for s in self:
            if "selective_dynamics" in s.properties and not any(
                s.properties["selective_dynamics"]
            ):
                center_sites.append(copy.deepcopy(s))

        # Get the central slab boundaries
        center = EnumerationSlab.from_sites(center_sites)
        center_region = get_slab_regions(center)
        lower_limit, upper_limit = center_region[0]
        return lower_limit, upper_limit, center

    def layers_finder(self, precision: int = 2):
        """
        Find species layer populations by fractional coordinate.

        Parameters
        ----------
        precision : int, default=2
            Decimal places used to round fractional layer coordinates.

        Returns
        -------
        dict[str, dict[float, int]]
            Species mapped to fractional heights and atom populations.

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
                        round(s.frac_coords[self.direction], precision)
                    ] += 1
            partial_layers = self.group_atoms_by_layer(partial_layers)
            layers[species] = partial_layers
        return layers

    def layer_distinguisher(self):
        r"""
        Distinguish the layers that will be relaxed.

        Returns
        -------
        tuple[float, float, dict[str, list[float]]]
            Lower and upper bounds of the fixed region, followed by target
            layer heights for each enumerated species. Symmetric slabs contain
            top and bottom heights; asymmetric slabs contain only the top.

        """
        # Initialize three dictionaries to store layer information of Li,
        # TM, and O atoms
        layers = self.layers_finder()

        # Get the boundaries of the surface relaxed region in fractional
        # coordinates format (c-direction)
        target_layers = {}
        for species in self.enumerated_species:
            all_c_frac = list(layers[species])
            if self.symmetric:
                target_layers[species] = [
                    all_c_frac[self.num_enumerated_layers[species] - 1],
                    all_c_frac[-self.num_enumerated_layers[species]],
                ]
            else:
                target_layers[species] = [
                    all_c_frac[self.num_enumerated_layers[species] - 1]
                ]

        # Collect the central layers based on the "selective_dynamics"
        lower_limit, upper_limit, _ = self.get_center_sites()

        return lower_limit, upper_limit, target_layers

    def get_fixed_region_bounds(self):
        """Return fractional-coordinate bounds of the fixed slab region.

        Returns
        -------
        tuple[float, float]
            Lower and upper surface-normal bounds of the fixed central sites.
        """
        lower_limit, upper_limit, _ = self.get_center_sites()
        return lower_limit, upper_limit

    def get_enumerated_site_indices(self, only_top: bool = False):
        """
        Get target species indices in the relaxed surface layers.

        Parameters
        ----------
        only_top : bool, default=False
            Return only top-surface indices when true.

        Returns
        -------
        dict[str, list[int]]
            Selected site indices grouped by enumerated species.

        """
        _, _, target_layers = self.layer_distinguisher()

        relaxed_index_by_species = {}
        for species in self.enumerated_species:
            relaxed_index = []
            for index, site in enumerate(self):
                # Surface Li atoms index extraction (both top and bottom)
                if species in site:
                    if (
                        site.frac_coords[self.direction]
                        >= target_layers[species][0] - self.tolerance
                    ):
                        relaxed_index.append(index)
                    if (
                        self.symmetric
                        and not only_top
                        and site.frac_coords[self.direction]
                        <= target_layers[species][1] + self.tolerance
                    ):
                        relaxed_index.append(index)
            relaxed_index_by_species[species] = relaxed_index

        return relaxed_index_by_species

    def _surface_substitute(self, dummy_species: list):
        """
        Substitute top-surface target species with dummy species.

        Parameters
        ----------
        dummy_species : list
            One pymatgen-compatible dummy species per enumerated species.

        Returns
        -------
        EnumerationSlab
            Deep copy with top-surface target sites replaced. The receiver is
            not mutated.

        """
        slab_surface_substitute = copy.deepcopy(self)
        relaxed_index = slab_surface_substitute.get_enumerated_site_indices(
            only_top=True
        )

        for i, species in enumerate(self.enumerated_species):
            for j in relaxed_index[species]:
                slab_surface_substitute.replace(
                    j,
                    dummy_species[i],
                    properties={"selective_dynamics": [True, True, True]},
                )

        return slab_surface_substitute

    def generate_supplemental_structures(
        self, subs_dict: dict, relaxed_index: dict
    ):
        """
        Generate supplemental structures.

        Supplemental structures include the cases that the to-be-enumerated
        structure is straight forward to generate and does not require the
        EnumWithComposition class to involve.

        Parameters
        ----------
        subs_dict : dict
            Target species mapped to replacement-species occupancies.
        relaxed_index : dict[str, list[int]]
            Target species mapped to surface-site indices.

        Returns
        -------
        EnumerationSlab
            Deep copy with sites removed where the target species has zero
            occupancy. The receiver is not mutated.

        """
        structure = copy.deepcopy(self)
        removed_index = []
        for species in self.enumerated_species:
            if subs_dict[species][species] == 0:
                removed_index += relaxed_index[species]
        structure.remove_sites(removed_index)
        return structure

    def symmetrize_top_base(self, symprec: float = 1e-4):
        """
        Symmetrize the top surface by the inversion symmetry center.

        This creates a fully enumerated slab model. The inversion symmetry
        center is
        determined based on the central fixed region only because the input
        slab model is not symmetric since on the one hand, the top surface
        of the slab is substituted by other dummy species,
        on the other hand, it has been enumerated with different Li and O
        compositions on the surface. Therefore, using the whole slab to
        detect the inversion symmetry does not work.

        Parameters
        ----------
        symprec : float, default=1e-4
            Cartesian symmetry tolerance in angstroms.

        Returns
        -------
        EnumerationSlab
            New slab with equivalent top and bottom surfaces.

        Raises
        ------
        NoInversionSymmetryError
            If the fixed central region has no inversion symmetry.

        Notes
        -----
        The receiver is wrapped in place before the new slab is constructed.
        On symmetry failure, the central region is written to
        ``debug-center.vasp`` in the current working directory.
        """
        # Generate the reference slab which is just the central fixed region
        _, _, slab_ref = self.get_center_sites()

        # Determine symmetry operations of the reference slab and make sure
        # the reference slab has an inversion center
        try:
            symmetric, origin, inversion = slab_ref.is_symmetry(
                symprec=symprec, return_isc=True
            )
        except TypeError:
            symmetric = slab_ref.is_symmetry(symprec, return_isc=False)

        if not symmetric:
            debug_filename = "debug-center.vasp"
            slab_ref.to(fmt="poscar", filename=debug_filename)
            logger.error(
                "Saved debug structure %s for central slab without inversion "
                "symmetry.",
                debug_filename,
            )
            raise NoInversionSymmetryError

        # New way to implement wrap_pbc
        self.wrap_pbc()

        # Symmetrized slab models based on top only
        top_sites = []
        bottom_sites = []
        center_sites = slab_ref.sites

        for s in [
            s
            for s in self
            if s.frac_coords[self.direction] > origin[self.direction]
        ]:
            if "selective_dynamics" in s.properties and any(
                s.properties["selective_dynamics"]
            ):
                top_sites.append(copy.deepcopy(s))
        for s in [
            s
            for s in self
            if s.frac_coords[self.direction] < origin[self.direction]
        ]:
            if "selective_dynamics" in s.properties and any(
                s.properties["selective_dynamics"]
            ):
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

        symmetrized_slab_top = EnumerationSlab.from_sites(sites)
        symmetrized_slab_top = symmetrized_slab_top.wrap_pbc()
        symmetrized_slab_top = symmetrized_slab_top.get_sorted_structure()

        return symmetrized_slab_top

    def get_max_min_c_frac(self):
        """
        Get the maximum and minimum fractional c coordinates.

        Returns
        -------
        tuple[float, float]
            Minimum and maximum fractional coordinates along ``direction``,
            in that order.

        """
        minimum, maximum = 1, 0
        for site in self:
            if site.frac_coords[self.direction] > maximum:
                maximum = site.frac_coords[self.direction]
            if site.frac_coords[self.direction] < minimum:
                minimum = site.frac_coords[self.direction]
        return minimum, maximum

    def is_symmetry(self, symprec: float = 1e-1, return_isc: bool = False):
        """
        Check whether the slab model is symmetric.

        If it is, the inversion symmetry center and the inversion operation
        are able to return.

        Parameters
        ----------
        symprec : float, default=1e-1
            Cartesian symmetry tolerance in angstroms.
        return_isc : bool, default=False
            Return inversion-center details when true.

        Returns
        -------
        bool or tuple[bool, numpy.ndarray, SymmOp] or None
            With ``return_isc=False``, whether the slab has Laue symmetry.
            With ``return_isc=True``, a true flag, fractional inversion-center
            coordinates, and the pymatgen inversion operation are returned
            when found. ``None`` is returned if Laue symmetry is reported but
            no explicit inversion operation is available.

        """
        sga = SpacegroupAnalyzer(self, symprec=symprec)
        if sga.is_laue():  # has Laue symmetry (centro-symmetry)
            if return_isc:
                ops = sga.get_symmetry_operations()
                for op in ops:
                    if np.all(op.rotation_matrix == -np.identity(3)):
                        inversion = op
                        origin = inversion.translation_vector / 2
                        return True, origin, inversion
                # assert (np.all(inversion.rotation_matrix == -np.identity(3)))
                # origin = inversion.translations_vector / 2
                # return True, origin, inversion
            else:
                return True
        else:
            return False

    def check_rotate(self, criteria: float):
        """
        Check if the slab needs rotation to satisfy the cuboid requirement.

        Parameters
        ----------
        criteria : float
            Target lattice length in angstroms for the surface-normal axis.

        Returns
        -------
        EnumerationSlab or None
            A deep copy, rotated when another lattice axis matches
            ``criteria``. ``None`` is returned when rotation is needed but no
            axis matches.

        """
        directions = [0, 1, 2]
        directions.remove(self.direction)
        slab_rotated = copy.deepcopy(self)

        if (
            round(self.lattice.abc[self.direction], 4) - 0.0001 >= criteria
            or round(self.lattice.abc[self.direction], 4) + 0.0001 <= criteria
        ):
            for d in directions:
                if (
                    round(self.lattice.abc[d], 4) - 0.0001
                    <= criteria
                    <= round(self.lattice.abc[d], 4) + 0.0001
                ):
                    slab_rotated = copy.deepcopy(self)
                    matrix = np.zeros((3, 3))
                    matrix[self.direction, d] = 1
                    matrix[d, self.direction] = 1
                    directions.remove(d)
                    rest = directions[0]
                    matrix[rest, rest] = 1
                    slab_rotated.make_supercell(matrix)
                    return slab_rotated
        else:
            return slab_rotated

    def calculate_num_sites(
        self, composition_list: list, relaxed_index: dict, max_cell_size: int
    ):
        """
        Calculate the number of sites for the user-defined composition.

        Parameters
        ----------
        composition_list : list[float]
            Target occupancy fraction for each enumerated species.
        relaxed_index : dict[str, list[int]]
            Enumerated species mapped to their relaxed surface-site indices.
        max_cell_size : int
            Cell-size multiple used for enumeration.

        Returns
        -------
        int
            Expected number of sites after applying the composition.

        Raises
        ------
        NonIntegerError
            If the calculated site count is not effectively integral.

        """
        num_relaxed_atoms = {}
        total_target_atoms = {}
        num_rest_atoms = {}
        num_curr_atoms = {}
        for i, species in enumerate(self.enumerated_species):
            num_relaxed_atoms[species] = (
                len(relaxed_index[species]) * max_cell_size
            )
            total_target_atoms[species] = (
                self.composition[species] * max_cell_size
            )
            num_rest_atoms[species] = (
                total_target_atoms[species] - num_relaxed_atoms[species]
            )
            num_curr_atoms[species] = (
                num_relaxed_atoms[species] * composition_list[i]
                + num_rest_atoms[species]
            )

        total_rest_atoms = self.num_sites * max_cell_size - sum(
            list(total_target_atoms.values())
        )

        curr_num_sites = check_int(
            sum(list(num_curr_atoms.values())) + total_rest_atoms
        )
        return curr_num_sites

    def tune_isc(self, origin, shift_isc_back: bool = True):
        """
        Shift or shift back the inversion symmetry center.

        Parameters
        ----------
        origin : array-like of shape (3,)
            Fractional translation applied when ``shift_isc_back`` is true.
            Otherwise this argument is ignored and the center is detected.
        shift_isc_back : bool, default=True
            Add ``origin`` when true; detect and subtract the current inversion
            center when false.

        Returns
        -------
        EnumerationSlab
            This slab after shifting and wrapping its coordinates.

        Notes
        -----
        This method mutates the receiver and returns it for fluent use.

        """
        if shift_isc_back:
            for site in self:
                site.frac_coords = site.frac_coords + origin
            self.wrap_pbc()
        else:
            symmetric, origin, _ = self.is_symmetry(
                symprec=0.1, return_isc=True
            )
            for site in self:
                site.frac_coords -= origin
            self.wrap_pbc()
        return self

    def tune_c(self, target_min_c: float):
        """
        Slightly adjust the slab along the slab direction defined before.

        This is to make sure that the central region is still located at
        around the same place.

        Parameters
        ----------
        target_min_c : float
            Desired minimum fractional coordinate along ``direction``.

        Returns
        -------
        EnumerationSlab
            This slab after translating its sites. Displacements smaller than
            ``tolerance`` are suppressed.

        Notes
        -----
        This method mutates the receiver and returns it for fluent use.

        """
        _min_c, _ = self.get_max_min_c_frac()
        displace = target_min_c - _min_c

        if abs(displace) < self.tolerance:
            displace = 0

        for site in self:
            site.frac_coords[self.direction] += displace
        return self

    def add_selective_dynamics(self, lower_limit: float, upper_limit: float):
        """
        Add selective dynamics to a refined slab model.

        Parameters
        ----------
        lower_limit : float
            Lower fractional-coordinate bound of the fixed central region.
        upper_limit : float
            Upper fractional-coordinate bound of the fixed central region.

        Returns
        -------
        EnumerationSlab
            Deep copy with fixed flags inside the central region and relaxed
            flags outside it. The receiver is not mutated.

        """
        direction = self.direction
        tolerance = self.tolerance
        structure = copy.deepcopy(self)
        for site in structure:
            if (
                lower_limit - tolerance
                <= site.frac_coords[direction]
                <= upper_limit + tolerance
            ):
                site.properties = {"selective_dynamics": [False, False, False]}
            else:
                site.properties = {"selective_dynamics": [True, True, True]}
        return structure
