"""Validate and prepare phase data for surface phase diagrams."""

import logging

import numpy as np
import pandas as pd
from pymatgen.core import Element

from surface_pd.plot.reference_energies import ReferenceEnergies
from surface_pd.plot.surface_energy import SurfaceEnergy

logger = logging.getLogger(__name__)

_AREA_RELATIVE_TOLERANCE = 5e-3
_AREA_ABSOLUTE_TOLERANCE_ANGSTROM2 = 1e-8
_STOICHIOMETRY_TOLERANCE = 1e-8


class PdData:
    """Own, normalize, and calculate one surface phase dataset.

    Parameters
    ----------
    dataframe : pandas.DataFrame
        Phase rows containing atom counts, DFT energies, and in-plane lattice
        parameters. A defensive copy is stored; the caller's dataframe is
        never mutated.
    lithium_like_species : str
        Column containing lithium-like species counts.
    oxygen_like_species : str
        Column containing oxygen-like species counts.
    reference_energies : ReferenceEnergies
        Validated user-provided Li, O2, and bulk LiTMO2 reference energies.

    Raises
    ------
    TypeError
        If ``dataframe`` is not a pandas dataframe or ``reference_energies``
        is not a :class:`ReferenceEnergies` object.
    ValueError
        If the lithium-like and oxygen-like column names are identical.

    Notes
    -----
    Call :meth:`standardize_pd_data` before alignment or energy methods.
    """

    def __init__(
        self,
        dataframe: pd.DataFrame,
        lithium_like_species: str,
        oxygen_like_species: str,
        reference_energies: ReferenceEnergies,
    ):
        if not isinstance(dataframe, pd.DataFrame):
            raise TypeError("dataframe must be a pandas DataFrame")
        if lithium_like_species == oxygen_like_species:
            raise ValueError(
                "lithium_like_species and oxygen_like_species must be distinct"
            )
        if not isinstance(reference_energies, ReferenceEnergies):
            raise TypeError(
                "reference_energies must be a ReferenceEnergies object"
            )
        self.dataframe = dataframe.copy(deep=True)
        self.lithium_like_species = lithium_like_species
        self.oxygen_like_species = oxygen_like_species
        self.reference_energies = reference_energies
        self._transition_metal = None
        self._is_standardized = False

    @property
    def transition_metal(self) -> str:
        """Return the detected transition-metal element symbol.

        Raises
        ------
        RuntimeError
            If the data has not yet been standardized.
        """
        self._require_standardized()
        return self._transition_metal

    @staticmethod
    def _is_transition_metal_symbol(value: object) -> bool:
        """Return whether *value* names a transition-metal element."""
        if not isinstance(value, str):
            return False
        try:
            return Element(value).is_transition_metal
        except ValueError:
            return False

    def _detect_transition_metal(self) -> str:
        """Return the sole transition-metal column or raise ``ValueError``."""
        detected = sorted(
            str(column)
            for column in self.dataframe.columns
            if self._is_transition_metal_symbol(column)
        )
        if len(detected) != 1:
            description = ", ".join(detected) if detected else "none"
            raise ValueError(
                "Phase data must contain exactly one transition-metal "
                f"column; detected: {description}."
            )
        return detected[0]

    def _normalize_energy_column(self) -> None:
        """Normalize ``E``/``dft_energy`` aliases without ambiguity."""
        has_short = "E" in self.dataframe
        has_canonical = "dft_energy" in self.dataframe
        if has_short and has_canonical:
            try:
                short_values = pd.to_numeric(
                    self.dataframe["E"], errors="raise"
                ).to_numpy(dtype=float)
                canonical_values = pd.to_numeric(
                    self.dataframe["dft_energy"], errors="raise"
                ).to_numpy(dtype=float)
            except (TypeError, ValueError) as error:
                raise ValueError(
                    "E and dft_energy must contain numeric values."
                ) from error
            if not np.array_equal(short_values, canonical_values):
                raise ValueError(
                    "Phase data contains conflicting E and dft_energy columns."
                )
            self.dataframe.drop(columns="E", inplace=True)
        elif has_short:
            self.dataframe.rename(columns={"E": "dft_energy"}, inplace=True)

    def _validate_columns_and_values(self, transition_metal: str) -> None:
        """Validate required columns and their physical numeric domains."""
        required = [
            self.lithium_like_species,
            transition_metal,
            self.oxygen_like_species,
            "dft_energy",
            "a",
            "b",
            "gamma",
        ]
        missing = [name for name in required if name not in self.dataframe]
        if missing:
            raise ValueError(
                "Missing required columns: " + ", ".join(missing) + "."
            )
        if self.dataframe.empty:
            raise ValueError("Phase data must contain at least one row.")

        for column in required:
            try:
                numeric = pd.to_numeric(self.dataframe[column], errors="raise")
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f"{column} must contain only numeric values"
                ) from error
            if not np.isfinite(numeric.to_numpy(dtype=float)).all():
                raise ValueError(f"{column} must contain only finite values")
            self.dataframe[column] = numeric.astype(float)

        if (self.dataframe[self.lithium_like_species] < 0).any():
            raise ValueError(
                f"{self.lithium_like_species} counts must be nonnegative"
            )
        if (self.dataframe[self.oxygen_like_species] < 0).any():
            raise ValueError(
                f"{self.oxygen_like_species} counts must be nonnegative"
            )
        if (self.dataframe[transition_metal] <= 0).any():
            raise ValueError(f"{transition_metal} counts must be positive")
        for lattice_length in ("a", "b"):
            if (self.dataframe[lattice_length] <= 0).any():
                raise ValueError(f"{lattice_length} must be positive")
        if (
            (self.dataframe["gamma"] <= 0)
            | (self.dataframe["gamma"] >= 180)
        ).any():
            raise ValueError("gamma must be between 0 and 180 degrees")

    def standardize_pd_data(self) -> "PdData":
        """Validate and normalize phase rows to their maximum TM count.

        Atom counts, DFT energy, and surface area are extensive quantities and
        are multiplied by the same row-specific factor. Area is scaled by
        changing ``b`` while keeping ``a`` and ``gamma`` fixed. Rows are sorted
        by oxygen-like and lithium-like counts in descending order.

        Returns
        -------
        PdData
            This instance, enabling explicit fluent use. Repeated calls are
            idempotent.
        """
        if self._is_standardized:
            return self

        self._normalize_energy_column()
        transition_metal = self._detect_transition_metal()
        self._validate_columns_and_values(transition_metal)

        if not isinstance(self.dataframe.index, pd.RangeIndex) or (
            self.dataframe.index.name is not None
        ):
            index_name = self.dataframe.index.name or "index"
            if index_name in self.dataframe:
                index_name = "source_index"
            self.dataframe.index = self.dataframe.index.rename(index_name)
            self.dataframe.reset_index(inplace=True)

        maximum_tm_count = self.dataframe[transition_metal].max()
        factors = maximum_tm_count / self.dataframe[transition_metal]
        extensive_columns = [
            self.lithium_like_species,
            transition_metal,
            self.oxygen_like_species,
            "dft_energy",
            "b",
        ]
        self.dataframe.loc[:, extensive_columns] = self.dataframe[
            extensive_columns
        ].mul(factors, axis=0)
        self.dataframe.sort_values(
            by=[self.oxygen_like_species, self.lithium_like_species],
            ascending=[False, False],
            kind="stable",
            inplace=True,
        )
        self.dataframe.reset_index(drop=True, inplace=True)
        self._transition_metal = transition_metal
        self._is_standardized = True
        return self

    def _require_standardized(self) -> None:
        """Raise if normalized-data methods are called too early."""
        if not self._is_standardized:
            raise RuntimeError(
                "Call standardize_pd_data() before this operation."
            )

    def get_alignment_reference(self) -> pd.Series:
        """Return one deterministic endpoint phase for dataset alignment.

        For a fully lithiated maximum composition, the minimum Li/O endpoint
        is selected; otherwise the maximum Li/O endpoint is selected. Energy
        breaks composition ties, with the lowest-energy row returned.

        Returns
        -------
        pandas.Series
            A copy of the selected standardized phase row.
        """
        self._require_standardized()
        lithium = self.dataframe[self.lithium_like_species]
        oxygen = self.dataframe[self.oxygen_like_species]
        if np.isclose(lithium.max(), oxygen.max() / 2):
            target_lithium = lithium.min()
            target_oxygen = oxygen.min()
        else:
            target_lithium = lithium.max()
            target_oxygen = oxygen.max()
        candidates = self.dataframe.loc[
            np.isclose(lithium, target_lithium)
            & np.isclose(oxygen, target_oxygen)
        ]
        selected_index = candidates["dft_energy"].idxmin()
        return self.dataframe.loc[selected_index].copy()

    @staticmethod
    def _surface_area(phase: pd.Series) -> float:
        """Return one phase's in-plane area in square angstroms."""
        return float(
            np.sin(np.deg2rad(phase["gamma"])) * phase["a"] * phase["b"]
        )

    def calculate_shift_energy(self, other: "PdData") -> float:
        """Return the energy-density shift to add to ``other``.

        The selected phases must differ by a whole number of LiTMO2 formula
        units and represent the same relaxed surface-cell area within 0.5%.

        Parameters
        ----------
        other : PdData
            The second standardized dataset whose surface energies will be
            shifted.

        Returns
        -------
        float
            Constant shift in eV per square angstrom.

        Raises
        ------
        TypeError
            If ``other`` is not :class:`PdData`.
        RuntimeError
            If either dataset has not been standardized.
        ValueError
            If reference energies, transition metals, surface areas,
            stoichiometry, or zero-energy masks are incompatible.
        """
        if not isinstance(other, PdData):
            raise TypeError("other must be a PdData object")
        self._require_standardized()
        other._require_standardized()
        if self.reference_energies != other.reference_energies:
            raise ValueError("Datasets use different reference energies.")
        if self.transition_metal != other.transition_metal:
            raise ValueError("Datasets use different transition metals.")

        first = self.get_alignment_reference()
        second = other.get_alignment_reference()
        first_area = self._surface_area(first)
        second_area = self._surface_area(second)
        if not np.isclose(
            first_area,
            second_area,
            rtol=_AREA_RELATIVE_TOLERANCE,
            atol=_AREA_ABSOLUTE_TOLERANCE_ANGSTROM2,
        ):
            raise ValueError(
                "Alignment phases have incompatible surface areas: "
                f"{first_area} and {second_area} angstrom^2."
            )

        differences = np.array(
            [
                second[self.lithium_like_species]
                - first[self.lithium_like_species],
                second[self.transition_metal] - first[self.transition_metal],
                (
                    second[self.oxygen_like_species]
                    - first[self.oxygen_like_species]
                )
                / 2,
            ],
            dtype=float,
        )
        number_of_bulk_units = round(differences[0])
        if not np.allclose(
            differences,
            number_of_bulk_units,
            rtol=0,
            atol=_STOICHIOMETRY_TOLERANCE,
        ):
            raise ValueError(
                "Alignment phases must differ by whole LiTMO2 formula units."
            )

        first_energy = float(first["dft_energy"])
        second_energy = float(second["dft_energy"])
        if first_energy == 0 and second_energy == 0:
            return 0.0
        if (first_energy == 0) != (second_energy == 0):
            raise ValueError(
                "Cannot align datasets when only one alignment energy is zero."
            )
        bulk_energy = (
            self.reference_energies.bulk_litmo2_ev_per_formula_unit
        )
        shift = (
            first_energy
            - second_energy
            + number_of_bulk_units * bulk_energy
        ) / (2 * first_area)
        logger.info(
            "Surface area: %s\nBulk units: %s\nShift energy: %s",
            first_area,
            number_of_bulk_units,
            shift,
        )
        return float(shift)

    def get_surface_energy(
        self,
        V: np.ndarray,
        T: np.ndarray,
    ) -> np.ndarray:
        """Calculate surface energies over matching voltage/temperature meshes.

        Parameters
        ----------
        V : array-like
            Voltage values in volts.
        T : array-like
            Absolute temperatures in kelvin, with the same shape as ``V``.

        Returns
        -------
        numpy.ndarray
            Surface energies in eV per square angstrom with shape
            ``(number_of_phases, *V.shape)``.

        Raises
        ------
        RuntimeError
            If the dataset has not been standardized.
        ValueError
            If ``V`` and ``T`` do not have identical shapes.
        """
        self._require_standardized()
        voltage = np.asarray(V)
        temperature = np.asarray(T)
        if voltage.shape != temperature.shape:
            raise ValueError("V and T must have identical shapes.")

        energies = []
        for _, phase in self.dataframe.iterrows():
            energies.append(
                SurfaceEnergy(
                    V=voltage,
                    T=temperature,
                    nLi=phase[self.lithium_like_species],
                    nTM=phase[self.transition_metal],
                    nO=phase[self.oxygen_like_species],
                    dft_energy=phase["dft_energy"],
                    a=phase["a"],
                    b=phase["b"],
                    gamma=phase["gamma"],
                    reference_energies=self.reference_energies,
                ).get_gibbs_free_energy()
            )
        return np.stack(energies, axis=0)
