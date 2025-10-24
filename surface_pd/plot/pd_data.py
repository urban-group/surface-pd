"""
Phase diagram data management module.

This module provides the PdData class for loading, processing, and normalizing
DFT energy data for surface phase diagram construction.
"""

import numpy as np
import pandas as pd

from surface_pd.plot.surface_energy import SurfaceEnergy

# Constants for bulk energies by functional
E_O2_by_funtional = {"PBE": -9.86, "SCAN": -10.45}
E_bulk_by_funtional = {
    "Ni": {"PBE": -5.55, "SCAN": -5.77},
    "Co": {"PBE": -7.11, "SCAN": -7.39},
    "Mn": {"PBE": -9.00, "SCAN": -9.20},
}
E_Li_by_funtional = {"PBE": -1.90, "SCAN": -2.00}


class PdData:
    """
    Phase diagram data container and processor.

    Args:
        dataframe:
        lithium_like_species:
        oxygen_like_species:
        functional:

    """

    def __init__(
        self,
        dataframe: pd.DataFrame,
        lithium_like_species: str,
        oxygen_like_species: str,
        functional: str,
    ):
        self.dataframe = dataframe
        self.lithium_like_species = lithium_like_species
        self.oxygen_like_species = oxygen_like_species
        self.functional = functional

    def tm_species(self):
        """
        Find the transition metal species in the material.

        Returns
        -------
            transition metal species
        """
        # Define transition metal from periodic table
        TM_dataset = ["Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"]
        column_names = self.dataframe.columns
        for name in column_names:
            if name in TM_dataset:
                return name

    def standardize_pd_data(self):
        """
        Normalize the pandas dataframe based on the number of TM species.

        This is to make the DFT energies directly comparable.

        Returns
        -------

        """
        self.dataframe.reset_index(inplace=True)

        # Standardize column names for compatibility
        if 'E' in self.dataframe.columns and 'dft_energy' not in self.dataframe.columns:
            self.dataframe.rename(columns={'E': 'dft_energy'}, inplace=True)

        max_TM_species = max(self.dataframe.iloc[:][self.tm_species()])
        for i, num in enumerate(self.dataframe.iloc[:][self.tm_species()]):
            if num < max_TM_species:
                multiple = max_TM_species / num
                self.dataframe.loc[i, self.lithium_like_species] = (
                    self.dataframe.loc[i, self.lithium_like_species] * multiple
                )
                self.dataframe.loc[i, self.oxygen_like_species] = (
                    self.dataframe.loc[i, self.oxygen_like_species] * multiple
                )
                self.dataframe.loc[i, "dft_energy"] = (
                    self.dataframe.loc[i, "dft_energy"] * multiple
                )
                self.dataframe.loc[i, "a"] = (
                    self.dataframe.loc[i, "a"] * multiple
                )
                self.dataframe.loc[i, "b"] = (
                    self.dataframe.loc[i, "b"] * multiple
                )
        return self.dataframe

    def get_check_phases(self):
        """
        Get phases to check for alignment energy calculation.

        Returns
        -------

        """
        max_Li = max(self.dataframe.iloc[:][self.lithium_like_species])
        max_O = max(self.dataframe.iloc[:][self.oxygen_like_species])
        if max_Li == max_O / 2:
            min_Li = min(self.dataframe.iloc[:][self.lithium_like_species])
            min_O = min(self.dataframe.iloc[:][self.oxygen_like_species])
            phases_set = self.dataframe.loc[
                (self.dataframe[self.lithium_like_species] == min_Li)
                & (self.dataframe[self.oxygen_like_species] == min_O)
            ]
            phases_list = []
            for i, phase in phases_set.iterrows():
                phases_list.append(phase)
            return phases_list

    def get_the_shift_energy(self, checked_phases):
        """
        Calculate shift energy between different slab thicknesses.

        Shift energy is calculated using the slab models with different
        number of layers. Though they have different number of layers,
        but proper normalization will make them same. This energy will be
        used to eliminate the energy difference between these models.

        Args:
            checked_phases:

        Returns
        -------
            shift energy
        """
        num_Li, num_O, E = [], [], []
        for phase in checked_phases:
            num_Li.append(int(phase.iloc[:][self.lithium_like_species]))
            num_O.append(int(phase.iloc[:][self.oxygen_like_species]))
            E.append(float(phase.iloc[:]["dft_energy"]))

        if len(checked_phases) > 1:
            num_bulk = num_Li[0] - num_Li[1] + (num_O[0] - num_O[1]) / 2
            a = checked_phases[0]["a"]
            b = checked_phases[0]["b"]
            gamma = checked_phases[0]["gamma"]
        # Calculate the slab surface area
        A = float(np.sin(gamma * np.pi / 180) * a * b)
        E_bulk = E_bulk_by_funtional[self.tm_species()][self.functional]
        E_shift = (1 / (2 * A)) * (E[0] - E[1] + num_bulk * E_bulk)
        if not all(E):
            E_shift = 0
        print(
            "Surface area:",
            A,
            "\n",
            "E_bulk:",
            E_bulk,
            "\n",
            "Shift energy:",
            E_shift,
        )
        return E_shift

    def get_surface_energy(self, V: np.ndarray, T: np.ndarray):
        """
        Calculate surface free energy over voltage and temperature meshes.

        Using the user defined voltage and temperature meshes to calculate
        the surface free energy.

        Args:
            V: voltage mesh
            T: temperature mesh

        Returns
        -------
            list of surface energies
        """
        E = []
        for i in range(len(self.dataframe)):
            E = np.append(
                E,
                SurfaceEnergy(
                    V=V,
                    T=T,
                    nLi=self.dataframe.loc[i, self.lithium_like_species],
                    nTM=self.dataframe.loc[i, self.tm_species()],
                    nO=self.dataframe.loc[i, self.oxygen_like_species],
                    dft_energy=self.dataframe.loc[i, "dft_energy"],
                    a=self.dataframe.loc[i, "a"],
                    b=self.dataframe.loc[i, "b"],
                    gamma=self.dataframe.loc[i, "gamma"],
                    TM_species=self.tm_species(),
                    functional=self.functional,
                ).get_gibbs_free_energy(),
            )
        return E
