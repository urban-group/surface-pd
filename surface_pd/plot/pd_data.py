import numpy as np
import pandas as pd

from surface_pd.plot.surface_energy import SurfaceEnergy
from surface_pd.error import InvalidPhasesAlignError

E_bulk_by_funtional = {
    'Ni': {"PBE+U": -19.92283375,
           "SCAN+rVV10+U": -36.8133525},
    'Co':
        {"PBE+U": -22.69242,
         "SCAN+rVV10+U": -37.2001966667,
         "r2SCAN+rVV10+U": -32.5698933333},
    'Mn':
        {"PBE+U": -26.319605,
         "SCAN+rVV10+U": 0,
         "r2SCAN+rVV10+U": -36.4717}
    }


class PdData(object):
    """

    Args:
        dataframe:
        lithium_like_species:
        oxygen_like_species:
        functional:

    """
    def __init__(self,
                 dataframe: pd.DataFrame,
                 lithium_like_species: str,
                 oxygen_like_species: str,
                 functional: str):
        self.dataframe = dataframe
        self.lithium_like_species = lithium_like_species
        self.oxygen_like_species = oxygen_like_species
        self.functional = functional

    def tm_species(self):
        """
        Find the transition metal species in the material

        Returns:
            transition metal species
        """

        # Define transition metal from periodic table
        TM_dataset = ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']
        column_names = self.dataframe.columns
        for name in column_names:
            if name in TM_dataset:
                return name

    def standardize_pd_data(self):
        """
        Normalized the pandas dataframe based on the number of TM species.
        This is to make the DFT energies directly comparable.

        Returns:

        """
        self.dataframe.reset_index(inplace=True)
        max_TM_species = max(self.dataframe.iloc[:][self.tm_species()])
        for i, num in enumerate(self.dataframe.iloc[:][self.tm_species()]):
            if num < max_TM_species:
                multiple = max_TM_species / num
                self.dataframe.at[i, self.lithium_like_species] = \
                    self.dataframe.iloc[i][self.lithium_like_species] \
                    * multiple
                self.dataframe.at[i, self.tm_species()] = \
                    self.dataframe.iloc[i][self.tm_species()] * multiple
                self.dataframe.at[i, self.oxygen_like_species] = \
                    self.dataframe.iloc[i][self.oxygen_like_species] * multiple
                self.dataframe.at[i, "E"] = \
                    self.dataframe.iloc[i]["E"] * multiple
                self.dataframe.at[i, "a"] = self.dataframe.iloc[i]["a"]
                self.dataframe.at[i, "b"] = \
                    self.dataframe.iloc[i]["b"] * multiple
                self.dataframe.at[i, "gamma"] = self.dataframe.iloc[i]["gamma"]
            else:
                pass
            self.dataframe = self.dataframe.sort_values(
                by=[self.oxygen_like_species,
                    self.lithium_like_species],
                ascending=[False, False])
            self.dataframe = self.dataframe.reset_index(drop=True)
        return self

    def get_check_phases(self):
        """

        Returns:

        """
        max_Li = max(self.dataframe.iloc[:][self.lithium_like_species])
        max_O = max(self.dataframe.iloc[:][self.oxygen_like_species])
        if max_Li == max_O / 2:
            min_Li = min(self.dataframe.iloc[:][self.lithium_like_species])
            min_O = min(self.dataframe.iloc[:][self.oxygen_like_species])
            phases_set = self.dataframe.loc[
                (self.dataframe[self.lithium_like_species] == min_Li)
                & (self.dataframe[self.oxygen_like_species] == min_O)]
            return phases_set
        else:
            phases_set = self.dataframe.loc[
                (self.dataframe[self.lithium_like_species] == max_Li)
                & (self.dataframe[self.oxygen_like_species] == max_O)]
            return phases_set

    def get_the_shift_energy(self,
                             checked_phases):
        """
        Shift energy is calculated using the slab models with different
        number of layers. Though they have different number of layers,
        but proper normalization will make them same. This energy will be
        used to eliminate the energy difference between these models.

        Args:
            checked_phases:

        Returns:
            shift energy
        """
        num_Li, num_O, E = [], [], []
        for phase in checked_phases:
            num_Li.append(int(phase.iloc[:][self.lithium_like_species]))
            num_O.append(int(phase.iloc[:][self.oxygen_like_species]))
            E.append(float(phase.iloc[:]["E"]))

        # Check number of extra bulk structures in the structure model
        num_bulk = num_Li[1] - num_Li[0]
        if num_bulk != (num_O[1] - num_O[0]) / 2:
            raise InvalidPhasesAlignError

        # Calculate the shift energy
        [a, b, gamma] = [checked_phases[0]["a"],
                         checked_phases[0]["b"],
                         checked_phases[0]["gamma"]]
        # Calculate the slab surface area
        A = float(np.sin(gamma * np.pi / 180) * a * b)
        E_bulk = E_bulk_by_funtional[self.tm_species()][self.functional]
        E_shift = (1 / (2 * A)) * (E[0] - E[1] + num_bulk * E_bulk)
        if not all(E):
            E_shift = 0
        print('Surface area:', A, '\n',
              'E_bulk:', E_bulk, '\n',
              'E_shift:', E_shift)
        return E_shift

    def get_surface_energy(self,
                           V: np.ndarray,
                           T: np.ndarray):
        """
        Using the user defined voltage and temperature meshes to calculate
        the surface free energy.

        Args:
            V: voltage mesh
            T: temperature mesh

        Returns:
            list of surface energies
        """
        E = []
        for i in range(len(self.dataframe)):
            E = np.append(
                E,
                SurfaceEnergy(
                    V=V,
                    T=T,
                    nLi=self.dataframe.iloc[i][self.lithium_like_species],
                    nTM=self.dataframe.iloc[i][self.tm_species()],
                    nO=self.dataframe.iloc[i][self.oxygen_like_species],
                    dft_energy=self.dataframe.iloc[i]['E'],
                    a=self.dataframe.iloc[i]['a'],
                    b=self.dataframe.iloc[i]['b'],
                    gamma=self.dataframe.iloc[i]['gamma'],
                    TM_species=self.tm_species(),
                    functional=self.functional
                ).gibbs_free_energy()
            )
        return E
