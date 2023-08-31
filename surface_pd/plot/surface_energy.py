import numpy as np

# DFT energies for pure BCC Li, calculated using different functionals
E_Li_by_funtional = {"PBE+U": -1.89965,
                     "SCAN+rVV10+U": -2.33333,
                     "r2SCAN+rVV10+U": -2.32338}

# DFT energies for isolated O2 molecule, calculated using different functionals
# oxygen correction: -1.36
E_O2_by_funtional = {"PBE+U": -9.86018 + 1.36,
                     "SCAN+rVV10+U": -12.00701,
                     "r2SCAN+rVV10+U": -11.54833}

# DFT energies for different LiTMO2, calculated using different functionals
E_bulk_by_funtional = {
    'Ni':
        {"PBE+U": -19.92283375,
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


class SurfaceEnergy(object):
    """
    Surface free energy class.

    Args:
        V: Voltage.
        T: Temperature.
        nLi: Number of Li atoms.
        nTM: Number of transition metal atoms.
        nO: Number of O atoms.
        dft_energy: DFT energy.
        a: Lattice parameter a, used to calculate the surface area.
        b: Lattice parameter b, used to calculate the surface area.
        gamma: Lattice parameter gamma, used to calculate the surface area.
        TM_species: Transition metal species in the material.
        functional: Functional used to perform the calculations.

    """
    def __init__(self,
                 V, T, nLi, nTM, nO, dft_energy,
                 a, b, gamma, TM_species, functional):
        self.V = V
        self.T = T
        self.nLi = nLi
        self.nTM = nTM
        self.nO = nO
        self.dft_energy = dft_energy
        self.a = a
        self.b = b
        self.gamma = gamma
        self.TM_species = TM_species
        self.functional = functional

    @property
    def e_li(self):
        return -2.3

    @e_li.setter
    def e_li(self, a):
        self.e_li = a

    def g_oxygen(self):
        """
        Equation to calculate the temperature dependent oxygen chemical
        potential.

        Returns:
            Oxygen chemical potential.
        """
        E_O2 = E_O2_by_funtional[self.functional]
        T0 = 298  # K
        kB = 0.008314463  # kJ/(mol K)
        H0 = 8.683  # kJ/mol
        S0 = 205.147 * 1e-3  # kJ/(mol K)
        Cp = 3.5 * kB  # kJ/(mol K)
        delta_H = Cp * (self.T - T0)  # kJ/mol
        delta_S = Cp * np.log(self.T / T0)  # kJ/(mol K)
        delta_mu = (H0 + delta_H) - self.T * (S0 + delta_S)  # kJ/mol
        # print(delta_mu/96.487)
        mu_oxygen = (1 / 2) * (E_O2 + (delta_mu / 96.487))
        return mu_oxygen

    def gibbs_free_energy(self):
        """
        Equation to calculate the surface Gibbs free energy.

        Returns:
            Surface Gibbs free energy.
        """
        E_bulk = E_bulk_by_funtional[self.TM_species][self.functional]
        E_Li = E_Li_by_funtional[self.functional]
        A = np.sin(self.gamma * np.pi / 180) * self.a * self.b
        G = (self.dft_energy + (self.nTM - self.nLi) * (E_Li - self.V)
             + (2 * self.nTM - self.nO) *
             self.g_oxygen() - self.nTM * E_bulk) / (2 * A)
        return G
