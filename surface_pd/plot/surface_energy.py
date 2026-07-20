"""
Surface energy calculation module for phase diagram construction.

This module provides the SurfaceEnergy class for calculating temperature and
voltage-dependent surface Gibbs free energies from DFT data.
"""

import numpy as np

from surface_pd.plot._reference_energies import (
    _BULK_ENERGY_BY_FUNCTIONAL,
    _LI_ENERGY_BY_FUNCTIONAL,
    _O2_ENERGY_BY_FUNCTIONAL,
)


class SurfaceEnergy:
    """
    Surface energy calculator for phase diagram construction.

    This class calculates the Gibbs free energy of a surface configuration
    as a function of electrochemical potential and temperature.

    Args:
        V: Electrochemical potential (voltage).
        T: Temperature in Kelvin.
        nLi: Number of lithium atoms.
        nTM: Number of transition metal atoms.
        nO: Number of oxygen atoms.
        dft_energy: DFT-calculated total energy.
        a: Lattice parameter a.
        b: Lattice parameter b.
        gamma: Angle gamma in degrees.
        TM_species: Transition metal species symbol.
        functional: DFT functional used (e.g., "PBE" or "SCAN").
    """

    def __init__(
        self,
        V,
        T,
        nLi,
        nTM,
        nO,
        dft_energy,
        a,
        b,
        gamma,
        TM_species,
        functional,
    ):
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
        """Return the reference lithium energy."""
        return -2.3

    def g_oxygen(self):
        """
        Calculate the temperature-dependent oxygen chemical potential.

        Equation to calculate the temperature dependent oxygen chemical
        potential.

        Returns
        -------
            Oxygen chemical potential.
        """
        E_O2 = _O2_ENERGY_BY_FUNCTIONAL[self.functional]
        T0 = 298  # K
        kB = 0.008314463  # kJ/(mol K)
        H0 = 8.683  # kJ/mol
        S0 = 205.147 * 1e-3  # kJ/(mol K)
        Cp = 3.5 * kB  # kJ/(mol K)
        delta_H = Cp * (self.T - T0)  # kJ/mol
        delta_S = Cp * np.log(self.T / T0)  # kJ/(mol K)
        delta_mu = (H0 + delta_H) - self.T * (S0 + delta_S)  # kJ/mol
        # print(delta_mu/96.487)
        g_O2 = (E_O2 + delta_mu) / 96.487  # eV/O2
        return g_O2 / 2

    def get_gibbs_free_energy(self):
        """
        Calculate surface Gibbs free energy.

        Returns
        -------
            Surface Gibbs free energy.
        """
        E_bulk = _BULK_ENERGY_BY_FUNCTIONAL[self.TM_species][self.functional]
        E_Li = _LI_ENERGY_BY_FUNCTIONAL[self.functional]
        A = np.sin(self.gamma * np.pi / 180) * self.a * self.b
        G = (
            self.dft_energy
            + (self.nTM - self.nLi) * (E_Li - self.V)
            + (2 * self.nTM - self.nO) * self.g_oxygen()
            - self.nTM * E_bulk
        ) / (2 * A)
        return G
