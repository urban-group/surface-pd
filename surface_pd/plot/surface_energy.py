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

_STANDARD_TEMPERATURE_K = 298.15
_GAS_CONSTANT_KJ_PER_MOL_K = 0.008314463
_O2_ENTHALPY_INCREMENT_KJ_PER_MOL = 8.683
_O2_STANDARD_ENTROPY_KJ_PER_MOL_K = 205.147e-3
_O2_HEAT_CAPACITY_KJ_PER_MOL_K = 3.5 * _GAS_CONSTANT_KJ_PER_MOL_K
_KJ_PER_MOL_PER_EV = 96.487


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
        Return the temperature-dependent oxygen chemical potential.

        The DFT reference is combined with ideal-gas enthalpy and entropy
        corrections at standard pressure. Molar thermal corrections are
        converted from kJ/mol to eV before being added to the DFT energy.

        Returns
        -------
        numpy.ndarray or numpy.float64
            Oxygen chemical potential in eV per O atom. The result has the
            same shape as ``T``.

        Raises
        ------
        ValueError
            If any absolute temperature is zero or negative.

        Notes
        -----
        Pressure dependence is neglected. The reference temperature is
        298.15 K, and the standard-state O2 entropy and 0-to-298.15 K
        enthalpy increment are taken from NIST-JANAF thermochemical data.
        """
        E_O2 = _O2_ENERGY_BY_FUNCTIONAL[self.functional]
        temperature = np.asarray(self.T)
        if np.any(temperature <= 0):
            raise ValueError("Temperature must be positive in Kelvin.")

        delta_H = _O2_HEAT_CAPACITY_KJ_PER_MOL_K * (
            temperature - _STANDARD_TEMPERATURE_K
        )
        delta_S = _O2_HEAT_CAPACITY_KJ_PER_MOL_K * np.log(
            temperature / _STANDARD_TEMPERATURE_K
        )
        delta_mu_kj_per_mol = (
            _O2_ENTHALPY_INCREMENT_KJ_PER_MOL + delta_H
        ) - temperature * (_O2_STANDARD_ENTROPY_KJ_PER_MOL_K + delta_S)
        return 0.5 * (
            E_O2 + delta_mu_kj_per_mol / _KJ_PER_MOL_PER_EV
        )

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
