"""
Surface energy calculation module for phase diagram construction.

This module provides the SurfaceEnergy class for calculating temperature and
voltage-dependent surface Gibbs free energies from DFT data.
"""

import numpy as np

from surface_pd.plot._reference_energies import (
    _get_reference_energies,
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
        functional: Complete DFT protocol label: ``"PBE+U"``,
            ``"SCAN+rVV10+U"``, or ``"r2SCAN+rVV10+U"``. The selected
            transition metal must have a confirmed bulk reference for it.

    Raises
    ------
    ValueError
        If the functional or transition-metal/functional combination has no
        confirmed reference-energy set.
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
        self._li_energy, self._o2_energy, self._bulk_energy = (
            _get_reference_energies(TM_species, functional)
        )

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
            self._o2_energy + delta_mu_kj_per_mol / _KJ_PER_MOL_PER_EV
        )

    def get_gibbs_free_energy(self):
        """
        Calculate surface Gibbs free energy.

        Returns
        -------
            Surface Gibbs free energy.
        """
        A = np.sin(self.gamma * np.pi / 180) * self.a * self.b
        G = (
            self.dft_energy
            + (self.nTM - self.nLi) * (self._li_energy - self.V)
            + (2 * self.nTM - self.nO) * self.g_oxygen()
            - self.nTM * self._bulk_energy
        ) / (2 * A)
        return G
