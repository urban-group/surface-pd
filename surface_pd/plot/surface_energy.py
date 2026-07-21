"""
Surface energy calculation module for phase diagram construction.

This module provides the SurfaceEnergy class for calculating temperature and
voltage-dependent surface Gibbs free energies from DFT data.
"""

import numpy as np

from surface_pd.plot.reference_energies import ReferenceEnergies
from surface_pd.thermodynamics import (
    FixedPressureOxygenChemicalPotential,
    ThermodynamicState,
)


class SurfaceEnergy:
    """Calculate the Gibbs free energy of one surface configuration.

    Parameters
    ----------
    V : float or array-like
        Electrochemical potential in volts.
    T : float or array-like
        Absolute temperature in kelvin. Must be positive and broadcastable
        with ``V``.
    nLi : float
        Number of lithium-like atoms in the slab model.
    nTM : float
        Number of transition-metal atoms in the slab model.
    nO : float
        Number of oxygen-like atoms in the slab model.
    dft_energy : float
        DFT total energy of the slab model in eV.
    a : float
        Positive first in-plane lattice length in angstroms.
    b : float
        Positive second in-plane lattice length in angstroms.
    gamma : float
        Angle between ``a`` and ``b`` in degrees, strictly between 0 and 180.
    reference_energies : ReferenceEnergies
        Validated Li, O2, and bulk LiTMO2 reference energies.

    Raises
    ------
    TypeError
        If ``reference_energies`` is not a :class:`ReferenceEnergies` object.

    Notes
    -----
    The constructor stores inputs without validating geometric or atom-count
    domains. :meth:`g_oxygen` validates temperature when energy is evaluated.
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
        reference_energies,
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
        if not isinstance(reference_energies, ReferenceEnergies):
            raise TypeError(
                "reference_energies must be a ReferenceEnergies object"
            )
        self.reference_energies = reference_energies

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
            same shape as ``T``; scalar input produces a NumPy scalar.

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
        model = FixedPressureOxygenChemicalPotential(
            raw_o2_energy_ev_per_molecule=(
                self.reference_energies.o2_raw_ev_per_molecule
            ),
            correction_ev_per_molecule=(
                self.reference_energies.o2_correction_ev_per_molecule
            ),
        )
        return model.evaluate(ThermodynamicState({"temperature": self.T}))

    def get_gibbs_free_energy(self):
        """Calculate the surface Gibbs free energy.

        Returns
        -------
        numpy.ndarray or numpy.float64
            Surface Gibbs free energy in eV per square angstrom. The result
            follows the broadcast shape of ``V`` and ``T``.

        Notes
        -----
        The energy is divided by twice the in-plane area, corresponding to the
        two equivalent surfaces of a symmetric slab model.
        """
        A = np.sin(self.gamma * np.pi / 180) * self.a * self.b
        G = (
            self.dft_energy
            + (self.nTM - self.nLi)
            * (self.reference_energies.li_ev_per_atom - self.V)
            + (2 * self.nTM - self.nO) * self.g_oxygen()
            - self.nTM
            * self.reference_energies.bulk_litmo2_ev_per_formula_unit
        ) / (2 * A)
        return G
