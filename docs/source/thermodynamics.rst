==============
Thermodynamics
==============

The :mod:`surface_pd.thermodynamics` package represents named thermodynamic
states and chemical-potential models independently of phase data and plotting.
All microscopic energies and chemical potentials use electronvolts (eV),
temperature uses kelvin (K), and voltage uses volts (V).

Model interface
===============

``ChemicalPotentialModel`` is a structural interface: a model declares the
state variables it needs and provides an ``evaluate(state)`` method. A custom
Python class may satisfy this interface without inheriting from a package base
class. Serialized configuration remains restricted to models implemented by
the package.

.. autoclass:: surface_pd.thermodynamics.ChemicalPotentialModel
    :members:

Thermodynamic state
===================

``ThermodynamicState`` holds finite, named scalar or array variables. Values
are defensively copied, broadcast to a common shape, and exposed as read-only
arrays. Models return chemical-potential arrays with that same shape.

.. autoclass:: surface_pd.thermodynamics.ThermodynamicState
    :members:

Constant chemical potential
===========================

A fixed reservoir has a state-independent chemical potential

.. math::

    \mu(\mathbf{x}) = \mu_0.

.. autoclass:: surface_pd.thermodynamics.ConstantChemicalPotential
    :members:

Direct chemical-potential coordinate
====================================

A direct coordinate varies a chemical potential relative to a configured
reference:

.. math::

    \mu = \mu_\mathrm{ref} + \Delta\mu.

A zero reference represents an absolute chemical-potential coordinate. A
nonzero reference gives the conventional chemical-potential deviation from a
named state.

.. autoclass:: surface_pd.thermodynamics.DirectChemicalPotential
    :members:

Intercalation-electrode voltage
===============================

For an intercalant transferring :math:`z` electrons per component, the
voltage-dependent chemical potential is

.. math::

    \mu(V) = \mu_\mathrm{ref} - zV.

The voltage-reference label is retained as scientific provenance and is not
interpreted as a lookup key.

.. autoclass:: surface_pd.thermodynamics.IntercalationChemicalPotential
    :members:

Fixed-pressure oxygen approximation
===================================

The oxygen model preserves the analytical approximation used by the legacy
Li-TM-O implementation. It returns the chemical potential per O atom:

.. math::

    \mu_\mathrm{O}(T) = \frac{1}{2}\left[
    E_\mathrm{O_2}^\mathrm{DFT} + E_\mathrm{corr}
    + \frac{\Delta\mu_\mathrm{O_2}^\mathrm{thermal}(T)}
    {96.487\ \mathrm{kJ\ mol^{-1}\ eV^{-1}}}\right],

where

.. math::

    \Delta\mu_\mathrm{O_2}^\mathrm{thermal}(T) =
    \Delta H^\circ + C_p(T-T_0)
    - T\left[S^\circ + C_p\ln(T/T_0)\right].

The constants are :math:`T_0=298.15` K,
:math:`\Delta H^\circ=8.683` kJ/mol O2,
:math:`S^\circ=205.147\times10^{-3}` kJ/(mol O2 K), and
:math:`C_p=3.5R` with
:math:`R=0.008314463` kJ/(mol K). The raw DFT energy and correction are supplied
in eV per O2 molecule.

This is a constant-heat-capacity approximation anchored by NIST-JANAF values,
not interpolation of a JANAF table. Pressure dependence is neglected. The
implementation is retained explicitly so the generalized API reproduces the
legacy package's corrected numerical results.

.. autoclass:: surface_pd.thermodynamics.FixedPressureOxygenChemicalPotential
    :members:
