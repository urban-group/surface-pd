==============
Thermodynamics
==============

The :mod:`surface_pd.thermodynamics` package represents named thermodynamic
states and chemical-potential models independently of phase data and plotting.
All microscopic energies and chemical potentials use electronvolts (eV),
temperature uses kelvin (K), and voltage uses volts (V).

Phase and reference data
========================

``Phase`` stores the absolute elemental counts :math:`n_{s,i}`, total DFT
energy :math:`E_s`, surface-unit-cell area :math:`A_s`, and number of
equivalent surfaces :math:`m_s` for a calculated phase :math:`s`. Counts are
not reduced to a formula unit. This retains enough information to compare
cells containing different numbers of atoms and to support adsorption as well
as substitutional defects.

For chemical potentials :math:`\mu_i`, a later grand-potential evaluator uses

.. math::

    \Omega_s = E_s - \sum_i n_{s,i}\mu_i,

and compares phases using the surface-normalized value

.. math::

    \gamma_s = \frac{\Omega_s}{m_s A_s}.

``surface_area_angstrom2`` is the area of one surface unit cell, while
``surface_multiplicity`` records how many equivalent surfaces the calculated
cell represents. Keeping the two quantities explicit avoids hiding the
factor-of-two convention commonly used for symmetric slabs.

.. autoclass:: surface_pd.thermodynamics.Phase
    :members:

``PhaseDataset`` groups phases calculated with one declared method and an
explicit ordered component basis. Phase identifiers are local to a dataset;
the canonical qualified identifier is ``dataset_id:phase_id``. The dataset
does not read files or alter the supplied energies.

.. autoclass:: surface_pd.thermodynamics.PhaseDataset
    :members:

``ReferencePhase`` represents one bulk thermodynamic equality. For reference
composition :math:`n_{r,i}` and energy per formula unit :math:`G_r`, the
constraint is

.. math::

    \sum_i n_{r,i}\mu_i = G_r.

The present object records this scientific input without choosing dependent
chemical potentials or solving the equality; that responsibility belongs to
the generalized grand-potential model.

.. autoclass:: surface_pd.thermodynamics.ReferencePhase
    :members:

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
