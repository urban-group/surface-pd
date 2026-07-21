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

    \gamma_s = \frac{\Omega_s}{N_{s,\mathrm{surfaces}} A_s}.

``surface_area_angstrom2`` is the area of one surface unit cell, while
``number_of_surfaces`` records how many equivalent surfaces the calculated
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

Constrained grand-potential model
=================================

``GrandPotentialModel`` combines independent chemical-potential laws with
bulk reference equalities. The user chooses which components are independent;
all remaining components are dependent. The model never chooses this partition
automatically.

Let :math:`C` be the reference-composition matrix, with one row per
``ReferencePhase`` and columns in the model's declared component order. After
partitioning its columns into dependent and independent blocks, every reference
equality can be written as

.. math::

    C_D\boldsymbol{\mu}_D + C_I\boldsymbol{\mu}_I
    = \mathbf{g}_\mathrm{ref}.

The dependent potentials are therefore obtained from the exact linear solve

.. math::

    C_D\boldsymbol{\mu}_D
    = \mathbf{g}_\mathrm{ref} - C_I\boldsymbol{\mu}_I.

The initial implementation deliberately requires :math:`C_D` to be square and
full-rank. Underdetermined, overdetermined, and singular systems are rejected;
there is no implicit choice of reservoirs and no least-squares approximation.
All reference phases in one model must use identical calculation-method
provenance. A dataset evaluated with those references must use the same method.

For example, an ``AB2`` reference gives

.. math::

    \mu_A + 2\mu_B = G_\mathrm{AB2}.

If :math:`\mu_A` is supplied independently, the model solves the directly
auditable expression

.. math::

    \mu_B = \frac{G_\mathrm{AB2} - \mu_A}{2}.

.. autoclass:: surface_pd.thermodynamics.GrandPotentialModel
    :members:

Grand-potential results
=======================

For each phase :math:`s`, the model first calculates the total excess grand
potential of the supplied atomistic cell,

.. math::

    \Omega_s = E_s - \sum_i n_{s,i}\mu_i.

It then reports both supported normalizations:

.. math::

    \Omega_s^\mathrm{cell} = \frac{\Omega_s}{m_s},
    \qquad
    \gamma_s = \frac{\Omega_s}{m_s A_s}.

Here :math:`m_s` is the number of equivalent surfaces represented by the
atomistic cell and :math:`A_s` is the area of one surface unit cell. Thus
``grand_potential_ev_per_surface_cell`` is useful for adsorption and other
surface-cell processes, while
``surface_grand_potential_ev_per_angstrom2`` is the quantity used to compare
surface stability.

Each result tensor has shape ``(number_of_phases, *state_shape)``. Phase order
matches the dataset, and all arrays are read-only. The result retains every
phase; stable-phase selection and tie handling belong to the later diagram
evaluation layer.

.. autoclass:: surface_pd.thermodynamics.GrandPotentialResult
    :members:

.. _dataset-alignment:

Explicit dataset alignment
==========================

Slab datasets with different numbers of bulk-like layers can share the same
surface termination but use different total-energy zeros. ``DatasetAlignment``
places one target dataset on the energy convention of one root dataset using
explicitly named anchor phases and exactly one bulk ``ReferencePhase``. It
never guesses anchors from compositions or energies.

The signed integer number of bulk units :math:`k` is defined by

.. math::

    \mathbf{n}_\mathrm{target}
    - \mathbf{n}_\mathrm{reference}
    = k\mathbf{n}_\mathrm{bulk}.

The constant total-energy offset applied to the target dataset is

.. math::

    \Delta E_\mathrm{target}
    = E_\mathrm{reference} - E_\mathrm{target}
    + kG_\mathrm{bulk}.

Consequently, the aligned anchors obey

.. math::

    E_\mathrm{target} + \Delta E_\mathrm{target}
    = E_\mathrm{reference} + kG_\mathrm{bulk}.

The sign of :math:`k` follows target minus reference and may be positive,
negative, or zero. ``energy_offset_ev`` is a total cell energy in eV, not an
energy per area. ``GrandPotentialModel`` adds this offset before applying its
ordinary surface-cell and surface-area normalizations.

Alignment requires identical ordered component bases and exact
calculation-method provenance. The explicitly selected anchors must have equal
numbers of represented surfaces, and their surface areas must agree using a relative
tolerance of 0.5% and an absolute tolerance of
:math:`10^{-8}` square angstroms. Their composition difference must be exactly
an integer multiple of the bulk formula unit across every component.

.. autoclass:: surface_pd.thermodynamics.DatasetAlignment
    :members:

``AlignedPhaseDataset`` is a non-mutating view created by
``DatasetAlignment.create_aligned_dataset()``. It retains the root and target
datasets, both explicit anchor identities, the bulk reference, signed bulk-unit
count, and energy offset through its ``alignment`` property. Source
``Phase.dft_energy_ev`` values remain unchanged.

Only direct-to-root alignment is supported initially. Both inputs to
``DatasetAlignment`` must be ordinary ``PhaseDataset`` objects, so aligned
views cannot be chained and cyclic alignment graphs cannot be constructed.

.. autoclass:: surface_pd.thermodynamics.AlignedPhaseDataset
    :members:

Two-dimensional numerical diagrams
===================================

``DiagramAxis`` defines one named state variable independently of its display
label and unit. Coordinates are finite, one-dimensional, strictly increasing,
and contain at least two points. Reversing an axis is a rendering operation and
does not change the thermodynamic state.

.. autoclass:: surface_pd.thermodynamics.DiagramAxis
    :members:

``PhaseDiagramSpecification`` combines distinct x and y axes with finite
scalar fixed conditions. Together, those inputs must match the state variables
required by the ``GrandPotentialModel`` exactly. Missing variables and unused
conditions are rejected before any phase energy is evaluated.

The state mesh uses conventional plotting orientation equivalent to
``numpy.meshgrid(x, y, indexing="xy")``. Every numerical diagram array
therefore has mesh shape

.. math::

    (N_y, N_x),

where x varies along the last dimension and y varies along the first. This
shape is retained for voltage/temperature, voltage/direct-chemical-potential,
temperature/direct-chemical-potential, and future state-variable pairs.

One specification can compare multiple ordinary or aligned datasets in the
declared order. Dataset IDs must be unique. An aligned target is accepted only
when its ordinary root dataset is also supplied, ensuring that every compared
energy uses the declared root convention.

.. autoclass:: surface_pd.thermodynamics.PhaseDiagramSpecification
    :members:

``PhaseDiagramResult`` retains the complete surface-grand-potential tensor
with shape ``(number_of_phases, N_y, N_x)`` and every underlying
``GrandPotentialResult``. Stability is always determined from the intensive
surface grand potential in eV per square angstrom.

At every mesh point, phase :math:`s` is co-stable when

.. math::

    \left|\gamma_s - \min_j\gamma_j\right|
    \le 10^{-10}\ \mathrm{eV/angstrom^2}.

The relative tolerance is zero. This threshold represents floating-point
equality, not uncertainty in the calculated energies. ``stable_phase_mask``
preserves every co-stable phase. ``representative_phase_indices`` chooses the
first stable phase in declared input order solely for deterministic
single-color rendering; it does not replace or alter the tie mask.

The numerical classes do not import Matplotlib. Rendering is a separate layer.

.. autoclass:: surface_pd.thermodynamics.PhaseDiagramResult
    :members:
