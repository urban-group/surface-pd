====
plot
====

Generalized phase-diagram rendering
***********************************

The generalized renderer consumes an already evaluated
:class:`~surface_pd.thermodynamics.PhaseDiagramResult`. It does not construct a
thermodynamic state, reevaluate energies, select a different stable phase, show
a GUI window, or save a file. It returns the Matplotlib figure, axes, and
colorbar so applications retain control of presentation and output.

With no composition coloring, each deterministic representative phase is
colored by its qualified ``dataset_id:phase_id`` identity. At a numerical tie,
the renderer uses the first phase in declared input order, while the source
result's complete ``stable_phase_mask`` remains unchanged.

Axis text comes only from each
:class:`~surface_pd.thermodynamics.DiagramAxis` label and unit. Optional axis
inversion changes presentation only and does not mutate coordinates or
thermodynamic results.

.. autofunction:: surface_pd.plot.plot_phase_diagram

Explicit composition coloring
*****************************

``CompositionColoring`` supports two explicit, chemistry-independent
normalizations. Atomic-fraction coloring uses

.. math::

    x_i = \frac{n_i}{\sum_j n_j},

where the denominator includes every component in the phase basis.
Component-ratio coloring instead uses

.. math::

    r_{i/k} = \frac{n_i}{n_k},

where both numerator component :math:`i` and denominator component :math:`k`
are named by the user. A zero denominator is rejected. The renderer never
infers a host species, transition metal, occupancy scale, or grouping from row
positions.

Composition coloring follows the same first-in-order representative policy at
ties. The numerical result continues to retain all co-stable phases and the
complete energy tensor.

.. autoclass:: surface_pd.plot.CompositionColoring
    :members:
