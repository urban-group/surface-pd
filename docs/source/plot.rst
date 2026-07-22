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

By default, the renderer uses a continuous atomic-fraction gradient for the
first independent component retained by the numerical result. Independent
component order follows the model's declared component order, making this
choice deterministic. Pass ``coloring="phase_identity"`` to color each
deterministic representative by its qualified ``dataset_id:phase_id``
identity. At a numerical tie, either rendering uses the first phase in declared
input order, while the source result's complete ``stable_phase_mask`` remains
unchanged.

Axis text comes only from each
:class:`~surface_pd.thermodynamics.DiagramAxis` label and unit. Optional axis
inversion changes presentation only and does not mutate coordinates or
thermodynamic results.

Thin black boundaries between differently rendered representative phases are
enabled by default. They follow the cell edges of the same nearest-shaded
representative field used for the colors, so phases with nonconsecutive
internal indices remain unambiguous. ``boundary_color`` accepts a Matplotlib
color, and ``boundary_linewidth`` sets its width in points. Pass
``boundary_color=None`` to draw only the colored mesh. These are presentation
choices: boundary construction neither changes the energy tensor nor discards
the complete tie-preserving stability mask.

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

The named constructors keep common Python calls concise:

.. code-block:: python

    li_fraction = CompositionColoring.atomic_fraction("Li")
    oxygen_fraction = CompositionColoring.atomic_fraction("O")
    li_per_ni = CompositionColoring.component_ratio("Li", "Ni")

Custom labels remain available through the constructors or full dataclass
initialization. Atomic fraction is the general default because it requires no
inferred host component. A component ratio is preferable when the user has an
explicit, scientifically justified invariant reference such as Ni.

Composition coloring follows the same first-in-order representative policy at
ties. The numerical result continues to retain all co-stable phases and the
complete energy tensor.

.. autoclass:: surface_pd.plot.CompositionColoring
    :members:
