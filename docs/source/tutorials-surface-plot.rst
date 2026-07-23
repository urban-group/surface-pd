=======================================
Surface phase-diagram command reference
=======================================

Overview
========

``surface-pd-plot`` evaluates and renders a generalized two-dimensional phase
diagram from one versioned JSON configuration. The configuration declares the
thermodynamic model, axis state variables, phase tables, explicit column
mappings, reference energies, and optional dataset alignments. No chemical
species has a special command-line role.

See :doc:`configuration` for the complete version-1 configuration format and
:doc:`thermodynamics` for the evaluation and alignment equations.

The runnable `phase-diagram command-line walkthrough
<https://github.com/urban-group/surface-pd/blob/main/examples/phase-diagram-cli.md>`_
lives with the maintained configurations and tables under ``examples/``.
Commands in that guide assume ``examples/`` is the current directory.

Usage
=====

.. code-block:: text

    surface-pd-plot CONFIG \
        --x-range MIN MAX \
        --y-range MIN MAX \
        [--output FIGURE] [--show]

Both ranges are required because the generalized command cannot infer
scientifically meaningful voltage, temperature, or chemical-potential
domains. At least one of ``--output`` or ``--show`` is required. The command
does not invent an output filename or open a window merely because the
configuration is valid.

``--output``
    Saves the figure. Its extension selects any format supported by
    Matplotlib, such as PDF, PNG, or SVG.

``--show``
    Displays the figure interactively. It may be combined with ``--output``.

``--mesh-points N``
    Selects one shared linear density for quick inspection. The default is 201
    points on each axis. Unequal or nonuniform coordinate arrays belong in the
    Python API.

``--condition STATE_VARIABLE=VALUE``
    Fixes a configured state variable that is not assigned to an axis. Repeat
    the option for multiple fixed conditions.

Use ``surface-pd-plot --help`` for the installed command's complete current
option summary.

Workflow guarantees
===================

The command performs this sequence:

1. Read and strictly validate the versioned JSON configuration.
2. Resolve and load every explicitly mapped phase table.
3. Construct declared direct-to-reference aligned dataset views.
4. Evaluate the configured thermodynamic model on the requested grid.
5. Render the selected composition or phase-identity coloring.
6. Save or display only after every preceding stage succeeds.

Errors retain configuration, dataset, table-row, or alignment context and are
reported as command-line errors without a Python traceback. A failed
configuration or calculation creates no requested output and opens no display.

Coloring
========

The default is an atomic-fraction gradient for the first independent component
in declared component order. ``--color-component COMPONENT`` selects another
independent component. ``--phase-identity-colors`` instead assigns discrete
colors to stable phases. These options are mutually exclusive.

Detailed Matplotlib styling intentionally remains in the Python
:func:`surface_pd.plot.plot_phase_diagram` interface. Colormaps, component
ratios, axis inversion, boundary styling, and figure dimensions are
presentation choices; they do not belong to the scientific JSON configuration
or alter stable-phase evaluation.

Dataset accessibility
=====================

Irreversible oxygen loss or another degradation history can restrict which
phases remain accessible in a later calculation. Such accessibility is
system-specific user policy, not a modification of DFT energies. A discharge
table can therefore contain a filtered subset of charge candidates while
retaining their original energies.

Surface-pd does not infer this filtering rule from chemical species or column
names. A zero DFT energy is an ordinary numeric value, not an exclusion marker.
The versioned configuration and explicit phase table determine exactly which
candidates enter evaluation.

Python API
==========

Use the Python API for unequal grids, reusable thermodynamic results, detailed
rendering control, or application-specific preprocessing. The maintained
`phase-diagram Python API notebook
<https://github.com/urban-group/surface-pd/blob/main/examples/phase-diagram-python-api.ipynb>`_
evaluates and renders a committed dataset step by step.
