==================================
Surface phase-diagram command line
==================================

Overview
========

``surface-pd-plot`` evaluates and renders a generalized two-dimensional phase
diagram from one versioned JSON configuration. The configuration names the
thermodynamic models, axes, fixed conditions, phase tables, explicit column
mappings, optional dataset alignments, and rendering choices. No chemical
species has a special command-line role.

See :doc:`configuration` for the complete version-1 format and
:doc:`thermodynamics` for the equations used in evaluation and alignment.

Usage
=====

Saving a figure is explicit:

.. code-block:: bash

    surface-pd-plot CONFIG.json --output diagram.pdf

The filename extension selects any output format supported by Matplotlib, such
as PDF, PNG, or SVG. Interactive display is also explicit:

.. code-block:: bash

    surface-pd-plot CONFIG.json --show

Both actions may be requested together:

.. code-block:: bash

    surface-pd-plot CONFIG.json --output diagram.png --show

At least one of ``--output`` or ``--show`` is required. The command never
chooses an output filename and never displays a window merely because a
configuration was validated.

Workflow
========

The command performs the following reviewed sequence:

1. Read and strictly validate the versioned JSON configuration.
2. Resolve and load every explicitly mapped phase table.
3. Construct any declared direct-to-root aligned dataset views.
4. Evaluate the configured thermodynamic model on the two-dimensional grid.
5. Render phase identity or the configured composition quantity.
6. Save and/or display only after every preceding stage succeeds.

Errors retain their configuration, dataset, table-row, or alignment context
and are reported as command-line errors without a Python traceback. A failed
configuration or calculation creates no requested output and opens no display.

Rendering configuration
=======================

``rendering.coloring.mode`` selects discrete ``phase_identity`` coloring,
``atomic_fraction`` of an explicit component, or an explicit
``component_ratio``. The JSON configuration also owns the colormap and any
x- or y-axis inversion. Axis inversion changes presentation only; it does not
change the thermodynamic state or stable-phase calculation.

Legacy migration
================

The older plotting command accepted one or two Li/O-specific ``.dat`` files,
inferred scientific meanings from column names, reconstructed areas from
lattice parameters, and read reference energies from leading comments. Those
arguments and inference rules are not part of the generalized command.

The old data and implementation remain temporarily for numerical regression in
issue 42. They are not a second supported command-line format and will be
removed before the initial release. Issue 42 will provide the maintained
Li--transition-metal--O configuration and document migration of the committed
legacy datasets after reproducing their validated results.

See also
========

* :doc:`configuration` -- JSON, table, and alignment format
* :doc:`thermodynamics` -- generalized thermodynamic equations
* :doc:`theory` -- surface phase-diagram background
* ``surface-enumeration`` -- enumerate surface structures
