==================================
Surface phase-diagram command line
==================================

Overview
========

``surface-pd-plot`` evaluates and renders a generalized two-dimensional phase
diagram from one versioned JSON configuration. The configuration names the
thermodynamic models, axis state variables, phase tables, explicit column
mappings, and optional dataset alignments. No chemical
species has a special command-line role.

See :doc:`configuration` for the complete version-1 format and
:doc:`thermodynamics` for the equations used in evaluation and alignment.

Usage
=====

Saving a figure is explicit:

.. code-block:: bash

    surface-pd-plot CONFIG.json \
        --x-range 0 5 --y-range 1 1500 \
        --output diagram.pdf

The filename extension selects any output format supported by Matplotlib, such
as PDF, PNG, or SVG. Interactive display is also explicit:

.. code-block:: bash

    surface-pd-plot CONFIG.json \
        --x-range 0 5 --y-range 1 1500 \
        --show

Both actions may be requested together:

.. code-block:: bash

    surface-pd-plot CONFIG.json \
        --x-range 0 5 --y-range 1 1500 \
        --output diagram.png --show

At least one of ``--output`` or ``--show`` is required. The command never
chooses an output filename and never displays a window merely because a
configuration was validated.

Both evaluation ranges are required because the generalized command cannot
infer scientifically meaningful voltage, temperature, or chemical-potential
domains. ``--mesh-points N`` selects one shared linear density for quick
inspection and defaults to 201 points on each axis. Repeat
``--condition STATE_VARIABLE=VALUE`` for any model variables not assigned to
an axis. Unequal or nonuniform coordinate arrays belong in the Python API.

Workflow
========

The command performs the following reviewed sequence:

1. Read and strictly validate the versioned JSON configuration.
2. Resolve and load every explicitly mapped phase table.
3. Construct any declared direct-to-reference aligned dataset views.
4. Evaluate the configured thermodynamic model on the two-dimensional grid.
5. Render the default composition gradient or a command-line composition
   choice.
6. Save and/or display only after every preceding stage succeeds.

Errors retain their configuration, dataset, table-row, or alignment context
and are reported as command-line errors without a Python traceback. A failed
configuration or calculation creates no requested output and opens no display.

Quick-inspection coloring
=========================

The command defaults to an atomic-fraction gradient for the first independent
component in declared component order. ``--color-component COMPONENT`` selects
a different independent component, while ``--phase-identity-colors`` selects
discrete stable-phase colors. The two options are mutually exclusive.

The command intentionally does not expose detailed Matplotlib styling. Use
the Python :func:`surface_pd.plot.plot_phase_diagram` interface for component
ratios, colormaps, axis inversion, boundary styling, figure dimensions, and
publication-quality output. None of these choices belongs to the scientific
JSON configuration or changes the stable-phase calculation.

Battery examples and legacy migration
=====================================

The older plotting command accepted one or two Li/O-specific ``.dat`` files,
inferred scientific meanings from column names, reconstructed areas from
lattice parameters, and read reference energies from leading comments. Those
arguments and inference rules are not part of the generalized command.

The maintained PBE+U and SCAN+rVV10+U inputs have been migrated under
``examples/plotting-examples``. Each directory contains version-1 JSON and
tables with explicit surface areas and numbers of represented surfaces. From
the repository's ``examples/`` directory, for example:

.. code-block:: bash

    surface-pd-plot \
        plotting-examples/lno-001-scan/charge.json \
        --x-range 0 5 --y-range 1 1500 \
        --output lno-001-charge.pdf

The chemistry-specific implementation and metadata reader have been removed.
Frozen numerical values preserve regression coverage without retaining a
second implementation or supported input format.

Irreversible oxygen loss during charge is a restriction on which phases remain
accessible during subsequent discharge. It does not change any phase's DFT
energy. The migrated discharge tables therefore contain a filtered subset of
the charge candidates with their original energies.

The accessibility rule is system-specific user policy. One transparent choice
is to normalize oxygen inventory to an explicitly selected invariant host
component. If Ni is invariant, a user might filter records before writing the
discharge table as follows:

.. code-block:: python

    retained_oxygen_per_host = 1.75
    accessible = [
        phase
        for phase in charge_phases
        if phase["O"] / phase["Ni"] <= retained_oxygen_per_host
    ]

The threshold must come from the scientific model of the degradation history;
the package does not choose it. A zero DFT energy is an ordinary numeric input,
not an exclusion marker. The discharge configuration uses the ordinary
evaluator, and optional voltage-axis inversion changes presentation only.

See also
========

* :doc:`configuration` -- JSON, table, and alignment format
* :doc:`thermodynamics` -- generalized thermodynamic equations
* :doc:`theory` -- surface phase-diagram background
* ``surface-enumeration`` -- enumerate surface structures
