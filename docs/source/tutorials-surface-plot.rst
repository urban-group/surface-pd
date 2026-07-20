==================================
Surface Phase Diagram Construction
==================================

Overview
========

The ``surface-pd-plot`` command generates voltage-temperature phase diagrams
for surface structures. These diagrams show which surface composition is
thermodynamically stable under different electrochemical conditions.

Data File Format
================

The input data file should be a whitespace-separated text file with the
following columns:

Required Columns
----------------

* ``structure`` - Structure identifier (string)
* Element columns - Number of each element (e.g., ``Li``, ``Ni``, ``O``)
* ``E`` - DFT total energy (eV)
* ``a``, ``b`` - In-plane lattice parameters (Å)
* ``gamma`` - Angle between a and b axes (degrees)

Example Data File
-----------------

.. code-block:: text

    structure                Li    Ni    O     E           a     b     gamma
    phase-1-symmetrized      16    12    32    -507.505   5.66  5.8   59.16
    phase-2-optimized        14    12    32    -497.526   5.66  5.8   59.16
    phase-3-relaxed          12    12    32    -486.452   5.66  5.8   59.16

.. note::
   The ``E`` column is automatically mapped to ``dft_energy`` internally.
   Both column names are supported for compatibility.

Basic Usage
===========

Command Line Interface
----------------------

.. code-block:: bash

    surface-pd-plot DATA_FILE -L LITHIUM_SPECIES -O OXYGEN_SPECIES

Required Arguments
------------------

* ``DATA_FILE`` - Path to surface phase diagram data file
* ``-L, --lithium-like-species`` - Name of lithium-like species column
* ``-O, --oxygen-like-species`` - Name of oxygen-like species column

Optional Arguments
------------------

* ``-lt, --low-T`` - Low temperature boundary (K)
* ``-ht, --high-T`` - High temperature boundary (K)
* ``-cl, --color-by-Li-content`` - Color the diagram by lithium-like species content
* ``-co, --color-by-O-content`` - Color the diagram by oxygen-like species content
* ``-d, --discharge`` - Treat the data as a discharge phase diagram
* ``-s, --save`` - Save plot to file instead of displaying

Examples
========

Interactive 3D Phase Diagram
-----------------------------

Generate an interactive 3D voltage-temperature phase diagram:

.. code-block:: bash

    surface-pd-plot \
        ./examples/plotting-examples/SCAN-Li-surface.dat \
        -L Li -O O

This will display an interactive matplotlib window showing stable phases
across voltage and temperature space.

Color by Lithium Content
------------------------

Generate a phase diagram colored by lithium content:

.. code-block:: bash

    surface-pd-plot \
        ./examples/plotting-examples/SCAN-Li-surface.dat \
        -L Li -O O --color-by-Li-content

Color by Oxygen Content
-----------------------

Generate a phase diagram colored by oxygen content:

.. code-block:: bash

    surface-pd-plot \
        ./examples/plotting-examples/SCAN-Li-surface.dat \
        -L Li -O O --color-by-O-content

Custom Temperature Range
------------------------

Specify custom temperature boundaries:

.. code-block:: bash

    surface-pd-plot \
        ./examples/plotting-examples/SCAN-Li-surface.dat \
        -L Li -O O \
        --low-T 200 \
        --high-T 400

Save Plot to File
-----------------

Save the generated plot instead of displaying it:

.. code-block:: bash

    surface-pd-plot \
        ./examples/plotting-examples/SCAN-Li-surface.dat \
        -L Li -O O -s

Understanding the Output
========================

The phase diagram shows:

* **X-axis** - Voltage (V vs Li/Li⁺) or composition
* **Y-axis** - Temperature (K) or composition
* **Z-axis** (3D plots) - Surface free energy (eV/Å²)
* **Colors** - Different stable phases at each condition

Stable phases are determined by minimizing the surface free energy, which
depends on:

* DFT total energy
* Chemical potentials (voltage-dependent)
* Temperature-dependent entropic contributions

Reference-Energy Metadata
=========================

The package does not select scientific reference energies from a DFT method
name. Each input file must begin with a complete metadata block before the
tabular header:

.. code-block:: text

    # method = SCAN+rVV10+U; U_Ni=? eV
    # reference_li_ev_per_atom = -2.33333
    # reference_o2_raw_ev_per_molecule = -12.00701
    # reference_o2_correction_ev_per_molecule = 0.0
    # reference_bulk_litmo2_ev_per_formula_unit = -36.8133525

The method is free-text provenance and is never interpreted as a lookup key.
Record relevant settings such as Hubbard U values there; use ``?`` rather than
guessing an unknown value. All four numerical fields are required and must be
finite. The O2 correction must be explicit, including when it is zero.

Troubleshooting
===============

Column Name Errors
------------------

If you encounter ``KeyError: 'dft_energy'``, ensure your data file has either:

* An ``E`` column (will be automatically converted), or
* A ``dft_energy`` column

Both are supported for backward compatibility.

Reference Metadata Errors
-------------------------

Missing, duplicated, unknown, nonnumeric, or non-finite metadata is rejected
before surface energies are calculated. When two files are supplied, their
complete method text and numerical reference values must match exactly.

Missing Required Arguments
---------------------------

Always specify the lithium-like and oxygen-like species:

.. code-block:: bash

    # Required
    -L Li -O O

Advanced Topics
===============

Multiple Data Files
-------------------

To compare charge and discharge behavior, you can process multiple data
files separately and overlay the results:

.. code-block:: bash

    # Generate charge phase diagram
    surface-pd-plot charge-data.dat -L Li -O O

    # Generate discharge phase diagram
    surface-pd-plot discharge-data.dat -L Li -O O --discharge

Custom Reference Energies
-------------------------

Supply references calculated with the same computational settings as the slab
energies. No package source change or method allowlist is involved. Python API
users provide the same values through :class:`surface_pd.plot.ReferenceEnergies`.

See Also
========

* :doc:`surface_energy` - Surface energy calculation details
* :doc:`theory` - Theoretical background
* ``surface-enumeration`` - Generate input structures
* ``generate-discharge-pd`` - Convert charge to discharge data
