======================================
Surface-enumeration command reference
======================================

Overview
========

``surface-enumeration`` enumerates vacancy orderings or substitutional
orderings in selected surface layers. Candidate generation uses enumlib
through pymatgen, filters out transformations that are not valid in-plane
surface cells, and optionally writes finalized structures in VASP POSCAR
format.

Enumlib is an external requirement for this command. Its ``enum.x`` and
``makestr.x`` executables, or a supported ``makeStr`` variant, must be
discoverable on ``PATH``.

The runnable `surface-enumeration walkthrough
<https://github.com/urban-group/surface-pd/blob/main/examples/surface-enumeration-cli.md>`_
lives with its committed inputs under ``examples/``. Commands in that guide
assume ``examples/`` is the current directory.

Slab requirements
=================

The target structure is read by pymatgen and is normally supplied as a VASP
POSCAR. One lattice-vector direction must contain the broken periodicity and
vacuum region; that vector does not need to be perpendicular to the surface
plane. Selective-dynamics flags identify fixed and relaxed regions. Finalized
structures preserve those constraints for sites retained from the parent and
complete flags omitted by enumlib-generated sites.

Symmetric enumeration requires a compatible inversion-symmetric slab with
relaxed outer layers and a fixed interior region. Asymmetric enumeration
modifies only the selected top-surface layers.

Input JSON
==========

The command accepts one JSON object with these required fields:

.. code-block:: json

    {
      "target_slab_path": "path/to/POSCAR",
      "replacements": {
        "Li": {"Li": [1.0, 0.75, 0.5]}
      },
      "num_enumerated_layers": {
        "Li": 1
      },
      "symmetric": false,
      "max_cell_size": 2
    }

``target_slab_path``
    Path to the parent slab. Relative paths are resolved from the command's
    working directory.

``replacements``
    Maps each selected parent species to replacement species and requested
    fractional occupancies. Retaining the parent species below occupancy one
    enumerates vacancies. Naming a different species enumerates substitutions.
    When multiple species or occupancy lists are present, the command evaluates
    their Cartesian product.

``num_enumerated_layers``
    Maps each selected parent species to an independent positive count of its
    outermost species-bearing layers. Counts are species-relative: one Li layer
    and two O layers do not necessarily describe the same geometric planes.

``symmetric``
    Selects whether corresponding layers on both surfaces are finalized as an
    inversion-symmetric model. False selects only the top surface.

``max_cell_size``
    Maximum in-plane area multiplier considered during candidate generation.
    It is not a guaranteed number of finalized structures.

Command usage
=============

.. code-block:: text

    surface-enumeration INPUT_JSON [--generate-poscar]

``--generate-poscar`` or ``-g`` writes accepted structures in
composition-specific directories. Without this option, the command reports
enumeration counts without writing POSCAR files.

Use ``surface-enumeration --help`` for the installed command's current option
summary. Errors in the input configuration, slab compatibility, or external
enumlib execution are reported before finalized structures are returned.

Python API
==========

For direct control over one composition, immutable result metadata, and
programmatic POSCAR export, use :class:`surface_pd.core.EnumerationSlab` and
:class:`surface_pd.core.SurfaceEnumerator`. The maintained `Python API
notebook
<https://github.com/urban-group/surface-pd/blob/main/examples/enumeration-python-api.ipynb>`_
demonstrates that workflow.
