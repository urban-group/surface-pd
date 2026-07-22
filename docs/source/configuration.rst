============================
Versioned JSON configuration
============================

``PhaseDiagramConfiguration`` is the serialization boundary for generalized
phase diagrams. Version 1 keeps thermodynamic model implementations in package
code: JSON selects a reviewed built-in model by name and supplies its validated
parameters. Import paths, Python expressions, and arbitrary callables are not
configuration values.

Every configuration contains exactly these top-level fields:

``schema_version``
    Integer ``1``.
``calculation_method``
    Nonempty, single-line, uninterpreted DFT-method provenance shared by phase
    datasets and reference phases.
``components``
    Ordered unique component names.
``independent_chemical_potentials``
    Mapping from independently modeled components to built-in model objects.
``reference_phases``
    Constant bulk equalities used to solve dependent chemical potentials.
``diagram``
    State-variable identities and default display metadata for two axes.
``datasets``
    Separate whitespace-delimited tables with canonical columns and optional
    name overrides.
``alignments``
    Explicit direct-to-reference dataset alignments.

Unknown and missing fields are rejected recursively. This strictness is part
of the version contract: a misspelled scientific parameter cannot silently
become unused metadata. Matplotlib coloring, colormaps, axis inversion, and
boundary styles are deliberately absent because they do not define the
thermodynamic calculation.

Built-in chemical-potential models
==================================

The model registry contains only:

``constant``
    ``value_ev_per_component`` in eV per component.
``direct``
    ``state_variable`` and ``reference_energy_ev_per_component`` in eV per
    component.
``intercalation_voltage``
    ``reference_energy_ev_per_component``, positive integer
    ``electrons_per_component``, and uninterpreted ``voltage_reference``.
``fixed_pressure_oxygen``
    ``raw_o2_energy_ev_per_molecule`` and ``correction_ev_per_molecule`` in eV
    per O2 molecule.

The ``diagram`` node assigns model state variables to x and y axes and stores
convenient default display metadata:

.. code-block:: json

    {
      "x_axis": {
        "state_variable": "voltage",
        "label": "Potential vs. Li/Li+",
        "unit": "V"
      },
      "y_axis": {
        "state_variable": "temperature",
        "label": "Temperature",
        "unit": "K"
      }
    }

The JSON does not own evaluation ranges, coordinate arrays, mesh density, or
fixed thermodynamic conditions. Python callers provide them for each numerical
evaluation:

.. code-block:: python

    specification = configuration.create_diagram_specification(
        x_values=numpy.linspace(0.0, 5.0, 501),
        y_values=numpy.linspace(1.0, 1500.0, 501),
        fixed_conditions={},
    )

Coordinates may be nonuniform but must be finite, strictly increasing, and
contain at least two values. Required model state variables not assigned to an
axis must be supplied in ``fixed_conditions``. The configured labels and units
initialize the Matplotlib axes and can be replaced through the returned
``Axes`` object.

The complete Draft 2020-12 schema is packaged as
``surface_pd/schemas/phase-diagram-config-v1.schema.json``. Runtime validation
is implemented directly and does not add a JSON-schema library dependency.

JSON ownership and round trips
==============================

``read_json`` retains the absolute source path so table paths can later resolve
relative to the configuration file. ``to_dict`` returns a deep mutable copy;
changing it cannot modify the validated configuration. ``write_json`` emits
the canonical version-1 representation with stable indentation.

Explicit phase-table loading
============================

``load_datasets`` reads each table named by ``datasets`` and returns immutable
datasets in the same order. An absolute table path is used directly. A relative
path is resolved against the directory containing the source JSON file, which
makes a configuration and its tables relocatable as one directory tree.
In-memory configurations therefore require absolute table paths.

Version-1 phase tables are deliberately simple whitespace-delimited text. The
first nonempty line is a unique column header and every later nonempty line is
one phase with the same number of fields. Quoted fields containing whitespace
and comment lines are not supported. In particular, legacy leading-comment
reference metadata is not read: the JSON file is the single source of
thermodynamic configuration and provenance.

The canonical table columns are ``phase_id``, one column named for every
declared component, ``dft_energy_ev``, and ``surface_area_angstrom2``. Extra
columns are ignored. For example:

.. code-block:: text

    phase_id  Li  Ni  O  dft_energy_ev  surface_area_angstrom2  note
    p0        1   1   2  -40.0          12.5                    pristine
    p1        0   1   2  -36.0          12.5                    delithiated

needs only the concise dataset definition

.. code-block:: json

    {
      "dataset_id": "facet_001",
      "path": "tables/facet-001.dat",
      "number_of_surfaces": 2
    }

``number_of_surfaces`` is a positive integer shared by the dataset. It is the
number of equivalent surfaces represented by each calculated cell and is the
explicit divisor in the surface-energy normalization. Two is appropriate for
a symmetric slab with two equivalent surfaces; the field describes the
normalization rather than merely asserting structural symmetry.

External tables with different headings may declare only the names that
differ from the canonical contract:

.. code-block:: json

    {
      "dataset_id": "facet_001",
      "path": "tables/external.dat",
      "number_of_surfaces": 2,
      "column_overrides": {
        "phase_id": "name",
        "composition": {
          "Li": "n_Li",
          "Ni": "n_Ni",
          "O": "n_O"
        },
        "dft_energy_ev": "energy_ev",
        "surface_area_angstrom2": "area_a2"
      }
    }

Counts must be nonnegative integers, total DFT energies must be finite in eV,
and areas must be positive in square angstroms. Resolved source columns must
be distinct. Phase identifiers must be unique within each dataset; their
stable global identities are ``dataset_id:phase_id``.

Configured alignment
====================

An alignment names an ordinary reference dataset, a target dataset, one anchor in
each dataset, and a bulk reference from ``reference_phases``:

.. code-block:: json

    {
      "reference_dataset_id": "facet_001",
      "target_dataset_id": "facet_104",
      "reference_anchor_phase_id": "p0",
      "target_anchor_phase_id": "q0",
      "bulk_reference_id": "bulk-LiNiO2"
    }

Loading replaces the declared target with a non-mutating
:class:`~surface_pd.thermodynamics.AlignedPhaseDataset` view. The raw phase
energies remain available unchanged through ``source_dataset``, while the view
retains the anchors, bulk reference, signed bulk-unit count, and calculated
energy offset. The alignment equation and compatibility conditions are given
in :ref:`dataset-alignment`. Alignments must be direct to an ordinary reference;
chains and cycles are rejected.

.. autoclass:: surface_pd.configuration.PhaseDiagramConfiguration
    :members:
