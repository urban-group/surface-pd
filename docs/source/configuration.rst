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
    Two axes and finite scalar fixed conditions.
``datasets``
    Separate whitespace-delimited table paths and explicit column mappings.
``alignments``
    Explicit direct-to-root dataset alignments.
``rendering``
    Identity or composition coloring and presentation options.

Unknown and missing fields are rejected recursively. This strictness is part
of the version contract: a misspelled scientific parameter cannot silently
become unused metadata.

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

Axes use either tagged linear coordinates

.. code-block:: json

    {
      "kind": "linear",
      "start": 0.0,
      "stop": 5.0,
      "number": 501
    }

or explicit, strictly increasing values:

.. code-block:: json

    {
      "kind": "values",
      "values": [0.0, 0.25, 1.0, 2.0]
    }

The complete Draft 2020-12 schema is packaged as
``surface_pd/schemas/phase-diagram-config-v1.schema.json``. Runtime validation
is implemented directly and does not add a JSON-schema library dependency.

JSON ownership and round trips
==============================

``read_json`` retains the absolute source path so table paths can later resolve
relative to the configuration file. ``to_dict`` returns a deep mutable copy;
changing it cannot modify the validated configuration. ``write_json`` emits
the canonical version-1 representation with stable indentation.

.. autoclass:: surface_pd.configuration.PhaseDiagramConfiguration
    :members:
