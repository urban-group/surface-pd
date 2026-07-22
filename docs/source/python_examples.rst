==========================
Executable Python examples
==========================

These small examples exercise the supported Python API during every
documentation build. They are intentionally deterministic and avoid graphical
windows, output files, and external executables.

Slab and enumeration setup
==========================

Create a toy slab, identify its surface-normal layers, and configure an
enumeration transformation:

.. testcode:: enumeration

    from pymatgen.core import Lattice

    from surface_pd.core import EnumerationSlab, EnumWithComposition

    slab = EnumerationSlab(
        Lattice.tetragonal(3.0, 15.0),
        ["Li", "Li"],
        [[0.0, 0.0, 0.2], [0.0, 0.0, 0.8]],
        enumerated_species=["Li"],
        num_enumerated_layers={"Li": 1},
        symmetric=True,
    )
    populations = slab.layers_finder()["Li"]
    print(sorted(populations.values()))

    enumerator = EnumWithComposition(
        {"Li": {"Li": 0.5}},
        min_cell_size=1,
        max_cell_size=1,
    )
    print((enumerator.min_cell_size, enumerator.max_cell_size))

.. testoutput:: enumeration

    [1, 1]
    (1, 1)

The example configures enumeration but does not call
:meth:`~surface_pd.core.EnumWithComposition.apply_enumeration`, because that
operation requires the external enumlib executables. See
:doc:`tutorials-surface-enum` for the complete command-line workflow.

Grand-potential calculation
===========================

Construct a one-phase dataset and evaluate its surface energy at one voltage
and temperature:

.. testcode:: phase-diagram

    from surface_pd.thermodynamics import (
        ConstantChemicalPotential,
        GrandPotentialModel,
        Phase,
        PhaseDataset,
        ThermodynamicState,
    )
    phase = Phase(
        "example",
        {"A": 2, "B": 1},
        dft_energy_ev=-8.0,
        surface_area_angstrom2=9.0,
        number_of_surfaces=2,
    )
    dataset = PhaseDataset(
        "surface",
        ("A", "B"),
        (phase,),
        "demonstration values (not for scientific use)",
    )
    model = GrandPotentialModel(
        ("A", "B"),
        {
            "A": ConstantChemicalPotential(-2.0),
            "B": ConstantChemicalPotential(-3.0),
        },
        (),
    )
    result = model.evaluate(dataset, ThermodynamicState({}))
    print(result.surface_grand_potential_ev_per_angstrom2.shape)
    print(result.surface_grand_potential_ev_per_angstrom2[0])

.. testoutput:: phase-diagram

    (1,)
    -0.05555555555555555

The numerical values above are deliberately synthetic. Real calculations must
use reference phases and chemical-potential parameters obtained with settings
consistent with the slab energies, as described in :doc:`configuration`.
