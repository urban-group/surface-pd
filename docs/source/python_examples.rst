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
        to_be_enumerated_species=["Li"],
        num_enumerated_layers={"Li": 1},
        symmetric=True,
    )
    populations = slab.layers_finder()["Li"]
    print(sorted(populations.values()))

    enumerator = EnumWithComposition(
        {"X": {"Li": 0.5, "O": 0.5}},
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

Surface-energy calculation
==========================

Construct a one-phase dataset and evaluate its surface energy at one voltage
and temperature:

.. testcode:: phase-diagram

    import numpy as np
    import pandas as pd

    from surface_pd.plot import PdData, ReferenceEnergies

    references = ReferenceEnergies(
        method="demonstration values (not for scientific use)",
        li_ev_per_atom=-2.0,
        o2_raw_ev_per_molecule=-10.0,
        o2_correction_ev_per_molecule=0.0,
        bulk_litmo2_ev_per_formula_unit=-20.0,
    )
    dataframe = pd.DataFrame(
        {
            "Li": [4],
            "Ni": [4],
            "O": [8],
            "E": [-80.0],
            "a": [3.0],
            "b": [3.0],
            "gamma": [90.0],
        }
    )
    phases = PdData(dataframe, "Li", "O", references)
    phases.standardize_pd_data()
    energies = phases.get_surface_energy(
        V=np.array([0.0]),
        T=np.array([298.15]),
    )
    print(energies.shape)
    print(bool(np.isfinite(energies).all()))

.. testoutput:: phase-diagram

    (1, 1)
    True

The numerical references above are deliberately synthetic. Real calculations
must use reference energies obtained with computational settings consistent
with the slab energies, as described in :doc:`tutorials-surface-plot`.
