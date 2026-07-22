"""Tests for the public surface-constrained enumeration workflow."""

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure

from surface_pd.core import EnumerationSlab, EnumWithComposition
from surface_pd.error import IncompatibleSymmError


def _configured_slab():
    structure = Structure(
        Lattice.tetragonal(3.0, 15.0),
        ["Li", "Li", "O"],
        [[0, 0, 0.2], [0, 0, 0.8], [0.5, 0.5, 0.5]],
        site_properties={
            "selective_dynamics": [
                [True, True, True],
                [True, True, True],
                [False, False, False],
            ]
        },
    )
    return EnumerationSlab.from_structure(
        structure,
        enumerated_species=["Li"],
        num_enumerated_layers={"Li": 1},
        symmetric=False,
    )


def test_selected_indices_are_separate_from_fixed_region_bounds():
    """Site selection and fixed-region metadata should have clear APIs."""
    slab = _configured_slab()

    assert slab.get_enumerated_site_indices() == {"Li": [1]}
    assert slab.get_fixed_region_bounds() == (0.5, 0.5)
    assert not hasattr(slab, "index_extraction")
    assert not hasattr(slab, "layers_finder")
    assert not hasattr(slab, "group_atoms_by_layer")
    assert not hasattr(slab, "layer_distinguisher")


def test_cartesian_layers_are_independent_of_vacuum_length():
    """Changing only the vacuum should not change physical layer grouping."""
    short = EnumerationSlab(
        Lattice.tetragonal(3, 10),
        ["Li", "Li", "O"],
        [[0, 0, 0.2], [0.5, 0.5, 0.22], [0, 0, 0.6]],
        layer_tolerance_angstrom=0.25,
    )
    long = EnumerationSlab(
        Lattice.tetragonal(3, 20),
        ["Li", "Li", "O"],
        [[0, 0, 0.1], [0.5, 0.5, 0.11], [0, 0, 0.3]],
        layer_tolerance_angstrom=0.25,
    )

    assert [layer.site_indices for layer in short.layers] == [
        layer.site_indices for layer in long.layers
    ]
    assert short.layers[0].coordinate == pytest.approx(2.1)
    assert long.layers[0].coordinate == pytest.approx(2.1)
    assert isinstance(short.layers[0].coordinate, float)


def test_layers_use_plane_height_for_a_slanted_vacuum_vector():
    """Layer distances should not use the length of a slanted cell vector."""
    lattice = Lattice([[3, 0, 0], [0, 3, 0], [4, 0, 10]])
    slab = EnumerationSlab(
        lattice,
        ["Li", "Li"],
        [[0, 0, 0.2], [0.5, 0, 0.22]],
        layer_tolerance_angstrom=0.25,
    )

    assert len(slab.layers) == 1
    assert slab.layers[0].coordinate == pytest.approx(2.1)
    assert slab.layers[0].species_counts == {"Li": 2}
    with pytest.raises(TypeError):
        slab.layers[0].species_counts["Li"] = 1


def test_layer_detection_handles_the_periodic_cell_boundary():
    """Sites in one layer can occupy opposite sides of the unit cell."""
    slab = EnumerationSlab(
        Lattice.tetragonal(3, 10),
        ["Li", "Li", "O"],
        [[0, 0, 0.98], [0.5, 0.5, 0.01], [0, 0, 0.5]],
        layer_tolerance_angstrom=0.31,
    )

    assert {layer.site_indices for layer in slab.layers} == {(0, 1), (2,)}


def test_layer_clustering_does_not_chain_beyond_the_tolerance():
    """A series of close neighbors must not create an over-wide layer."""
    slab = EnumerationSlab(
        Lattice.tetragonal(3, 10),
        ["Li", "Li", "Li"],
        [[0, 0, 0.02], [0, 0, 0.04], [0, 0, 0.06]],
        layer_tolerance_angstrom=0.25,
    )

    assert [layer.site_indices for layer in slab.layers] == [(0, 1), (2,)]


def test_surface_enumerator_rejects_normal_multiplication_and_mixing(
    monkeypatch,
):
    """Only strict in-plane derivatives should reach public results."""
    slab = _configured_slab()
    valid = slab.copy()
    valid.make_supercell([[2, 0, 0], [0, 1, 0], [0, 0, 1]])
    doubled_normal = slab.copy()
    doubled_normal.make_supercell([[1, 0, 0], [0, 1, 0], [0, 0, 2]])
    tilted = slab.copy()
    tilted.make_supercell([[1, 0, 1], [0, 2, 0], [0, 0, 1]])

    def fake_raw_enumeration(*args, **kwargs):
        return [
            {"structure": doubled_normal},
            {"structure": tilted},
            {"structure": valid},
        ]

    monkeypatch.setattr(
        "surface_pd.core.enum._apply_raw_enumeration",
        fake_raw_enumeration,
    )
    enumerator = EnumWithComposition(
        {"Li": {"Li": 0.5}}, min_cell_size=2, max_cell_size=2
    )

    results = enumerator.apply_enumeration(slab)

    assert len(results) == 1
    assert np.allclose(results[0].lattice.matrix[2], slab.lattice.matrix[2])
    assert all("X" not in str(species) for species in results[0].species)


def test_symmetric_enumeration_rejects_one_sided_relaxation(monkeypatch):
    """Symmetric enumeration should validate both relaxed surfaces first."""
    slab = _configured_slab()
    slab[0].properties["selective_dynamics"] = [False, False, False]
    slab.symmetric = True
    raw_called = False

    def fail_if_called(*args, **kwargs):
        nonlocal raw_called
        raw_called = True
        return []

    monkeypatch.setattr(
        "surface_pd.core.enum._apply_raw_enumeration", fail_if_called
    )
    enumerator = EnumWithComposition(
        {"Li": {"Li": 0.5}}, min_cell_size=2, max_cell_size=2
    )

    with pytest.raises(IncompatibleSymmError):
        enumerator.apply_enumeration(slab)

    assert raw_called is False
