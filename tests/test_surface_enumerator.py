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
