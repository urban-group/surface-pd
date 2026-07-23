"""Tests for the public surface-constrained enumeration workflow."""

from pathlib import Path

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure
from pymatgen.io.vasp import Poscar

from surface_pd.core import EnumerationSlab, SurfaceEnumerator
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
        direction=2,
        enumerated_species=["Li"],
        num_enumerated_layers={"Li": 1},
        symmetric=False,
    )


def test_selected_indices_are_separate_from_fixed_region_bounds():
    """Site selection and fixed-region metadata should have clear APIs."""
    slab = _configured_slab()
    analysis = slab.analyze()

    assert analysis.enumerated_site_indices == {"Li": (1,)}
    assert analysis.fixed_site_indices == (2,)
    assert analysis.fixed_region_bounds_angstrom == (7.5, 7.5)
    assert not hasattr(slab, "index_extraction")
    assert not hasattr(slab, "get_enumerated_site_indices")
    assert not hasattr(slab, "get_fixed_region_bounds")
    assert not hasattr(slab, "get_center_sites")
    assert not hasattr(slab, "layers_finder")
    assert not hasattr(slab, "group_atoms_by_layer")
    assert not hasattr(slab, "layer_distinguisher")
    assert not hasattr(slab, "layers")


def test_enumeration_slab_string_extends_the_structure_summary():
    """Configured slab strings should include surface-enumeration settings."""
    summary = str(_configured_slab())

    assert summary.startswith("Full Formula (Li2 O1)\n")
    assert "Sites (3)" in summary
    assert summary.endswith(
        "Surface enumeration\n"
        "  Vacuum-bearing direction: 2\n"
        "  Layer tolerance: 0.500 Å\n"
        "  Enumerated species: Li\n"
        "  Enumerated layers: Li: 1\n"
        "  Symmetric: no"
    )


def test_unconfigured_enumeration_slab_string_reports_missing_settings():
    """Unset optional enumeration settings should be visible, not invented."""
    slab = EnumerationSlab(
        Lattice.tetragonal(3, 10),
        ["Li"],
        [[0, 0, 0.5]],
    )

    summary = str(slab)

    assert "Enumerated species: not configured" in summary
    assert "Enumerated layers: not configured" in summary
    assert "Symmetric: not configured" in summary


def test_analysis_has_a_readable_string_summary():
    """Analysis should provide a concise report for interactive inspection."""
    summary = str(_configured_slab().analyze())

    assert summary.startswith("SlabAnalysis\n")
    assert "Layers (3):" in summary
    assert "0: 3.000 Å | sites: 0 | species: Li: 1" in summary
    assert "Enumerated sites: Li: 1" in summary
    assert "Fixed sites: 2" in summary
    assert "Fixed region: 7.500–7.500 Å" in summary


def test_analysis_ignores_missing_selective_dynamics_flags():
    """Missing pymatgen site flags should not be mistaken for fixed sites."""
    slab = _configured_slab()
    slab[0].properties["selective_dynamics"] = None

    analysis = slab.analyze()

    assert analysis.fixed_site_indices == (2,)


def test_analysis_snapshot_does_not_change_with_the_source_structure():
    """A snapshot should remain valid after the mutable slab changes."""
    slab = _configured_slab()
    original = slab.analyze()

    slab.translate_sites([0, 1, 2], [0, 0, 0.1], frac_coords=True)
    updated = slab.analyze()

    assert original.fixed_region_bounds_angstrom == (7.5, 7.5)
    assert updated.fixed_region_bounds_angstrom == (9.0, 9.0)


def test_cartesian_analysis_reproduces_symmetric_reference_selection():
    """The Cartesian analysis should retain the reviewed Li2O selection."""
    root = Path(__file__).resolve().parents[1]
    slab = EnumerationSlab.from_file(
        root
        / "examples/enumeration-examples/structure/electrode"
        / "POSCAR_Li2O_110.vasp",
        direction=2,
        layer_tolerance_angstrom=0.5,
        enumerated_species=["Li", "O"],
        num_enumerated_layers={"Li": 1, "O": 2},
        symmetric=True,
    )

    assert slab.analyze().enumerated_site_indices == {
        "Li": (0, 1, 24, 25),
        "O": (26, 27, 37, 38),
    }


def test_cartesian_analysis_reproduces_asymmetric_reference_selection():
    """The Cartesian analysis should retain the reviewed Li selection."""
    root = Path(__file__).resolve().parents[1]
    slab = EnumerationSlab.from_file(
        root
        / "examples/enumeration-examples/structure/electrode"
        / "POSCAR_Li_100.vasp",
        direction=2,
        layer_tolerance_angstrom=0.5,
        enumerated_species=["Li"],
        num_enumerated_layers={"Li": 1},
        symmetric=False,
    )

    analysis = slab.analyze()

    assert analysis.enumerated_site_indices == {"Li": (8,)}
    assert analysis.fixed_site_indices == (0, 1, 2, 3)
    assert analysis.fixed_region_bounds_angstrom == pytest.approx(
        (1.163226969138, 6.295090732845)
    )


def test_analysis_without_fixed_sites_has_no_fixed_bounds():
    """Missing fixed flags should be represented without invented bounds."""
    slab = EnumerationSlab(
        Lattice.tetragonal(3, 10),
        ["Li"],
        [[0, 0, 0.5]],
    )

    analysis = slab.analyze()

    assert analysis.fixed_site_indices == ()
    assert analysis.fixed_region_bounds_angstrom is None


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

    assert [layer.site_indices for layer in short.analyze().layers] == [
        layer.site_indices for layer in long.analyze().layers
    ]
    assert short.analyze().layers[0].coordinate == pytest.approx(2.1)
    assert long.analyze().layers[0].coordinate == pytest.approx(2.1)
    assert isinstance(short.analyze().layers[0].coordinate, float)


def test_layers_use_plane_height_for_a_slanted_vacuum_vector():
    """Layer distances should not use the length of a slanted cell vector."""
    lattice = Lattice([[3, 0, 0], [0, 3, 0], [4, 0, 10]])
    slab = EnumerationSlab(
        lattice,
        ["Li", "Li"],
        [[0, 0, 0.2], [0.5, 0, 0.22]],
        layer_tolerance_angstrom=0.25,
    )

    layers = slab.analyze().layers
    assert len(layers) == 1
    assert layers[0].coordinate == pytest.approx(2.1)
    assert layers[0].species_counts == {"Li": 2}
    with pytest.raises(TypeError):
        layers[0].species_counts["Li"] = 1


def test_layer_detection_handles_the_periodic_cell_boundary():
    """Sites in one layer can occupy opposite sides of the unit cell."""
    slab = EnumerationSlab(
        Lattice.tetragonal(3, 10),
        ["Li", "Li", "O"],
        [[0, 0, 0.98], [0.5, 0.5, 0.01], [0, 0, 0.5]],
        layer_tolerance_angstrom=0.31,
    )

    assert {layer.site_indices for layer in slab.analyze().layers} == {
        (0, 1),
        (2,),
    }


def test_layer_clustering_does_not_chain_beyond_the_tolerance():
    """A series of close neighbors must not create an over-wide layer."""
    slab = EnumerationSlab(
        Lattice.tetragonal(3, 10),
        ["Li", "Li", "Li"],
        [[0, 0, 0.02], [0, 0, 0.04], [0, 0, 0.06]],
        layer_tolerance_angstrom=0.25,
    )

    assert [layer.site_indices for layer in slab.analyze().layers] == [
        (0, 1),
        (2,),
    ]


def test_analysis_rejects_incomplete_enumeration_configuration():
    """Species and layer counts must be configured together."""
    slab = EnumerationSlab(
        Lattice.tetragonal(3, 10),
        ["Li"],
        [[0, 0, 0.5]],
        enumerated_species=["Li"],
    )

    with pytest.raises(ValueError, match="num_enumerated_layers"):
        slab.analyze()


def test_analysis_rejects_missing_species_and_excessive_layer_counts():
    """Layer requests must be realizable on the configured slab."""
    missing = EnumerationSlab(
        Lattice.tetragonal(3, 10),
        ["Li"],
        [[0, 0, 0.5]],
        enumerated_species=["O"],
        num_enumerated_layers={"O": 1},
        symmetric=False,
    )
    excessive = EnumerationSlab(
        Lattice.tetragonal(3, 10),
        ["Li"],
        [[0, 0, 0.5]],
        enumerated_species=["Li"],
        num_enumerated_layers={"Li": 2},
        symmetric=False,
    )

    with pytest.raises(ValueError, match="O.*no detected layers"):
        missing.analyze()
    with pytest.raises(ValueError, match="requests 2.*only 1"):
        excessive.analyze()


def test_symmetric_analysis_rejects_overlapping_surface_layers():
    """Top and bottom selections must not overlap for symmetric slabs."""
    slab = EnumerationSlab(
        Lattice.tetragonal(3, 10),
        ["Li", "Li", "Li"],
        [[0, 0, 0.2], [0, 0, 0.5], [0, 0, 0.8]],
        enumerated_species=["Li"],
        num_enumerated_layers={"Li": 2},
        symmetric=True,
    )

    with pytest.raises(ValueError, match="non-overlapping layers per surface"):
        slab.analyze()


def test_public_symmetry_check_has_a_boolean_contract():
    """Public inversion-symmetry inspection should always return a Boolean."""
    slab = _configured_slab()

    assert isinstance(slab.has_inversion_symmetry(), bool)
    assert not hasattr(slab, "is_symmetry")


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

    requested_raw_limit = None

    def fake_raw_enumeration(*args, **kwargs):
        nonlocal requested_raw_limit
        requested_raw_limit = args[-1]
        return [
            {"structure": doubled_normal},
            {"structure": tilted},
            {"structure": valid},
        ]

    analyze_calls = 0
    original_analyze = EnumerationSlab.analyze

    def counted_analyze(structure):
        nonlocal analyze_calls
        analyze_calls += 1
        return original_analyze(structure)

    monkeypatch.setattr(
        "surface_pd.core.enum._apply_raw_enumeration",
        fake_raw_enumeration,
    )
    monkeypatch.setattr(EnumerationSlab, "analyze", counted_analyze)
    enumerator = SurfaceEnumerator(
        {"Li": {"Li": 0.5}}, min_cell_size=2, max_cell_size=2
    )

    results = enumerator.apply_enumeration(slab, max_structures=7)

    assert len(results) == 1
    assert requested_raw_limit == 7
    assert analyze_calls == 1
    assert np.allclose(results[0].lattice.matrix[2], slab.lattice.matrix[2])
    assert all("X" not in str(species) for species in results[0].species)
    assert results[0].enumeration_metadata.transformation_matrix == (
        (2, 0, 0),
        (0, 1, 0),
        (0, 0, 1),
    )
    assert results[0].enumeration_metadata.area_multiplier == 2
    assert results[0].enumeration_metadata.raw_candidate_rank == 2
    assert results[0].enumeration_metadata.symmetric is False
    with pytest.raises(AttributeError):
        results[0].enumeration_metadata.area_multiplier = 3


def test_finalized_candidates_retain_raw_order_without_deduplication(
    monkeypatch,
):
    """Accepted raw candidates should remain separate and retain rank."""
    slab = _configured_slab()
    candidate = slab.copy()
    candidate.make_supercell([[2, 0, 0], [0, 1, 0], [0, 0, 1]])

    monkeypatch.setattr(
        "surface_pd.core.enum._apply_raw_enumeration",
        lambda *args, **kwargs: [
            {"structure": candidate.copy()},
            {"structure": candidate.copy()},
        ],
    )
    results = SurfaceEnumerator(
        {"Li": {"Li": 0.5}}, min_cell_size=2, max_cell_size=2
    ).apply_enumeration(slab)

    assert len(results) == 2
    assert [
        result.enumeration_metadata.raw_candidate_rank for result in results
    ] == [0, 1]


def test_asymmetric_finalization_restores_missing_selective_dynamics(
    monkeypatch, tmp_path
):
    """Final results should remain analyzable and writable as POSCAR files."""
    slab = _configured_slab()
    candidate = slab.copy()
    candidate.make_supercell([[2, 0, 0], [0, 1, 0], [0, 0, 1]])
    for site in candidate:
        site.properties["selective_dynamics"] = None
    monkeypatch.setattr(
        "surface_pd.core.enum._apply_raw_enumeration",
        lambda *args, **kwargs: [{"structure": candidate}],
    )

    result = SurfaceEnumerator(
        {"Li": {"Li": 0.5}}, min_cell_size=2, max_cell_size=2
    ).apply_enumeration(slab)[0]

    assert all(
        site.properties["selective_dynamics"] is not None for site in result
    )
    assert sum(
        not any(site.properties["selective_dynamics"]) for site in result
    ) == 2
    assert len(result.analyze().fixed_site_indices) == 2
    Poscar(result).write_file(tmp_path / "POSCAR")


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
    enumerator = SurfaceEnumerator(
        {"Li": {"Li": 0.5}}, min_cell_size=2, max_cell_size=2
    )

    with pytest.raises(IncompatibleSymmError):
        enumerator.apply_enumeration(slab)

    assert raw_called is False


def test_symmetric_enumeration_rejects_incomplete_selective_dynamics(
    monkeypatch,
):
    """Every symmetric-slab site needs complete Boolean relaxation flags."""
    slab = _configured_slab()
    slab.symmetric = True
    del slab[1].properties["selective_dynamics"]
    raw_called = False

    def fail_if_called(*args, **kwargs):
        nonlocal raw_called
        raw_called = True
        return []

    monkeypatch.setattr(
        "surface_pd.core.enum._apply_raw_enumeration", fail_if_called
    )

    with pytest.raises(ValueError, match="selective_dynamics"):
        SurfaceEnumerator(
            {"Li": {"Li": 0.5}}, min_cell_size=2, max_cell_size=2
        ).apply_enumeration(slab)

    assert raw_called is False


@pytest.mark.parametrize(
    "kwargs",
    [
        {"replacements": {}},
        {"replacements": {"Li": {"Li": -0.1}}},
        {"replacements": {"Li": {"Li": 1.1}}},
        {"replacements": {"Li": {"Li": float("nan")}}},
        {"replacements": {"Li": {"Li": 0.6, "Na": 0.5}}},
        {"replacements": {"Li": {"Li": 1}}, "min_cell_size": 0},
        {
            "replacements": {"Li": {"Li": 1}},
            "min_cell_size": 2,
            "max_cell_size": 1,
        },
        {
            "replacements": {"Li": {"Li": 1}},
            "enum_precision_parameter": 0,
        },
    ],
)
def test_surface_enumerator_rejects_invalid_configuration(kwargs):
    """Malformed occupancy and cell-size requests should fail immediately."""
    with pytest.raises((TypeError, ValueError)):
        SurfaceEnumerator(**kwargs)


def test_surface_enumerator_owns_validated_replacements():
    """Caller mutation must not alter validated enumeration settings."""
    replacements = {"Li": {"Li": 0.5}}
    enumerator = SurfaceEnumerator(replacements)

    replacements["Li"]["Li"] = 1.0

    assert enumerator.replacements == {"Li": {"Li": 0.5}}
    with pytest.raises(TypeError):
        enumerator.replacements["Li"]["Li"] = 1.0


def test_surface_enumerator_rejects_unrealizable_occupancy(monkeypatch):
    """No allowed area multiplier can realize one-third of one site."""
    slab = _configured_slab()
    enumerator = SurfaceEnumerator(
        {"Li": {"Li": 1 / 3}}, min_cell_size=1, max_cell_size=2
    )
    raw_called = False

    def fail_if_called(*args, **kwargs):
        nonlocal raw_called
        raw_called = True
        return []

    monkeypatch.setattr(
        "surface_pd.core.enum._apply_raw_enumeration", fail_if_called
    )

    with pytest.raises(ValueError, match="cannot be realized"):
        enumerator.apply_enumeration(slab)

    assert raw_called is False


@pytest.mark.parametrize("max_structures", [True, 0, -1, 1.5])
def test_surface_enumerator_validates_raw_candidate_limit(max_structures):
    """The raw enumlib candidate limit must be a positive integer."""
    enumerator = SurfaceEnumerator({"Li": {"Li": 1.0}})

    with pytest.raises((TypeError, ValueError)):
        enumerator.apply_enumeration(
            _configured_slab(), max_structures=max_structures
        )
