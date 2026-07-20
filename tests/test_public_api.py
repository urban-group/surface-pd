"""Contract tests for the intentionally supported Python API."""

from pymatgen.core.lattice import Lattice

import surface_pd
from surface_pd import analysis, core, error, plot, util
from surface_pd.core.post_check import PostCheck
from surface_pd.core.slab import Slab
from surface_pd.plot import pd_data, surface_energy
from surface_pd.plot.reference_energies import ReferenceEnergies


def test_package_root_remains_lightweight():
    """The package root should expose metadata, not eager imports."""
    assert surface_pd.__all__ == ["__version__"]


def test_subpackage_exports_define_the_supported_api():
    """Only domain-facing classes and raised errors are public exports."""
    assert core.__all__ == ["Slab", "EnumWithComposition"]
    assert plot.__all__ == ["PdData", "ReferenceEnergies", "SurfaceEnergy"]
    assert analysis.__all__ == []
    assert util.__all__ == []
    assert error.__all__ == [
        "NoInversionSymmetryError",
        "SlabOrientationError",
        "NonCentralInversionSymmetryError",
        "PrimitiveStructureFinderError",
        "NonDefinedSelectiveDynamicsError",
        "NonSlabError",
        "TooLargeSlabError",
        "InvalidCompositionError",
        "IncompatibleSymmError",
        "NonIntegerError",
        "InvalidInputFormatError",
    ]


def test_legacy_internal_exports_are_not_public():
    """Workflow and parsing helpers are not part of the compatibility API."""
    assert not hasattr(core, "PreCheck")
    assert not hasattr(core, "PostCheck")
    assert not hasattr(plot, "find_stable_phases")
    assert not hasattr(util, "csv2dict")
    assert not hasattr(analysis, "structure_filter")


def test_slab_uses_corrected_enumerated_layer_name():
    """The public Slab property should use the complete spelling."""
    slab = Slab(
        Lattice.cubic(3),
        ["Li"],
        [[0, 0, 0]],
        _num_enumerated_layers={"Li": 1},
    )

    assert slab.num_enumerated_layers == {"Li": 1}
    assert not hasattr(slab, "num_layers_enumed")


def test_slab_uses_corrected_supplemental_structure_method_name():
    """The public Slab method should use an action-oriented name."""
    assert hasattr(Slab, "generate_supplemental_structures")
    assert not hasattr(Slab, "supplemental_structures_gene")


def test_removed_post_check_compatibility_wrapper_is_not_available():
    """The internal post-check pipeline has one supported repair method."""
    assert hasattr(PostCheck, "repair_refined_slab_geometry")
    assert not hasattr(PostCheck, "slab_size_check")


def test_energy_references_are_private_and_shared():
    """Reference energies should be public values, not built-in tables."""
    assert not hasattr(surface_energy, "E_O2_by_funtional")
    assert not hasattr(surface_energy, "E_bulk_by_funtional")
    assert not hasattr(surface_energy, "E_Li_by_funtional")
    assert plot.ReferenceEnergies is ReferenceEnergies
    assert pd_data.ReferenceEnergies is ReferenceEnergies
    assert surface_energy.ReferenceEnergies is ReferenceEnergies
