"""Contract tests for the intentionally supported Python API."""

import inspect

from pymatgen.core.lattice import Lattice

import surface_pd
from surface_pd import analysis, core, error, plot, util
from surface_pd.core.post_check import PostCheck
from surface_pd.core.slab import Slab
from surface_pd.plot import pd_data, surface_energy
from surface_pd.plot.reference_energies import ReferenceEnergies

PUBLIC_CLASSES = (
    core.Slab,
    core.EnumWithComposition,
    plot.PdData,
    plot.ReferenceEnergies,
    plot.SurfaceEnergy,
)


def _owned_public_members(cls):
    """Yield public methods and properties defined directly on *cls*."""
    for name, member in cls.__dict__.items():
        if name.startswith("_"):
            continue
        if isinstance(member, property):
            yield name, member.fget
        elif inspect.isfunction(member):
            yield name, member


def _section_entries(docstring, section):
    """Return top-level entry names from one NumPy-style section."""
    lines = inspect.cleandoc(docstring).splitlines()
    try:
        start = lines.index(section)
    except ValueError:
        return set()
    if start + 1 >= len(lines) or set(lines[start + 1]) != {"-"}:
        return set()

    entries = set()
    for line in lines[start + 2 :]:
        if line and not line.startswith(" "):
            if entries and ":" not in line:
                break
            names = line.split(":", 1)[0].split(",")
            entries.update(name.strip() for name in names)
    return entries


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


def test_supported_api_has_consistent_structural_docstrings():
    """Every supported object should expose a useful NumPy-style docstring."""
    for cls in PUBLIC_CLASSES:
        class_doc = inspect.getdoc(cls)
        assert class_doc, cls.__qualname__
        assert "Args:" not in class_doc, cls.__qualname__

        constructor_parameters = {
            name
            for name in inspect.signature(cls).parameters
            if name not in {"args", "kwargs"}
        }
        assert constructor_parameters <= _section_entries(
            class_doc, "Parameters"
        ), cls.__qualname__

        for name, member in _owned_public_members(cls):
            member_doc = inspect.getdoc(member)
            qualified_name = f"{cls.__qualname__}.{name}"
            assert member_doc, qualified_name
            assert "Args:" not in member_doc, qualified_name
            parameters = {
                parameter_name
                for parameter_name in inspect.signature(member).parameters
                if parameter_name != "self"
            }
            assert parameters <= _section_entries(
                member_doc, "Parameters"
            ), qualified_name

    for name in error.__all__:
        exception_doc = inspect.getdoc(getattr(error, name))
        assert exception_doc, name
        assert exception_doc.startswith("Raised when"), name
