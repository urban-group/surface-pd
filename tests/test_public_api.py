"""Contract tests for the intentionally supported Python API."""

import ast
import inspect
from pathlib import Path

import pytest
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.surface import Slab as PymatgenSlab

import surface_pd
from surface_pd import (
    analysis,
    configuration,
    core,
    error,
    plot,
    thermodynamics,
    util,
)
from surface_pd.core.enumeration_slab import EnumerationSlab
from surface_pd.core.post_check import PostCheck

PUBLIC_CLASSES = (
    core.EnumerationSlab,
    core.SurfaceEnumerationMetadata,
    core.SurfaceEnumerator,
    plot.CompositionColoring,
    configuration.PhaseDiagramConfiguration,
    thermodynamics.ThermodynamicState,
    thermodynamics.ConstantChemicalPotential,
    thermodynamics.DirectChemicalPotential,
    thermodynamics.IntercalationChemicalPotential,
    thermodynamics.FixedPressureOxygenChemicalPotential,
    thermodynamics.Phase,
    thermodynamics.PhaseDataset,
    thermodynamics.ReferencePhase,
    thermodynamics.GrandPotentialModel,
    thermodynamics.GrandPotentialResult,
    thermodynamics.DatasetAlignment,
    thermodynamics.AlignedPhaseDataset,
    thermodynamics.DiagramAxis,
    thermodynamics.PhaseDiagramSpecification,
    thermodynamics.PhaseDiagramResult,
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
    assert core.__all__ == [
        "EnumerationSlab",
        "SlabAnalysis",
        "SlabLayer",
        "SurfaceEnumerationMetadata",
        "SurfaceEnumerator",
    ]
    assert plot.__all__ == [
        "CompositionColoring",
        "plot_phase_diagram",
    ]
    assert configuration.__all__ == ["PhaseDiagramConfiguration"]
    assert thermodynamics.__all__ == [
        "AlignedPhaseDataset",
        "ChemicalPotentialModel",
        "ConstantChemicalPotential",
        "DirectChemicalPotential",
        "FixedPressureOxygenChemicalPotential",
        "IntercalationChemicalPotential",
        "DatasetAlignment",
        "DiagramAxis",
        "GrandPotentialModel",
        "GrandPotentialResult",
        "Phase",
        "PhaseDataset",
        "PhaseDiagramResult",
        "PhaseDiagramSpecification",
        "ReferencePhase",
        "ThermodynamicState",
    ]
    assert analysis.__all__ == []
    assert not hasattr(core, "EnumWithComposition")
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
    """The public slab property should use the complete spelling."""
    slab = EnumerationSlab(
        Lattice.cubic(3),
        ["Li"],
        [[0, 0, 0]],
        num_enumerated_layers={"Li": 1},
    )

    assert slab.num_enumerated_layers == {"Li": 1}
    assert not hasattr(slab, "num_layers_enumed")


def test_enumeration_slab_pipeline_helpers_are_private():
    """Legacy CLI machinery should not appear in the supported slab API."""
    helper_names = {
        "wrap_pbc",
        "generate_supplemental_structures",
        "symmetrize_top_base",
        "get_max_min_c_frac",
        "check_rotate",
        "calculate_num_sites",
        "tune_isc",
        "tune_c",
        "add_selective_dynamics",
    }

    for name in helper_names:
        assert not hasattr(EnumerationSlab, name)
        assert hasattr(EnumerationSlab, f"_{name}")
    assert not hasattr(EnumerationSlab, "supplemental_structures_gene")


def test_enumeration_slab_is_distinct_from_pymatgen_slab():
    """The public slab type should not imply pymatgen slab metadata."""
    assert hasattr(core, "EnumerationSlab")
    assert not hasattr(core, "Slab")
    assert issubclass(core.EnumerationSlab, Structure)
    assert not issubclass(core.EnumerationSlab, PymatgenSlab)


def test_package_constructors_do_not_expose_private_parameter_names():
    """Explicit constructors should not leak underscore-prefixed storage."""
    package_root = Path(surface_pd.__file__).parent
    for path in package_root.rglob("*.py"):
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            is_constructor = (
                isinstance(node, ast.FunctionDef) and node.name == "__init__"
            )
            if not is_constructor:
                continue
            parameters = (
                *node.args.posonlyargs,
                *node.args.args,
                *node.args.kwonlyargs,
            )
            private = [
                parameter.arg
                for parameter in parameters
                if parameter.arg.startswith("_")
                and parameter.arg != "self"
            ]
            assert private == [], f"{path}:{node.lineno}: {private}"


def test_enumeration_slab_has_explicit_validated_signature():
    """Surface-specific configuration should be public and keyword-only."""
    signature = inspect.signature(core.EnumerationSlab)
    assert list(signature.parameters) == [
        "lattice",
        "species",
        "coords",
        "charge",
        "validate_proximity",
        "to_unit_cell",
        "coords_are_cartesian",
        "site_properties",
        "labels",
        "properties",
        "direction",
        "layer_tolerance_angstrom",
        "enumerated_species",
        "num_enumerated_layers",
        "symmetric",
    ]
    for name in (
        "direction",
        "layer_tolerance_angstrom",
        "enumerated_species",
        "num_enumerated_layers",
        "symmetric",
    ):
        assert (
            signature.parameters[name].kind
            is inspect.Parameter.KEYWORD_ONLY
        )


@pytest.mark.parametrize(
    "old_name",
    [
        "_direction",
        "_layer_tolerance_angstrom",
        "_to_be_enumerated_species",
        "_num_enumerated_layers",
        "_symmetric",
    ],
)
def test_enumeration_slab_rejects_private_constructor_names(old_name):
    """Former private-looking names should not remain as compatibility API."""
    with pytest.raises(TypeError, match="unexpected keyword argument"):
        EnumerationSlab(
            Lattice.cubic(3),
            ["Li"],
            [[0, 0, 0]],
            **{old_name: None},
        )


@pytest.mark.parametrize(
    ("name", "value"),
    [
        ("direction", True),
        ("direction", -1),
        ("direction", 3),
        ("layer_tolerance_angstrom", True),
        ("layer_tolerance_angstrom", 0),
        ("layer_tolerance_angstrom", float("inf")),
        ("enumerated_species", "Li"),
        ("enumerated_species", []),
        ("enumerated_species", ["Li", "Li"]),
        ("num_enumerated_layers", {}),
        ("num_enumerated_layers", {"Li": 0}),
        ("num_enumerated_layers", {"Li": True}),
        ("symmetric", 1),
    ],
)
def test_enumeration_slab_rejects_invalid_configuration(name, value):
    """Construction and property assignment should share validation."""
    arguments = {
        "lattice": Lattice.cubic(3),
        "species": ["Li"],
        "coords": [[0, 0, 0]],
    }
    with pytest.raises((TypeError, ValueError)):
        core.EnumerationSlab(**arguments, **{name: value})

    slab = core.EnumerationSlab(**arguments)
    with pytest.raises((TypeError, ValueError)):
        setattr(slab, name, value)


def test_enumeration_slab_owns_configuration_collections():
    """Caller mutation should not bypass validated slab configuration."""
    species = ["Li"]
    layers = {"Li": 1}
    slab = core.EnumerationSlab(
        Lattice.cubic(3),
        ["Li"],
        [[0, 0, 0]],
        labels=["surface"],
        properties={"source": "test"},
        direction=1,
        layer_tolerance_angstrom=0.1,
        enumerated_species=species,
        num_enumerated_layers=layers,
        symmetric=False,
    )
    species.append("O")
    layers["O"] = 2

    assert slab.direction == 1
    assert slab.layer_tolerance_angstrom == 0.1
    assert not hasattr(slab, "tolerance")
    assert slab.enumerated_species == ["Li"]
    assert not hasattr(slab, "to_be_enumerated_species")
    assert slab.num_enumerated_layers == {"Li": 1}
    assert slab.symmetric is False
    assert slab.labels == ["surface"]
    assert slab.properties == {"source": "test"}


def test_enumeration_slab_from_structure_preserves_structure_data():
    """The factory should preserve pymatgen data and add surface options."""
    structure = Structure(
        Lattice.tetragonal(3.0, 15.0),
        ["Li", "O"],
        [[0, 0, 0.2], [0, 0, 0.8]],
        charge=0,
        site_properties={"tag": ["bottom", "top"]},
        labels=["li", "oxygen"],
        properties={"source": "unit test"},
    )

    slab = EnumerationSlab.from_structure(
        structure,
        direction=2,
        enumerated_species=["Li"],
        num_enumerated_layers={"Li": 1},
        symmetric=False,
    )

    assert slab == structure
    assert slab.charge == structure.charge
    assert slab.site_properties == structure.site_properties
    assert slab.labels == structure.labels
    assert slab.properties == structure.properties
    assert slab.enumerated_species == ["Li"]


def test_enumeration_slab_from_file_preserves_data_and_surface_options(
    tmp_path,
):
    """The file factory should parse with pymatgen and configure the slab."""
    structure = Structure(
        Lattice.tetragonal(3.0, 15.0),
        ["Li", "O"],
        [[0, 0, 0.2], [0, 0, 0.8]],
        site_properties={"selective_dynamics": [[False] * 3, [True] * 3]},
    )
    filename = tmp_path / "POSCAR"
    structure.to(filename=filename, fmt="poscar")

    slab = EnumerationSlab.from_file(
        filename,
        direction=1,
        layer_tolerance_angstrom=0.05,
        enumerated_species=["Li"],
        num_enumerated_layers={"Li": 1},
        symmetric=False,
    )

    assert isinstance(slab, EnumerationSlab)
    assert slab == Structure.from_file(filename)
    assert slab.direction == 1
    assert slab.layer_tolerance_angstrom == 0.05
    assert slab.enumerated_species == ["Li"]
    assert slab.num_enumerated_layers == {"Li": 1}
    assert slab.symmetric is False


def test_enumeration_slab_from_file_uses_default_surface_options(tmp_path):
    """A file call should retain defaults other than explicit direction."""
    filename = tmp_path / "structure.json"
    Structure(Lattice.cubic(3), ["Li"], [[0, 0, 0]]).to(filename=filename)

    slab = EnumerationSlab.from_file(filename, direction=2)

    assert isinstance(slab, EnumerationSlab)
    assert slab.direction == 2
    assert slab.layer_tolerance_angstrom == 0.5
    assert slab.enumerated_species is None
    assert slab.num_enumerated_layers is None
    assert slab.symmetric is None


def test_high_level_slab_factories_require_direction():
    """Factories must not guess which lattice direction contains vacuum."""
    structure_signature = inspect.signature(EnumerationSlab.from_structure)
    file_signature = inspect.signature(EnumerationSlab.from_file)

    assert (
        structure_signature.parameters["direction"].default
        is inspect.Parameter.empty
    )
    assert (
        file_signature.parameters["direction"].default
        is inspect.Parameter.empty
    )
    with pytest.raises(TypeError, match="direction"):
        EnumerationSlab.from_structure(
            Structure(Lattice.cubic(3), ["Li"], [[0, 0, 0]])
        )


def test_layer_map_keys_must_match_enumerated_species():
    """Each selected species should have exactly one layer count."""
    with pytest.raises(ValueError, match="same species"):
        EnumerationSlab(
            Lattice.cubic(3),
            ["Li", "O"],
            [[0, 0, 0], [0, 0, 0.5]],
            enumerated_species=["Li", "O"],
            num_enumerated_layers={"Li": 1},
        )


def test_removed_post_check_compatibility_wrapper_is_not_available():
    """The internal post-check pipeline has one supported repair method."""
    assert hasattr(PostCheck, "repair_refined_slab_geometry")
    assert not hasattr(PostCheck, "slab_size_check")


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
