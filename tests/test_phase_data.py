"""Tests for immutable generalized thermodynamic phase data."""

from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from surface_pd.thermodynamics import Phase, PhaseDataset, ReferencePhase


def _phase(**overrides):
    """Return a valid phase with selected fields overridden."""
    values = {
        "phase_id": "phase-1",
        "composition": {"Li": 4, "Ni": 4, "O": 8},
        "dft_energy_ev": -80.0,
        "surface_area_angstrom2": 9.0,
        "surface_multiplicity": 2,
    }
    values.update(overrides)
    return Phase(**values)


def _reference(**overrides):
    """Return a valid bulk reference with selected fields overridden."""
    values = {
        "reference_id": "bulk-LiNiO2",
        "composition": {"Li": 1, "Ni": 1, "O": 2},
        "energy_ev_per_formula_unit": -20.0,
        "calculation_method": "test DFT method; U_Ni=? eV",
    }
    values.update(overrides)
    return ReferencePhase(**values)


def test_phase_owns_immutable_absolute_composition():
    """Phase construction should own counts without normalizing them."""
    composition = {"Li": 4.0, "Ni": np.int64(4), "O": 8}

    phase = _phase(composition=composition)
    composition["Li"] = 999

    assert dict(phase.composition) == {"Li": 4, "Ni": 4, "O": 8}
    assert phase.dft_energy_ev == -80.0
    assert phase.surface_area_angstrom2 == 9.0
    assert phase.surface_multiplicity == 2
    with pytest.raises(TypeError):
        phase.composition["Li"] = 3
    with pytest.raises(FrozenInstanceError):
        phase.dft_energy_ev = 0.0


@pytest.mark.parametrize(
    "phase_id",
    ["", "   ", "dataset:phase", "phase\nother"],
)
def test_phase_rejects_invalid_identifier(phase_id):
    """Local phase IDs should compose unambiguously with dataset IDs."""
    with pytest.raises(ValueError, match="phase_id"):
        _phase(phase_id=phase_id)


@pytest.mark.parametrize(
    "composition",
    [
        {},
        {"Li": 0, "Ni": 0, "O": 0},
        {"Li": -1, "Ni": 1, "O": 2},
        {"Li": 0.5, "Ni": 1, "O": 2},
        {"Li": True, "Ni": 1, "O": 2},
        {"Li": np.nan, "Ni": 1, "O": 2},
        {"Li": 1, "Ni": 1, "O atom": 2},
    ],
)
def test_phase_rejects_invalid_composition(composition):
    """Absolute elemental counts should be integral and physically defined."""
    with pytest.raises((TypeError, ValueError), match="composition"):
        _phase(composition=composition)


@pytest.mark.parametrize("energy", [True, "-80", np.nan, np.inf])
def test_phase_rejects_invalid_dft_energy(energy):
    """DFT energy should be a finite real value in eV."""
    with pytest.raises((TypeError, ValueError), match="dft_energy_ev"):
        _phase(dft_energy_ev=energy)


@pytest.mark.parametrize("area", [True, "9", 0.0, -1.0, np.nan, np.inf])
def test_phase_rejects_invalid_surface_area(area):
    """Surface-cell area should be finite and positive."""
    with pytest.raises(
        (TypeError, ValueError), match="surface_area_angstrom2"
    ):
        _phase(surface_area_angstrom2=area)


@pytest.mark.parametrize("multiplicity", [True, 0, -1, 1.5])
def test_phase_rejects_invalid_surface_multiplicity(multiplicity):
    """Surface multiplicity should be a positive integer."""
    with pytest.raises(
        (TypeError, ValueError), match="surface_multiplicity"
    ):
        _phase(surface_multiplicity=multiplicity)


def test_reference_phase_owns_stoichiometry_and_provenance():
    """A bulk equality should retain immutable formula-unit data."""
    composition = {"Li": 1, "Ni": 1, "O": 2}

    reference = _reference(composition=composition)
    composition["O"] = 999

    assert dict(reference.composition) == {"Li": 1, "Ni": 1, "O": 2}
    assert reference.energy_ev_per_formula_unit == -20.0
    assert reference.calculation_method == "test DFT method; U_Ni=? eV"
    with pytest.raises(TypeError):
        reference.composition["O"] = 3


@pytest.mark.parametrize(
    "reference_id",
    ["", "   ", "bulk:reference", "bulk\nreference"],
)
def test_reference_phase_rejects_invalid_identifier(reference_id):
    """Reference identifiers should be nonempty and unambiguous."""
    with pytest.raises(ValueError, match="reference_id"):
        _reference(reference_id=reference_id)


@pytest.mark.parametrize("method", ["", "   ", "PBE\nSCAN", 1])
def test_reference_phase_rejects_invalid_method_provenance(method):
    """Calculation provenance should be nonempty, single-line free text."""
    with pytest.raises((TypeError, ValueError), match="calculation_method"):
        _reference(calculation_method=method)


def test_phase_dataset_owns_phases_and_component_order():
    """A dataset should preserve its basis and deterministic identities."""
    components = ["Li", "Ni", "O"]
    phases = [_phase(), _phase(phase_id="phase-2", dft_energy_ev=-79.0)]

    dataset = PhaseDataset(
        dataset_id="li-terminated",
        components=components,
        phases=phases,
        calculation_method="test DFT method; U_Ni=? eV",
    )
    components.append("H")
    phases.clear()

    assert dataset.components == ("Li", "Ni", "O")
    assert tuple(phase.phase_id for phase in dataset.phases) == (
        "phase-1",
        "phase-2",
    )
    assert dataset.qualified_phase_ids == (
        "li-terminated:phase-1",
        "li-terminated:phase-2",
    )
    assert dataset.get_phase("phase-2").dft_energy_ev == -79.0
    assert dataset.qualified_phase_id("phase-2") == "li-terminated:phase-2"
    with pytest.raises(KeyError, match="missing"):
        dataset.get_phase("missing")


@pytest.mark.parametrize(
    "dataset_id",
    ["", "   ", "dataset:child", "dataset\nchild"],
)
def test_phase_dataset_rejects_invalid_identifier(dataset_id):
    """Dataset IDs should form unambiguous composite phase identities."""
    with pytest.raises(ValueError, match="dataset_id"):
        PhaseDataset(
            dataset_id=dataset_id,
            components=("Li", "Ni", "O"),
            phases=(_phase(),),
            calculation_method="test method",
        )


@pytest.mark.parametrize(
    "components",
    [(), ("Li", "Li"), ("Li", "O atom"), "Li"],
)
def test_phase_dataset_rejects_invalid_component_basis(components):
    """The component basis should be explicit, ordered, and unique."""
    with pytest.raises((TypeError, ValueError), match="components"):
        PhaseDataset(
            dataset_id="dataset",
            components=components,
            phases=(_phase(),),
            calculation_method="test method",
        )


def test_phase_dataset_requires_nonempty_unique_compatible_phases():
    """Every phase should occur once and match the declared component basis."""
    with pytest.raises(ValueError, match="at least one phase"):
        PhaseDataset(
            dataset_id="dataset",
            components=("Li", "Ni", "O"),
            phases=(),
            calculation_method="test method",
        )

    with pytest.raises(ValueError, match="duplicate phase_id"):
        PhaseDataset(
            dataset_id="dataset",
            components=("Li", "Ni", "O"),
            phases=(_phase(), _phase()),
            calculation_method="test method",
        )

    with pytest.raises(ValueError, match="component basis"):
        PhaseDataset(
            dataset_id="dataset",
            components=("Li", "Ni", "O"),
            phases=(_phase(composition={"Li": 1, "Ni": 1, "H": 1}),),
            calculation_method="test method",
        )
