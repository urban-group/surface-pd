"""Tests for explicit, non-mutating slab-dataset alignment."""

from pathlib import Path

import numpy as np
import pytest

from surface_pd.configuration import PhaseDiagramConfiguration
from surface_pd.thermodynamics import (
    AlignedPhaseDataset,
    ConstantChemicalPotential,
    DatasetAlignment,
    GrandPotentialModel,
    Phase,
    PhaseDataset,
    ReferencePhase,
    ThermodynamicState,
)

_METHOD = "test DFT method"
def _phase(
    phase_id, composition, energy, area=2.0, number_of_surfaces=2
):
    """Return one concise phase for alignment tests."""
    return Phase(phase_id, composition, energy, area, number_of_surfaces)


def _dataset(dataset_id, phases, components=("A", "B"), method=_METHOD):
    """Return a dataset using the standard test basis and provenance."""
    return PhaseDataset(dataset_id, components, phases, method)


def _bulk(composition=None, energy=15.0, method=_METHOD):
    """Return the standard bulk reference used by synthetic tests."""
    return ReferencePhase(
        "bulk-A2B",
        composition or {"A": 2, "B": 1},
        energy,
        method,
    )


def _alignment_pair():
    """Return anchors where target minus reference is one A2B unit."""
    reference = _dataset(
        "reference",
        (
            _phase("anchor", {"A": 1, "B": 1}, 10.0),
            _phase("other", {"A": 2, "B": 1}, 12.0),
        ),
    )
    target = _dataset(
        "target",
        (
            _phase("anchor", {"A": 3, "B": 2}, 30.0),
            _phase("other", {"A": 4, "B": 2}, 33.0),
        ),
    )
    return reference, target


def test_alignment_uses_explicit_anchor_equation_and_sign():
    """The target offset should satisfy the documented anchor equality."""
    reference, target = _alignment_pair()

    alignment = DatasetAlignment(
        reference,
        target,
        "anchor",
        "anchor",
        _bulk(),
    )

    assert alignment.bulk_unit_count == 1
    assert alignment.energy_offset_ev == pytest.approx(-5.0)
    assert alignment.reference_anchor_id == "reference:anchor"
    assert alignment.target_anchor_id == "target:anchor"
    target_anchor = target.get_phase("anchor")
    reference_anchor = reference.get_phase("anchor")
    assert target_anchor.dft_energy_ev + alignment.energy_offset_ev == (
        reference_anchor.dft_energy_ev
        + alignment.bulk_unit_count
        * alignment.bulk_reference.energy_ev_per_formula_unit
    )


@pytest.mark.parametrize(
    ("reverse", "target_energy", "expected_units", "expected_offset"),
    [
        (True, None, -1, 5.0),
        (False, 11.0, 0, -1.0),
    ],
)
def test_alignment_supports_negative_and_zero_bulk_unit_counts(
    reverse, target_energy, expected_units, expected_offset
):
    """The integer bulk-unit count should retain its physical sign."""
    reference, target = _alignment_pair()
    if reverse:
        reference, target = target, reference
    else:
        target = _dataset(
            "same-composition",
            (_phase("anchor", {"A": 1, "B": 1}, target_energy),),
        )

    alignment = DatasetAlignment(
        reference, target, "anchor", "anchor", _bulk()
    )

    assert alignment.bulk_unit_count == expected_units
    assert alignment.energy_offset_ev == pytest.approx(expected_offset)


def test_aligned_view_preserves_source_data_and_alignment_provenance():
    """Alignment should be an immutable energy view, not a data rewrite."""
    reference, target = _alignment_pair()
    alignment = DatasetAlignment(
        reference, target, "anchor", "anchor", _bulk()
    )

    aligned = alignment.create_aligned_dataset()

    assert isinstance(aligned, AlignedPhaseDataset)
    assert aligned.source_dataset is target
    assert aligned.alignment is alignment
    assert aligned.dataset_id == "target"
    assert aligned.root_dataset_id == "reference"
    assert aligned.components == target.components
    assert aligned.phases == target.phases
    assert aligned.qualified_phase_ids == target.qualified_phase_ids
    assert aligned.energy_offset_ev == pytest.approx(-5.0)
    assert target.get_phase("anchor").dft_energy_ev == 30.0


def test_grand_potential_model_applies_offset_only_to_aligned_view():
    """Every target phase should receive one state-independent total offset."""
    reference, target = _alignment_pair()
    aligned = DatasetAlignment(
        reference, target, "anchor", "anchor", _bulk()
    ).create_aligned_dataset()
    model = GrandPotentialModel(
        ("A", "B"),
        {
            "A": ConstantChemicalPotential(2.0),
            "B": ConstantChemicalPotential(3.0),
        },
        (),
    )
    state = ThermodynamicState({})

    raw_result = model.evaluate(target, state)
    aligned_result = model.evaluate(aligned, state)

    np.testing.assert_allclose(
        aligned_result.total_grand_potential_ev
        - raw_result.total_grand_potential_ev,
        [-5.0, -5.0],
    )
    assert target.get_phase("other").dft_energy_ev == 33.0


def test_alignment_rejects_missing_anchor_and_self_alignment():
    """Anchor identities and alignment direction should be unambiguous."""
    reference, target = _alignment_pair()
    with pytest.raises(ValueError, match="reference_anchor_phase_id"):
        DatasetAlignment(reference, target, "missing", "anchor", _bulk())
    with pytest.raises(ValueError, match="target_anchor_phase_id"):
        DatasetAlignment(reference, target, "anchor", "missing", _bulk())
    with pytest.raises(ValueError, match="distinct datasets"):
        DatasetAlignment(reference, reference, "anchor", "anchor", _bulk())


@pytest.mark.parametrize(
    ("target", "bulk", "message"),
    [
        (
            _dataset(
                "target",
                (_phase("anchor", {"B": 2, "A": 3}, 30.0),),
                components=("B", "A"),
            ),
            _bulk(),
            "component bases",
        ),
        (
            _dataset(
                "target",
                (_phase("anchor", {"A": 3, "B": 2}, 30.0),),
                method="other method",
            ),
            _bulk(),
            "calculation_method",
        ),
        (
            _dataset(
                "target",
                (_phase("anchor", {"A": 3, "B": 2}, 30.0),),
            ),
            _bulk(method="other method"),
            "bulk reference.*calculation_method",
        ),
        (
            _dataset(
                "target",
                (_phase("anchor", {"A": 3, "B": 2}, 30.0, area=2.1),),
            ),
            _bulk(),
            "surface areas",
        ),
        (
            _dataset(
                "target",
                (
                    _phase(
                        "anchor",
                        {"A": 3, "B": 2},
                        30.0,
                        number_of_surfaces=1,
                    ),
                ),
            ),
            _bulk(),
            "numbers of surfaces",
        ),
        (
            _dataset(
                "target",
                (_phase("anchor", {"A": 4, "B": 2}, 30.0),),
            ),
            _bulk(),
            "integer number.*bulk reference",
        ),
        (
            _dataset(
                "target",
                (_phase("anchor", {"A": 3, "B": 2}, 30.0),),
            ),
            _bulk(composition={"A": 2, "C": 1}),
            "bulk reference.*component basis",
        ),
    ],
)
def test_alignment_rejects_incompatible_physical_inputs(
    target, bulk, message
):
    """Every scientific compatibility assumption should be checked."""
    reference, _ = _alignment_pair()

    with pytest.raises(ValueError, match=message):
        DatasetAlignment(reference, target, "anchor", "anchor", bulk)


def test_alignment_rejects_chaining_an_aligned_view():
    """The initial API should permit only direct-to-root alignments."""
    reference, target = _alignment_pair()
    aligned = DatasetAlignment(
        reference, target, "anchor", "anchor", _bulk()
    ).create_aligned_dataset()

    with pytest.raises(TypeError, match="ordinary PhaseDataset"):
        DatasetAlignment(reference, aligned, "anchor", "anchor", _bulk())


def test_migrated_lno001_configuration_preserves_explicit_alignment():
    """The maintained LNO-001 input should retain reviewed provenance."""
    path = (
        Path(__file__).resolve().parents[1]
        / "examples/plotting-examples/lno-001-scan/charge.json"
    )
    _, aligned = PhaseDiagramConfiguration.read_json(path).load_datasets()

    assert aligned.alignment.bulk_unit_count == 4
    assert aligned.alignment.energy_offset_ev == pytest.approx(0.07722)
    assert aligned.alignment.reference_anchor_id == (
        "li:Li8Ni12O24_3-phase-25"
    )
    assert aligned.alignment.target_anchor_id == "ni:full-O-full-Li"
