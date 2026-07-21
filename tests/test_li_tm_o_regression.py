"""Frozen scientific regressions for the legacy Li--Ni--O workflow."""

import json
from pathlib import Path

import numpy as np
import pytest

from surface_pd.thermodynamics import (
    DatasetAlignment,
    FixedPressureOxygenChemicalPotential,
    GrandPotentialModel,
    IntercalationChemicalPotential,
    Phase,
    PhaseDataset,
    ReferencePhase,
    ThermodynamicState,
)

_FIXTURE = Path(__file__).parent / "data" / "li_tm_o_legacy_regression.json"


@pytest.fixture(scope="module")
def regression():
    """Return the reviewed, implementation-independent regression values."""
    return json.loads(_FIXTURE.read_text())


def _model(regression):
    """Construct the generalized Li--Ni--O model from frozen references."""
    references = regression["references"]
    method = regression["provenance"]["calculation_method"]
    bulk = ReferencePhase(
        "bulk-LiNiO2",
        {"Li": 1, "Ni": 1, "O": 2},
        references["bulk_lino2_ev_per_formula_unit"],
        method,
    )
    model = GrandPotentialModel(
        ("Li", "Ni", "O"),
        {
            "Li": IntercalationChemicalPotential(
                references["li_ev_per_atom"], 1, "Li/Li+"
            ),
            "O": FixedPressureOxygenChemicalPotential(
                references["o2_raw_ev_per_molecule"],
                references["o2_correction_ev_per_molecule"],
            ),
        },
        (bulk,),
    )
    return model, bulk


def test_generalized_chemical_potentials_match_frozen_legacy_values(
    regression,
):
    """The two reservoir laws should preserve the validated specialization."""
    model, _ = _model(regression)
    grid = regression["grid"]
    state = ThermodynamicState(
        {
            "voltage": np.array(grid["voltage_v"]),
            "temperature": np.array(grid["temperature_k"]),
        }
    )

    chemical_potentials = model.chemical_potentials(state)

    np.testing.assert_allclose(
        chemical_potentials["Li"],
        regression["chemical_potentials"]["li_ev_per_atom"],
        rtol=0.0,
        atol=1e-13,
    )
    np.testing.assert_allclose(
        chemical_potentials["O"],
        regression["chemical_potentials"]["o_ev_per_atom"],
        rtol=0.0,
        atol=1e-13,
    )


def test_generalized_surface_energy_matches_representative_legacy_phase(
    regression,
):
    """The generalized grand potential should reproduce the legacy equation."""
    model, _ = _model(regression)
    frozen_phase = regression["representative_phase"]
    lattice = frozen_phase["lattice"]
    area = (
        np.sin(np.deg2rad(lattice["gamma_degrees"]))
        * lattice["a"]
        * lattice["b"]
    )
    phase = Phase(
        frozen_phase["phase_id"],
        frozen_phase["composition"],
        frozen_phase["dft_energy_ev"],
        area,
        frozen_phase["surface_multiplicity"],
    )
    dataset = PhaseDataset(
        "representative",
        ("Li", "Ni", "O"),
        (phase,),
        regression["provenance"]["calculation_method"],
    )
    voltage, temperature = np.meshgrid(
        regression["grid"]["voltage_v"],
        regression["grid"]["temperature_k"],
        indexing="xy",
    )

    result = model.evaluate(
        dataset,
        ThermodynamicState(
            {"voltage": voltage, "temperature": temperature}
        ),
    )

    np.testing.assert_allclose(
        result.surface_grand_potential_ev_per_angstrom2[0],
        frozen_phase["surface_energy_ev_per_angstrom2"],
        rtol=0.0,
        atol=1e-13,
    )


def test_alignment_difference_is_bounded_and_scientifically_explicit(
    regression,
):
    """Area-normalization convention should explain the small legacy delta."""
    _, bulk = _model(regression)
    frozen = regression["lno001_alignment"]
    method = regression["provenance"]["calculation_method"]
    multiplicity = frozen["surface_multiplicity"]

    def dataset(dataset_id, definition):
        phase = Phase(
            definition["phase_id"],
            definition["composition"],
            definition["dft_energy_ev"],
            definition["surface_area_angstrom2"],
            multiplicity,
        )
        return PhaseDataset(
            dataset_id, ("Li", "Ni", "O"), (phase,), method
        )

    reference = dataset("li", frozen["reference_anchor"])
    target = dataset("ni", frozen["target_anchor"])
    alignment = DatasetAlignment(
        reference,
        target,
        frozen["reference_anchor"]["phase_id"],
        frozen["target_anchor"]["phase_id"],
        bulk,
    )
    generalized_density_shift = alignment.energy_offset_ev / (
        multiplicity
        * frozen["target_anchor"]["surface_area_angstrom2"]
    )

    assert alignment.bulk_unit_count == frozen["bulk_unit_count"]
    assert alignment.energy_offset_ev == pytest.approx(
        frozen["generalized_total_offset_ev"], abs=1e-13
    )
    assert generalized_density_shift == pytest.approx(
        frozen["generalized_shift_ev_per_target_angstrom2"], abs=1e-13
    )
    assert abs(
        generalized_density_shift
        - frozen["legacy_shift_ev_per_angstrom2"]
    ) == pytest.approx(
        frozen["maximum_grid_difference_ev_per_angstrom2"], abs=1e-13
    )


def test_stable_phase_baselines_cover_both_maintained_orientations(
    regression,
):
    """Migration must later reproduce both recorded stable-phase maps."""
    expected = regression["stable_phase_ids"]

    assert np.asarray(expected["lno001_charge"]).shape == (3, 3)
    assert np.asarray(expected["lno104_charge"]).shape == (3, 3)
    assert all(
        ":" in phase_id
        for diagram in expected.values()
        for row in diagram
        for phase_id in row
    )
