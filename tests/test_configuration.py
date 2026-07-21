"""Tests for versioned generalized phase-diagram configuration."""

import json
from pathlib import Path

import pytest

from surface_pd.configuration import PhaseDiagramConfiguration
from surface_pd.thermodynamics import (
    ConstantChemicalPotential,
    DirectChemicalPotential,
    FixedPressureOxygenChemicalPotential,
    IntercalationChemicalPotential,
)

_METHOD = "test DFT method; U_X=? eV"
_PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _configuration_data():
    """Return a valid configuration containing every built-in model."""
    return {
        "schema_version": 1,
        "calculation_method": _METHOD,
        "components": ["Li", "O", "X", "C"],
        "independent_chemical_potentials": {
            "Li": {
                "model": "intercalation_voltage",
                "reference_energy_ev_per_component": -2.0,
                "electrons_per_component": 1,
                "voltage_reference": "Li/Li+",
            },
            "O": {
                "model": "fixed_pressure_oxygen",
                "raw_o2_energy_ev_per_molecule": -10.0,
                "correction_ev_per_molecule": 1.0,
            },
            "X": {
                "model": "direct",
                "state_variable": "delta_mu",
                "reference_energy_ev_per_component": -1.0,
            },
            "C": {
                "model": "constant",
                "value_ev_per_component": 0.5,
            },
        },
        "reference_phases": [],
        "diagram": {
            "x_axis": {
                "state_variable": "voltage",
                "coordinates": {
                    "kind": "linear",
                    "start": 0.0,
                    "stop": 5.0,
                    "number": 6,
                },
                "label": "Potential vs. Li/Li+",
                "unit": "V",
            },
            "y_axis": {
                "state_variable": "temperature",
                "coordinates": {
                    "kind": "values",
                    "values": [298.15, 500.0, 1000.0],
                },
                "label": "Temperature",
                "unit": "K",
            },
            "fixed_conditions": {"delta_mu": -0.25},
        },
        "datasets": [
            {
                "dataset_id": "dataset",
                "path": "phases.dat",
                "number_of_surfaces": 2,
                "column_overrides": {
                    "phase_id": "structure",
                    "composition": {
                        "Li": "n_Li",
                        "O": "n_O",
                        "X": "n_X",
                        "C": "n_C",
                    },
                    "dft_energy_ev": "energy",
                    "surface_area_angstrom2": "area",
                },
            }
        ],
        "alignments": [],
        "rendering": {
            "coloring": {"mode": "phase_identity"},
            "colormap": "tab20",
            "invert_x_axis": False,
            "invert_y_axis": False,
        },
    }


def test_configuration_constructs_every_builtin_model_and_diagram():
    """Registry names should map only to the reviewed package models."""
    configuration = PhaseDiagramConfiguration(_configuration_data())

    models = configuration.model.independent_chemical_potentials
    assert isinstance(models["Li"], IntercalationChemicalPotential)
    assert isinstance(models["O"], FixedPressureOxygenChemicalPotential)
    assert isinstance(models["X"], DirectChemicalPotential)
    assert isinstance(models["C"], ConstantChemicalPotential)
    assert configuration.schema_version == 1
    assert configuration.calculation_method == _METHOD
    assert configuration.components == ("Li", "O", "X", "C")
    assert configuration.diagram_specification.x_axis.values.tolist() == [
        0.0,
        1.0,
        2.0,
        3.0,
        4.0,
        5.0,
    ]
    assert configuration.diagram_specification.y_axis.values.tolist() == [
        298.15,
        500.0,
        1000.0,
    ]


def test_configuration_constructs_dependent_reference_constraints():
    """Reference configuration should retain stoichiometry and provenance."""
    data = _configuration_data()
    del data["independent_chemical_potentials"]["C"]
    data["reference_phases"] = [
        {
            "reference_id": "bulk-CX",
            "composition": {"C": 1, "X": 1},
            "energy_ev_per_formula_unit": -3.0,
        }
    ]

    configuration = PhaseDiagramConfiguration(data)

    reference = configuration.model.reference_phases[0]
    assert reference.reference_id == "bulk-CX"
    assert dict(reference.composition) == {"C": 1, "X": 1}
    assert reference.calculation_method == _METHOD
    assert configuration.model.dependent_components == ("C",)


def test_configuration_json_round_trip_is_canonical_and_owned(tmp_path):
    """Read/write should preserve canonical data without shared mutations."""
    source = tmp_path / "input.json"
    source.write_text(json.dumps(_configuration_data()))

    configuration = PhaseDiagramConfiguration.read_json(source)
    copied = configuration.to_dict()
    copied["components"].append("mutated")
    output = tmp_path / "output.json"
    configuration.write_json(output)
    reread = PhaseDiagramConfiguration.read_json(output)

    assert configuration.source_path == source.resolve()
    assert reread.to_dict() == configuration.to_dict()
    assert configuration.components == ("Li", "O", "X", "C")
    assert output.read_text().endswith("\n")


@pytest.mark.parametrize(
    ("path", "value", "message"),
    [
        (("unknown",), 1, "configuration.*unknown fields"),
        (
            ("independent_chemical_potentials", "Li", "unknown"),
            1,
            "intercalation_voltage.*unknown fields",
        ),
        (
            ("diagram", "x_axis", "coordinates", "unknown"),
            1,
            "linear coordinates.*unknown fields",
        ),
        (
            ("datasets", 0, "column_overrides", "unknown"),
            "column",
            "column_overrides.*unknown fields",
        ),
        (
            ("rendering", "coloring", "unknown"),
            1,
            "phase_identity coloring.*unknown fields",
        ),
    ],
)
def test_configuration_rejects_unknown_fields_recursively(
    path, value, message
):
    """Typos should fail in their local schema context."""
    data = _configuration_data()
    target = data
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value

    with pytest.raises(ValueError, match=message):
        PhaseDiagramConfiguration(data)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda data: data.pop("components"), "configuration.*missing fields"),
        (
            lambda data: data.update(schema_version=2),
            "unsupported schema_version",
        ),
        (
            lambda data: data["independent_chemical_potentials"]["Li"].update(
                model="python:arbitrary.callable"
            ),
            "unknown chemical-potential model",
        ),
        (
            lambda data: data["diagram"]["x_axis"]["coordinates"].update(
                number=1
            ),
            "number",
        ),
        (
            lambda data: data["datasets"][0].update(number_of_surfaces=0),
            "number_of_surfaces",
        ),
    ],
)
def test_configuration_rejects_missing_versioned_or_unsafe_data(
    mutation, message
):
    """Unsupported or executable configuration should never be accepted."""
    data = _configuration_data()
    mutation(data)

    with pytest.raises((TypeError, ValueError), match=message):
        PhaseDiagramConfiguration(data)


def test_configuration_rejects_malformed_json(tmp_path):
    """JSON syntax errors should identify the configuration path."""
    path = tmp_path / "broken.json"
    path.write_text("{not JSON")

    with pytest.raises(ValueError, match="broken.json.*valid JSON"):
        PhaseDiagramConfiguration.read_json(path)


def test_version_one_json_schema_is_committed_and_parseable():
    """Users and editors should have a machine-readable schema document."""
    schema_path = (
        _PROJECT_ROOT
        / "surface_pd"
        / "schemas"
        / "phase-diagram-config-v1.schema.json"
    )

    schema = json.loads(schema_path.read_text())

    assert schema["$schema"] == "https://json-schema.org/draft/2020-12/schema"
    assert schema["properties"]["schema_version"] == {"const": 1}
    assert schema["additionalProperties"] is False
