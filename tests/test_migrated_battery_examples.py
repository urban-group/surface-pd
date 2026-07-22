"""Integration tests for migrated Li--Ni--O example configurations."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

from surface_pd.configuration import PhaseDiagramConfiguration
from surface_pd.thermodynamics import DiagramAxis, PhaseDiagramSpecification

_PROJECT_ROOT = Path(__file__).resolve().parents[1]
_EXAMPLE_ROOT = _PROJECT_ROOT / "examples" / "plotting-examples"
_REGRESSION = json.loads(
    (_PROJECT_ROOT / "tests/data/li_tm_o_legacy_regression.json").read_text()
)
_CONFIGURATIONS = {
    "lno-001-pbe/charge.json": 50,
    "lno-001-scan/charge.json": 50,
    "lno-001-scan/discharge.json": 5,
    "lno-104-pbe/charge.json": 25,
    "lno-104-scan/charge.json": 25,
    "lno-104-scan/discharge.json": 20,
}


@pytest.mark.parametrize(
    ("relative_path", "phase_count"), _CONFIGURATIONS.items()
)
def test_migrated_configuration_loads_explicit_phase_data(
    relative_path, phase_count
):
    """Every retained scientific example should use the generalized format."""
    configuration = PhaseDiagramConfiguration.read_json(
        _EXAMPLE_ROOT / relative_path
    )
    datasets = configuration.load_datasets()

    assert sum(len(dataset.phases) for dataset in datasets) == phase_count
    assert configuration.components == ("Li", "Ni", "O")
    assert configuration.calculation_method.endswith("U_Ni=? eV")
    assert all(
        phase.number_of_surfaces == 2
        and phase.surface_area_angstrom2 > 0
        for dataset in datasets
        for phase in dataset.phases
    )


@pytest.mark.parametrize(
    ("directory", "expected_count"),
    [("lno-001-scan", 5), ("lno-104-scan", 20)],
)
def test_discharge_candidates_preserve_original_charge_energies(
    directory, expected_count
):
    """Filtering should remove candidates without rewriting their energies."""
    charge = PhaseDiagramConfiguration.read_json(
        _EXAMPLE_ROOT / directory / "charge.json"
    ).load_datasets()
    discharge = PhaseDiagramConfiguration.read_json(
        _EXAMPLE_ROOT / directory / "discharge.json"
    ).load_datasets()
    charge_energies = {
        phase.phase_id: phase.dft_energy_ev
        for dataset in charge
        for phase in dataset.phases
    }
    discharge_phases = [
        phase for dataset in discharge for phase in dataset.phases
    ]

    assert len(discharge_phases) == expected_count
    assert all(
        phase.dft_energy_ev == charge_energies[phase.phase_id]
        for phase in discharge_phases
    )
    assert all(phase.dft_energy_ev != 0.0 for phase in discharge_phases)


@pytest.mark.parametrize(
    ("relative_path", "baseline_name"),
    [
        ("lno-001-scan/charge.json", "lno001_charge"),
        ("lno-104-scan/charge.json", "lno104_charge"),
    ],
)
def test_migrated_scan_data_reproduce_frozen_stable_phases(
    relative_path, baseline_name
):
    """Migrated explicit data should preserve the validated phase map."""
    configuration = PhaseDiagramConfiguration.read_json(
        _EXAMPLE_ROOT / relative_path
    )
    specification = PhaseDiagramSpecification(
        DiagramAxis(
            "voltage",
            _REGRESSION["grid"]["voltage_v"],
            "Potential vs. Li/Li+",
            "V",
        ),
        DiagramAxis(
            "temperature",
            _REGRESSION["grid"]["temperature_k"],
            "Temperature",
            "K",
        ),
        {},
    )

    result = specification.evaluate(
        configuration.model, configuration.load_datasets()
    )
    actual = [
        [result.phase_ids[index] for index in row]
        for row in result.representative_phase_indices
    ]

    assert actual == _REGRESSION["stable_phase_ids"][baseline_name]


def test_migrated_tables_have_no_legacy_metadata_or_inferred_geometry():
    """JSON should be the sole thermodynamic source of truth."""
    tables = sorted(_EXAMPLE_ROOT.glob("lno-*/*.dat"))

    assert len(tables) == 8
    for table in tables:
        first_line = table.read_text().splitlines()[0]
        assert not first_line.startswith("#")
        columns = first_line.split()
        assert columns == [
            "phase_id",
            "Li",
            "Ni",
            "O",
            "dft_energy_ev",
            "surface_area_angstrom2",
        ]


def test_migrated_cli_configuration_creates_a_diagram(tmp_path):
    """A maintained battery example should execute through the CLI."""
    output = tmp_path / "lno-104.png"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
                "surface_pd.cli.surface_pd_plot",
                str(_EXAMPLE_ROOT / "lno-104-scan/discharge.json"),
                "--x-range",
                "0",
                "5",
                "--y-range",
                "1",
                "1500",
                "--output",
            str(output),
        ],
        cwd=_PROJECT_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert output.is_file()
    assert output.stat().st_size > 0
