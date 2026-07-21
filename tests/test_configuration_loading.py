"""Tests for explicitly mapped version-1 phase-table loading."""

import json

import pytest

from surface_pd.configuration import PhaseDiagramConfiguration
from surface_pd.thermodynamics import AlignedPhaseDataset, PhaseDataset


def _configuration_data(dataset_path="phases.dat"):
    """Return a minimal configuration with two variable reservoirs."""
    return {
        "schema_version": 1,
        "calculation_method": "test DFT method",
        "components": ["A", "B", "C"],
        "independent_chemical_potentials": {
            "A": {
                "model": "direct",
                "state_variable": "mu_A",
                "reference_energy_ev_per_component": 0.0,
            },
            "B": {
                "model": "direct",
                "state_variable": "mu_B",
                "reference_energy_ev_per_component": 0.0,
            },
        },
        "reference_phases": [
            {
                "reference_id": "bulk-C",
                "composition": {"A": 0, "B": 0, "C": 1},
                "energy_ev_per_formula_unit": -5.0,
            }
        ],
        "diagram": {
            "x_axis": {
                "state_variable": "mu_A",
                "coordinates": {
                    "kind": "values",
                    "values": [-1.0, 0.0],
                },
                "label": "A chemical potential",
                "unit": "eV",
            },
            "y_axis": {
                "state_variable": "mu_B",
                "coordinates": {
                    "kind": "values",
                    "values": [-1.0, 0.0],
                },
                "label": "B chemical potential",
                "unit": "eV",
            },
            "fixed_conditions": {},
        },
        "datasets": [
            {
                "dataset_id": "root",
                "path": dataset_path,
                "columns": {
                    "phase_id": "name",
                    "composition": {"A": "nA", "B": "nB", "C": "nC"},
                    "dft_energy_ev": "energy",
                    "surface_area_angstrom2": {"constant": 10.0},
                    "surface_multiplicity": {"constant": 2},
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


def _write_configuration(tmp_path, data):
    """Write and read one configuration so relative paths have a base."""
    path = tmp_path / "config.json"
    path.write_text(json.dumps(data))
    return PhaseDiagramConfiguration.read_json(path)


def test_load_datasets_maps_columns_constants_and_ignores_extra_data(tmp_path):
    """Only declared columns should acquire thermodynamic meaning."""
    (tmp_path / "phases.dat").write_text(
        "name nA nB nC energy note\n"
        "alpha 1 2 3 -10.5 ignored\n"
        "beta 2.0 1 0 -8.0 also-ignored\n"
    )
    configuration = _write_configuration(tmp_path, _configuration_data())

    datasets = configuration.load_datasets()

    assert len(datasets) == 1
    assert isinstance(datasets[0], PhaseDataset)
    assert datasets[0].dataset_id == "root"
    assert datasets[0].components == ("A", "B", "C")
    assert datasets[0].calculation_method == "test DFT method"
    alpha = datasets[0].get_phase("alpha")
    assert dict(alpha.composition) == {"A": 1, "B": 2, "C": 3}
    assert alpha.dft_energy_ev == -10.5
    assert alpha.surface_area_angstrom2 == 10.0
    assert alpha.surface_multiplicity == 2


def test_load_datasets_accepts_area_and_multiplicity_columns(tmp_path):
    """Per-row surface geometry should be available through explicit maps."""
    data = _configuration_data()
    columns = data["datasets"][0]["columns"]
    columns["surface_area_angstrom2"] = "area"
    columns["surface_multiplicity"] = "surfaces"
    (tmp_path / "phases.dat").write_text(
        "name nA nB nC energy area surfaces\n"
        "alpha 1 2 3 -10.5 12.25 1\n"
    )

    dataset = _write_configuration(tmp_path, data).load_datasets()[0]

    phase = dataset.get_phase("alpha")
    assert phase.surface_area_angstrom2 == 12.25
    assert phase.surface_multiplicity == 1


def test_relative_paths_move_with_the_configuration_directory(tmp_path):
    """Relative table paths should not depend on the working directory."""
    project = tmp_path / "movable"
    project.mkdir()
    (project / "phases.dat").write_text(
        "name nA nB nC energy\nalpha 1 2 3 -10.5\n"
    )

    configuration = _write_configuration(project, _configuration_data())

    phase = configuration.load_datasets()[0].get_phase("alpha")
    assert phase.dft_energy_ev == -10.5


def test_in_memory_configuration_requires_absolute_dataset_paths(tmp_path):
    """Relative paths without a source JSON directory should be ambiguous."""
    configuration = PhaseDiagramConfiguration(_configuration_data())

    with pytest.raises(ValueError, match="relative.*source JSON"):
        configuration.load_datasets()


def test_in_memory_configuration_accepts_absolute_dataset_paths(tmp_path):
    """An absolute table path should not require source-file provenance."""
    path = tmp_path / "phases.dat"
    path.write_text("name nA nB nC energy\nalpha 1 2 3 -10.5\n")
    configuration = PhaseDiagramConfiguration(_configuration_data(str(path)))

    assert configuration.load_datasets()[0].dataset_id == "root"


@pytest.mark.parametrize(
    ("table", "message"),
    [
        (
            "name nA nB nC\nalpha 1 2 3\n",
            "root.*missing required columns.*energy",
        ),
        (
            "name nA nB nC energy\nalpha 1 2 3\n",
            "row 2.*5 fields",
        ),
        (
            "name nA nB nC energy\nalpha 1 bad 3 -1\n",
            "root.*row 2.*nB.*integer",
        ),
        (
            "name nA nB nC energy\nalpha 1 -1 3 -1\n",
            "root.*row 2.*nB.*nonnegative integer",
        ),
        (
            "name nA nB nC energy\nalpha 1 1.5 3 -1\n",
            "root.*row 2.*nB.*nonnegative integer",
        ),
        (
            "name nA nB nC energy\nalpha 1 2 3 not-a-number\n",
            "root.*row 2.*energy.*finite",
        ),
        (
            "name nA nB nC energy\nalpha 1 2 3 -1\nalpha 2 2 3 -2\n",
            "duplicate phase_id.*alpha",
        ),
        (
            "name nA nB nC energy energy\nalpha 1 2 3 -1 -1\n",
            "header columns must be unique",
        ),
        (
            "# legacy metadata = forbidden\n"
            "name nA nB nC energy\n"
            "alpha 1 2 3 -1\n",
            "comments.*not supported",
        ),
        ("name nA nB nC energy\n", "at least one phase row"),
    ],
)
def test_load_datasets_rejects_malformed_or_invalid_tables(
    tmp_path, table, message
):
    """Input failures should identify their dataset and row context."""
    (tmp_path / "phases.dat").write_text(table)
    configuration = _write_configuration(tmp_path, _configuration_data())

    with pytest.raises((TypeError, ValueError), match=message):
        configuration.load_datasets()


def test_configuration_rejects_reused_mapped_source_columns():
    """One source column should not silently acquire two meanings."""
    data = _configuration_data()
    data["datasets"][0]["columns"]["dft_energy_ev"] = "nA"

    with pytest.raises(ValueError, match="mapped source columns.*unique"):
        PhaseDiagramConfiguration(data)


def test_load_datasets_reports_an_unreadable_table_path(tmp_path):
    """Missing files should identify both the dataset and resolved path."""
    configuration = _write_configuration(tmp_path, _configuration_data())

    with pytest.raises(ValueError, match="root.*could not read.*phases.dat"):
        configuration.load_datasets()


@pytest.mark.parametrize(
    ("area", "multiplicity", "message"),
    [
        ("not-a-number", "2", "area.*finite"),
        ("0", "2", "surface_area_angstrom2 must be positive"),
        ("10", "1.5", "surfaces.*positive integer"),
        ("10", "0", "surfaces.*positive integer"),
    ],
)
def test_load_datasets_validates_mapped_surface_geometry(
    tmp_path, area, multiplicity, message
):
    """Area and multiplicity columns should retain their physical domains."""
    data = _configuration_data()
    columns = data["datasets"][0]["columns"]
    columns["surface_area_angstrom2"] = "area"
    columns["surface_multiplicity"] = "surfaces"
    (tmp_path / "phases.dat").write_text(
        "name nA nB nC energy area surfaces\n"
        f"alpha 1 2 3 -1 {area} {multiplicity}\n"
    )

    with pytest.raises(ValueError, match=message):
        _write_configuration(tmp_path, data).load_datasets()


def test_load_datasets_applies_direct_alignment_in_declared_order(tmp_path):
    """Configured targets should become provenance-preserving aligned views."""
    (tmp_path / "root.dat").write_text(
        "name nA nB nC energy\nanchor 1 1 1 10\nother 2 1 1 12\n"
    )
    (tmp_path / "target.dat").write_text(
        "name nA nB nC energy\nanchor 1 1 2 30\nother 2 1 2 33\n"
    )
    data = _configuration_data("root.dat")
    target = json.loads(json.dumps(data["datasets"][0]))
    target.update(dataset_id="target", path="target.dat")
    data["datasets"].append(target)
    data["alignments"] = [
        {
            "root_dataset_id": "root",
            "target_dataset_id": "target",
            "reference_anchor_phase_id": "anchor",
            "target_anchor_phase_id": "anchor",
            "bulk_reference_id": "bulk-C",
        }
    ]

    root, aligned = _write_configuration(tmp_path, data).load_datasets()

    assert isinstance(root, PhaseDataset)
    assert isinstance(aligned, AlignedPhaseDataset)
    assert aligned.dataset_id == "target"
    assert aligned.root_dataset_id == "root"
    assert aligned.energy_offset_ev == pytest.approx(-25.0)
    assert aligned.alignment.reference_anchor_id == "root:anchor"
    assert aligned.alignment.target_anchor_id == "target:anchor"
    assert aligned.source_dataset.get_phase("anchor").dft_energy_ev == 30.0


def test_configuration_rejects_alignment_chains():
    """Every alignment must point directly from an unaligned root."""
    data = _configuration_data()
    for dataset_id in ("middle", "target"):
        dataset = json.loads(json.dumps(data["datasets"][0]))
        dataset["dataset_id"] = dataset_id
        data["datasets"].append(dataset)
    data["alignments"] = [
        {
            "root_dataset_id": "root",
            "target_dataset_id": "middle",
            "reference_anchor_phase_id": "anchor",
            "target_anchor_phase_id": "anchor",
            "bulk_reference_id": "bulk-C",
        },
        {
            "root_dataset_id": "middle",
            "target_dataset_id": "target",
            "reference_anchor_phase_id": "anchor",
            "target_anchor_phase_id": "anchor",
            "bulk_reference_id": "bulk-C",
        },
    ]

    with pytest.raises(ValueError, match="root.*alignment target"):
        PhaseDiagramConfiguration(data)
