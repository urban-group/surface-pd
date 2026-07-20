"""Tests for self-contained phase-diagram input files."""

from io import StringIO
from pathlib import Path

import pandas as pd
import pytest

from surface_pd.plot import ReferenceEnergies
from surface_pd.plot._phase_data_io import (
    _read_phase_diagram_file,
    _require_matching_reference_energies,
    _write_phase_diagram_file,
)

VALID_INPUT = """\
# method = SCAN+rVV10+U; U_Ni=? eV
# reference_li_ev_per_atom = -2.33333
# reference_o2_raw_ev_per_molecule = -12.00701
# reference_o2_correction_ev_per_molecule = 0.0
# reference_bulk_litmo2_ev_per_formula_unit = -36.8133525
structure Li Ni O E a b gamma
phase-1 1 1 2 -40.0 2.0 3.0 90.0
"""
PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_phase_diagram_loader_returns_data_and_references(tmp_path):
    """One file should provide both its phase table and scientific inputs."""
    path = tmp_path / "phases.dat"
    path.write_text(VALID_INPUT)

    dataframe, references = _read_phase_diagram_file(path)

    assert list(dataframe.columns) == ["Li", "Ni", "O", "E", "a", "b", "gamma"]
    assert dataframe.loc["phase-1", "E"] == -40.0
    assert references == ReferenceEnergies(
        method="SCAN+rVV10+U; U_Ni=? eV",
        li_ev_per_atom=-2.33333,
        o2_raw_ev_per_molecule=-12.00701,
        o2_correction_ev_per_molecule=0.0,
        bulk_litmo2_ev_per_formula_unit=-36.8133525,
    )


@pytest.mark.parametrize(
    ("old", "new", "message"),
    [
        (
            "# reference_li_ev_per_atom = -2.33333\n",
            "",
            "Missing required metadata.*reference_li_ev_per_atom",
        ),
        (
            "# reference_li_ev_per_atom = -2.33333\n",
            "# reference_li_ev_per_atom = -2.33333\n" * 2,
            "Duplicate metadata key.*reference_li_ev_per_atom",
        ),
        (
            "# reference_li_ev_per_atom = -2.33333",
            "# reference_li_ev_per_atom = unknown",
            "reference_li_ev_per_atom.*finite number",
        ),
        (
            "# reference_li_ev_per_atom = -2.33333",
            "# reference_li_ev_per_atom = nan",
            "li_ev_per_atom.*finite real number",
        ),
        (
            "# method = SCAN+rVV10+U; U_Ni=? eV",
            "# method =   ",
            "method must not be empty",
        ),
        (
            "# method = SCAN+rVV10+U; U_Ni=? eV",
            "# typo_reference = -1.0",
            "Unknown metadata key.*typo_reference",
        ),
    ],
)
def test_phase_diagram_loader_rejects_invalid_metadata(
    tmp_path,
    old,
    new,
    message,
):
    """Malformed scientific metadata should fail before dataframe use."""
    path = tmp_path / "invalid.dat"
    path.write_text(VALID_INPUT.replace(old, new))

    with pytest.raises(ValueError, match=message):
        _read_phase_diagram_file(path)


def test_phase_diagram_loader_rejects_metadata_after_header(tmp_path):
    """Metadata is valid only in one leading block."""
    path = tmp_path / "misplaced.dat"
    path.write_text(VALID_INPUT + "# method = misplaced\n")

    with pytest.raises(ValueError, match="Metadata must precede"):
        _read_phase_diagram_file(path)


def test_phase_diagram_loader_rejects_empty_table(tmp_path):
    """A metadata-only file cannot define a phase diagram."""
    path = tmp_path / "empty.dat"
    path.write_text(VALID_INPUT.rsplit("phase-1", maxsplit=1)[0])

    with pytest.raises(ValueError, match="at least one phase row"):
        _read_phase_diagram_file(path)


def test_phase_diagram_writer_round_trips_metadata_and_table(tmp_path):
    """Generated discharge files should retain complete provenance."""
    dataframe = pd.read_csv(
        StringIO("structure Li Ni O E a b gamma\np 1 1 2 0 2 3 90\n"),
        sep=r"\s+",
        index_col=0,
    )
    references = ReferenceEnergies(
        "custom; U_Ni=? eV", -1.0, -2.0, 0.25, -3.0
    )
    path = tmp_path / "written.dat"

    _write_phase_diagram_file(path, dataframe, references)
    actual_dataframe, actual_references = _read_phase_diagram_file(path)

    pd.testing.assert_frame_equal(actual_dataframe, dataframe)
    assert actual_references == references


def test_reference_compatibility_reports_every_differing_field():
    """Multi-file errors should identify rather than conceal mismatches."""
    first = ReferenceEnergies("method one", -1.0, -2.0, 0.0, -3.0)
    second = ReferenceEnergies("method two", -1.1, -2.0, 0.0, -3.0)

    with pytest.raises(
        ValueError,
        match="differing fields: method, li_ev_per_atom",
    ):
        _require_matching_reference_energies(first, second)


def test_identical_reference_sets_are_compatible():
    """Exact equality should permit aligned multi-file calculations."""
    references = ReferenceEnergies("method", -1.0, -2.0, 0.0, -3.0)

    _require_matching_reference_energies(references, references)


def test_all_committed_phase_examples_are_self_contained():
    """Every maintained example should carry explicit, auditable references."""
    example_paths = sorted(
        path
        for path in (
            PROJECT_ROOT / "examples" / "plotting-examples"
        ).glob("**/*.dat")
        if ".ipynb_checkpoints" not in path.parts
    )

    assert len(example_paths) == 12
    for path in example_paths:
        dataframe, references = _read_phase_diagram_file(path)
        assert not dataframe.empty
        assert "U_Ni=? eV" in references.method
