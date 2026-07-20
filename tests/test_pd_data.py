"""Contract tests for phase-diagram data preparation and alignment."""

import numpy as np
import pandas as pd
import pytest

from surface_pd.plot import PdData, ReferenceEnergies


@pytest.fixture
def references():
    """Return a compact explicit reference set."""
    return ReferenceEnergies("test; U_Ni=? eV", -1.0, -2.0, 0.0, -3.0)


def _dataframe(**overrides):
    """Return two phase rows suitable for standardization tests."""
    values = {
        "Li": [1.0, 4.0],
        "Ni": [1.0, 2.0],
        "O": [2.0, 4.0],
        "E": [-10.0, -20.0],
        "a": [2.0, 2.0],
        "b": [3.0, 3.0],
        "gamma": [90.0, 90.0],
    }
    values.update(overrides)
    return pd.DataFrame(values, index=pd.Index(["s1", "s2"], name="structure"))


def _pd_data(dataframe, references):
    """Construct phase data with standard species names."""
    return PdData(dataframe, "Li", "O", references)


def test_pd_data_owns_a_defensive_dataframe_copy(references):
    """Normalization must not mutate the caller's dataframe."""
    original = _dataframe()
    expected = original.copy(deep=True)

    phase_data = _pd_data(original, references).standardize_pd_data()

    pd.testing.assert_frame_equal(original, expected)
    assert phase_data.dataframe is not original


def test_standardization_is_fluent_idempotent_and_scales_area(references):
    """Extensive values and area should share one TM normalization factor."""
    phase_data = _pd_data(_dataframe(), references)

    result = phase_data.standardize_pd_data()
    once = result.dataframe.copy(deep=True)
    result.standardize_pd_data()

    assert result is phase_data
    pd.testing.assert_frame_equal(result.dataframe, once)
    assert list(result.dataframe["structure"]) == ["s2", "s1"]
    assert list(result.dataframe["Li"]) == [4.0, 2.0]
    assert list(result.dataframe["Ni"]) == [2.0, 2.0]
    assert list(result.dataframe["O"]) == [4.0, 4.0]
    assert list(result.dataframe["dft_energy"]) == [-20.0, -20.0]
    assert list(result.dataframe["a"]) == [2.0, 2.0]
    assert list(result.dataframe["b"]) == [3.0, 6.0]
    assert phase_data.transition_metal == "Ni"


def test_equal_energy_aliases_collapse_to_canonical_column(references):
    """Identical E and dft_energy columns should not remain duplicated."""
    dataframe = _dataframe(dft_energy=[-10.0, -20.0])

    result = _pd_data(dataframe, references).standardize_pd_data().dataframe

    assert "E" not in result
    assert "dft_energy" in result


def test_conflicting_energy_aliases_are_rejected(references):
    """Two different energy columns are scientifically ambiguous."""
    dataframe = _dataframe(dft_energy=[-10.0, -21.0])

    with pytest.raises(ValueError, match="conflicting.*E.*dft_energy"):
        _pd_data(dataframe, references).standardize_pd_data()


@pytest.mark.parametrize(
    ("drop_column", "message"),
    [
        ("Li", "Missing required columns.*Li"),
        ("O", "Missing required columns.*O"),
        ("E", "Missing required columns.*dft_energy"),
        ("a", "Missing required columns.*a"),
        ("b", "Missing required columns.*b"),
        ("gamma", "Missing required columns.*gamma"),
    ],
)
def test_missing_required_columns_are_rejected(
    references,
    drop_column,
    message,
):
    """Incomplete tables should fail during explicit standardization."""
    dataframe = _dataframe().drop(columns=drop_column)

    with pytest.raises(ValueError, match=message):
        _pd_data(dataframe, references).standardize_pd_data()


def test_species_columns_must_be_distinct(references):
    """One column cannot represent both reservoir species."""
    with pytest.raises(ValueError, match="must be distinct"):
        PdData(_dataframe(), "Li", "Li", references)


def test_exactly_one_transition_metal_is_required(references):
    """Element detection should reject absent and ambiguous TM columns."""
    without_tm = _dataframe().drop(columns="Ni")
    with_two = _dataframe(Fe=[1.0, 2.0])

    with pytest.raises(ValueError, match="exactly one transition-metal.*none"):
        _pd_data(without_tm, references).standardize_pd_data()
    with pytest.raises(
        ValueError,
        match="exactly one transition-metal.*Fe, Ni",
    ):
        _pd_data(with_two, references).standardize_pd_data()


@pytest.mark.parametrize(
    ("column", "values", "message"),
    [
        ("Li", [-1.0, 2.0], "Li counts must be nonnegative"),
        ("O", [2.0, np.nan], "O must contain only finite"),
        ("Ni", [0.0, 2.0], "Ni counts must be positive"),
        ("E", [-10.0, np.inf], "dft_energy must contain only finite"),
        ("a", [0.0, 2.0], "a must be positive"),
        ("b", [3.0, -1.0], "b must be positive"),
        ("gamma", [90.0, 180.0], "gamma must be between 0 and 180"),
    ],
)
def test_invalid_numerical_data_is_rejected(
    references,
    column,
    values,
    message,
):
    """Invalid physical or non-finite values must fail before calculation."""
    dataframe = _dataframe(**{column: values})

    with pytest.raises(ValueError, match=message):
        _pd_data(dataframe, references).standardize_pd_data()


def test_alignment_reference_is_deterministic(references):
    """Tied endpoint compositions should select the lowest-energy phase."""
    dataframe = pd.DataFrame(
        {
            "Li": [2, 0, 0],
            "Ni": [1, 1, 1],
            "O": [4, 2, 2],
            "E": [-5.0, -1.0, -2.0],
            "a": [2.0] * 3,
            "b": [3.0] * 3,
            "gamma": [90.0] * 3,
        },
        index=pd.Index(["max", "tie-high", "tie-low"], name="structure"),
    )
    phase_data = _pd_data(dataframe, references).standardize_pd_data()

    reference = phase_data.get_alignment_reference()

    assert reference["structure"] == "tie-low"
    assert reference["dft_energy"] == -2.0


def test_normalized_methods_require_standardization(references):
    """State-dependent operations should fail with a corrective instruction."""
    phase_data = _pd_data(_dataframe(), references)

    with pytest.raises(RuntimeError, match="standardize_pd_data"):
        phase_data.get_alignment_reference()
    with pytest.raises(RuntimeError, match="standardize_pd_data"):
        phase_data.get_surface_energy(np.array([0.0]), np.array([300.0]))


def _alignment_pair(references):
    """Return the complementary committed LNO-001 reference phases."""
    first = pd.DataFrame(
        {
            "Li": [8.0],
            "Ni": [12.0],
            "O": [24.0],
            "E": [-415.07716],
            "a": [5.661],
            "b": [5.802],
            "gamma": [59.156],
        }
    )
    second = pd.DataFrame(
        {
            "Li": [12.0],
            "Ni": [16.0],
            "O": [32.0],
            "E": [-562.40779],
            "a": [5.517],
            "b": [5.098],
            "gamma": [89.994],
        }
    )
    return (
        _pd_data(first, references).standardize_pd_data(),
        _pd_data(second, references).standardize_pd_data(),
    )


def test_shift_energy_uses_bulk_difference_and_first_area():
    """The committed complementary references should reproduce their shift."""
    references = ReferenceEnergies(
        "SCAN+rVV10+U; U_Ni=? eV",
        -2.33333,
        -12.00701,
        0.0,
        -36.8133525,
    )
    first, second = _alignment_pair(references)

    shift = first.calculate_shift_energy(second)

    assert shift == pytest.approx(0.0013691625836902974)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (
            {"references": ReferenceEnergies("other", -1, -2, 0, -3)},
            "reference energies",
        ),
        ({"column": "Ni", "replacement": "Co"}, "transition metal"),
        ({"column": "a", "value": 7.0}, "surface areas"),
        ({"column": "O", "value": 31.0}, "whole LiTMO2"),
        (
            {"column": "dft_energy", "value": 0.0},
            "one alignment energy is zero",
        ),
    ],
)
def test_shift_energy_rejects_incompatible_datasets(
    references,
    mutation,
    message,
):
    """Alignment must reject mismatched scientific and structural inputs."""
    first, second = _alignment_pair(references)
    if "references" in mutation:
        second.reference_energies = mutation["references"]
    elif "replacement" in mutation:
        second.dataframe.rename(
            columns={mutation["column"]: mutation["replacement"]},
            inplace=True,
        )
        second._transition_metal = mutation["replacement"]
    else:
        second.dataframe.loc[0, mutation["column"]] = mutation["value"]

    with pytest.raises(ValueError, match=message):
        first.calculate_shift_energy(second)


def test_two_zero_alignment_energies_return_zero(references):
    """Discharge masks should align without manufacturing an energy offset."""
    first, second = _alignment_pair(references)
    first.dataframe.loc[0, "dft_energy"] = 0.0
    second.dataframe.loc[0, "dft_energy"] = 0.0

    assert first.calculate_shift_energy(second) == 0.0


@pytest.mark.parametrize(
    ("shape", "expected_shape"),
    [((), (2,)), ((3,), (2, 3)), ((2, 3), (2, 2, 3))],
)
def test_surface_energy_preserves_mesh_shape(
    references,
    shape,
    expected_shape,
):
    """The phase axis should be prepended without flattening input meshes."""
    phase_data = _pd_data(_dataframe(), references).standardize_pd_data()
    voltage = np.zeros(shape)
    temperature = np.full(shape, 300.0)

    result = phase_data.get_surface_energy(voltage, temperature)

    assert result.shape == expected_shape


def test_surface_energy_requires_matching_mesh_shapes(references):
    """Voltage and temperature meshes must address identical conditions."""
    phase_data = _pd_data(_dataframe(), references).standardize_pd_data()

    with pytest.raises(ValueError, match="identical shapes"):
        phase_data.get_surface_energy(np.zeros(2), np.zeros((2, 1)))
