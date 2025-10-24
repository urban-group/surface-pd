"""Tests for utility functions."""

import pytest

from surface_pd.error import NonIntegerError, TooLargeSlabError
from surface_pd.util import (
    all_int,
    check_int,
    csv2dict,
    define_scaling_matrix,
    get_values_nested_dict,
    have_zero,
    replace_dummy,
)


class TestDefineScalingMatrix:
    """Tests for define_scaling_matrix function."""

    def test_multiple_4(self):
        """Test with maximum cell size of 4."""
        result = define_scaling_matrix(10.0, 10.0, 4)
        assert result == [2, 2, 1]

    def test_multiple_3_ratio_2(self):
        """Test with cell size 3 and a/b ratio of 2."""
        result = define_scaling_matrix(20.0, 10.0, 3)
        assert result == [1, 3, 1]

    def test_multiple_3_ratio_half(self):
        """Test with cell size 3 and b/a ratio of 2."""
        result = define_scaling_matrix(10.0, 20.0, 3)
        assert result == [3, 1, 1]

    def test_multiple_3_normal(self):
        """Test with cell size 3 and normal ratio."""
        result = define_scaling_matrix(10.0, 12.0, 3)
        assert result == [3, 1, 1]

    def test_too_large_error(self):
        """Test that TooLargeSlabError is raised for multiple > 4."""
        with pytest.raises(TooLargeSlabError):
            define_scaling_matrix(10.0, 10.0, 5)


class TestCsv2dict:
    """Tests for csv2dict function."""

    def test_simple_substitution(self):
        """Test simple key-value substitution."""
        result = csv2dict(["Li='Na'"])
        assert result == {"Li": "Na"}

    def test_numeric_value(self):
        """Test numeric value conversion."""
        result = csv2dict(["Co=0.5"])
        assert result == {"Co": 0.5}

    def test_nested_dictionary(self):
        """Test nested dictionary creation."""
        result = csv2dict(["Co=Co:0.5&Ni:0.5"])
        assert result == {"Co": {"Co": 0.5, "Ni": 0.5}}

    def test_multiple_entries(self):
        """Test multiple entries."""
        result = csv2dict(["Co=Co:0.5&Ni:0.5", "Li='Na'"])
        assert result == {"Co": {"Co": 0.5, "Ni": 0.5}, "Li": "Na"}

    def test_expression_evaluation(self):
        """Test that pandas eval works for expressions."""
        result = csv2dict(["x=1+1"])
        assert result == {"x": 2}


class TestCheckInt:
    """Tests for check_int function."""

    def test_valid_integer(self):
        """Test valid integer within tolerance."""
        assert check_int(5.0) == 5.0
        assert check_int(5.001) == 5.0
        assert check_int(4.999) == 5.0

    def test_invalid_integer(self):
        """Test invalid integer outside tolerance."""
        with pytest.raises(NonIntegerError):
            check_int(5.02)
        with pytest.raises(NonIntegerError):
            check_int(4.97)


class TestGetValuesNestedDict:
    """Tests for get_values_nested_dict function."""

    def test_flat_dict(self):
        """Test with flat dictionary."""
        result = list(get_values_nested_dict({"a": 1, "b": 2}))
        assert set(result) == {1, 2}

    def test_nested_dict(self):
        """Test with nested dictionary."""
        result = list(get_values_nested_dict({"a": 1, "b": {"c": 2, "d": 3}}))
        assert set(result) == {1, 2, 3}

    def test_deeply_nested(self):
        """Test with deeply nested dictionary."""
        result = list(
            get_values_nested_dict(
                {"a": {"b": {"c": 1}}, "d": {"e": 2}, "f": 3}
            )
        )
        assert set(result) == {1, 2, 3}


class TestAllInt:
    """Tests for all_int function."""

    def test_all_integers(self):
        """Test list with all integers."""
        assert all_int([1.0, 2.0, 3.0]) is True

    def test_mixed_values(self):
        """Test list with mixed integers and non-integers."""
        assert all_int([1.0, 2.5, 3.0]) is False

    def test_all_floats(self):
        """Test list with all non-integers."""
        assert all_int([1.5, 2.5, 3.5]) is False


class TestHaveZero:
    """Tests for have_zero function."""

    def test_with_zero(self):
        """Test list containing zero."""
        assert have_zero([1, 0, 2]) is True
        assert have_zero([0]) is True

    def test_without_zero(self):
        """Test list without zero."""
        assert have_zero([1, 2, 3]) is False

    def test_empty_list(self):
        """Test empty list."""
        assert have_zero([]) is False


class TestReplaceDummy:
    """Tests for replace_dummy function."""

    def test_basic_replacement(self):
        """Test basic dummy species replacement."""
        subs_dict = {"Li": {"Li": 1.0}, "O": {"O": 1.0}}
        dummy_species = ["X", "Y"]
        result = replace_dummy(subs_dict, dummy_species)
        assert "X" in result
        assert "Y" in result
        assert "Li" not in result
        assert "O" not in result

    def test_zero_occupancy_removal(self):
        """Test that zero occupancy entries are removed."""
        subs_dict = {"Li": {"Li": 0.0}, "O": {"O": 1.0}}
        dummy_species = ["X", "Y"]
        result = replace_dummy(subs_dict, dummy_species)
        assert "X" not in result
        assert "Y" in result

    def test_nested_values(self):
        """Test with nested dictionary values."""
        subs_dict = {"Co": {"Co": 0.5, "Ni": 0.5}}
        dummy_species = ["X"]
        result = replace_dummy(subs_dict, dummy_species)
        assert "X" in result
        assert result["X"]["X"] == 0.5
        assert result["X"]["Ni"] == 0.5
