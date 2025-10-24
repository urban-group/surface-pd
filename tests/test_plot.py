"""Tests for plotting utilities."""

import numpy as np
import pandas as pd

from surface_pd.plot.plot import (
    convert_numbers,
    find_stable_phases,
    get_compositions,
    get_ticks_and_levels,
)


class TestFindStablePhases:
    """Tests for find_stable_phases function."""

    def test_simple_case(self):
        """Test with simple energy values."""
        # Create test data: 3 phases, 5 composition points
        G = np.array(
            [
                1.0,
                2.0,
                3.0,
                4.0,
                5.0,  # Phase 0
                2.0,
                1.5,
                2.5,
                3.5,
                4.5,  # Phase 1
                3.0,
                2.5,
                2.0,
                2.5,
                3.0,  # Phase 2
            ]
        )
        df = pd.DataFrame({"phase": [0, 1, 2]})

        result = find_stable_phases(G, df)

        # Phase 0 is stable at point 0 (1.0 < 2.0 < 3.0)
        # Phase 1 is stable at point 1 (1.5 < 2.0 < 2.5)
        # Phase 2 is stable at points 2, 3, 4 (2.0 < 2.5 < 3.0, etc.)
        expected = np.array([0, 1, 2, 2, 2])
        np.testing.assert_array_equal(result, expected)

    def test_single_phase(self):
        """Test with single dominant phase."""
        G = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
        df = pd.DataFrame({"phase": [0]})

        result = find_stable_phases(G, df)

        expected = np.array([0, 0, 0, 0, 0])
        np.testing.assert_array_equal(result, expected)


class TestConvertNumbers:
    """Tests for convert_numbers function."""

    def test_continuous_sequence(self):
        """Test with continuous sequence."""
        array = np.array([0, 0, 1, 1, 2, 2])
        unique_phases = np.array([0, 1, 2])

        result = convert_numbers(array, unique_phases)

        expected = np.array([0, 0, 1, 1, 2, 2])
        np.testing.assert_array_equal(result, expected)

    def test_discontinuous_sequence(self):
        """Test with discontinuous sequence."""
        array = np.array([0, 0, 2, 2, 5, 5])
        unique_phases = np.array([0, 2, 5])

        result = convert_numbers(array, unique_phases)

        # Should map [0, 2, 5] -> [0, 1, 2]
        expected = np.array([0, 0, 1, 1, 2, 2])
        np.testing.assert_array_equal(result, expected)

    def test_single_phase(self):
        """Test with single phase."""
        array = np.array([5, 5, 5])
        unique_phases = np.array([5])

        result = convert_numbers(array, unique_phases)

        expected = np.array([0, 0, 0])
        np.testing.assert_array_equal(result, expected)


class TestGetTicksAndLevels:
    """Tests for get_ticks_and_levels function."""

    def test_simple_phases(self):
        """Test with simple phase array."""
        converted_phase = np.array([0, 0, 1, 1, 2, 2])

        ticks, levels = get_ticks_and_levels(converted_phase)

        # Should have ticks for 3 phases (0, 1, 2)
        assert len(ticks) == 3
        # Levels should be from -1 to 2 (min-1 to max)
        assert len(levels) == 4
        np.testing.assert_array_almost_equal(
            levels, np.array([-1.0, 0.0, 1.0, 2.0])
        )

    def test_single_phase(self):
        """Test with single phase."""
        converted_phase = np.array([0, 0, 0])

        ticks, levels = get_ticks_and_levels(converted_phase)

        assert len(ticks) == 1
        assert len(levels) == 2


class TestGetCompositions:
    """Tests for get_compositions function."""

    def test_single_file(self):
        """Test with single data file."""
        df = pd.DataFrame(
            {
                "Li": [1, 2, 3],
                "O": [5, 5, 5],
                "Ni": [5, 5, 5],
            }
        )
        ticks = [0, 1, 2]

        labels = get_compositions(df, 1, "Li", ticks)

        # Should show Li compositions as percentages
        assert len(labels) == 3
        # Check that labels are formatted correctly
        for label in labels:
            assert "%" in label
            assert "Li" in label

    def test_multiple_files(self):
        """Test with multiple data files."""
        df = pd.DataFrame(
            {
                "Li": [2, 4, 6, 8, 10, 12],
                "O": [5, 5, 5, 5, 5, 5],
                "Ni": [5, 5, 5, 5, 5, 5],
            }
        )
        ticks = [0, 1, 2]

        labels = get_compositions(df, 2, "Li", ticks)

        # Should normalize by num_files
        assert len(labels) == 3
        # Check that labels are formatted correctly
        for label in labels:
            assert "%" in label
            assert "Li" in label
