"""
Plotting utilities for surface phase diagram construction.

This module provides functions for analyzing phase stability, converting phase
numbers, and preparing data for visualization of surface phase diagrams.
"""

import collections

import numpy as np
import pandas as pd


def find_stable_phases(G: list, dataframe: pd.DataFrame) -> np.ndarray:
    """
    Identify the most stable phases for each composition from free energy data.

    This function processes a list of Gibbs free energies to determine which
    phase is most stable at each composition. The algorithm:
    1. Splits the free energy list into sub-arrays (one per phase)
    2. Stacks these arrays horizontally for comparison
    3. Finds the minimum value (most stable phase) along each composition

    Example transformation:
        Input: [1, 2, 3, 4, 5, ..., 12] (energies for 3 phases, 4 compositions)
        Step 1: [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]]
        Step 2: [[1, 5, 9], [2, 6, 10], [3, 7, 11], [4, 8, 12]]
        Step 3: Find argmin along axis 1

    Args:
        G: List of Gibbs free energies for all phases and compositions.
            Length should be (number of phases) × (number of compositions).
        dataframe: DataFrame containing composition and energy data for each phase.

    Returns
    -------
        Array of indices indicating the most stable phase at each composition point.
    """
    # Split an array into multiple sub-arrays
    G_split = np.split(G, len(dataframe))

    # Stack arrays in sequence horizontally
    G_hstack = np.column_stack(G_split)

    # Get the indices of the minimum values along an axis
    stable_phases_index = np.argmin(G_hstack, axis=1)
    return stable_phases_index


def convert_numbers(
    array: np.ndarray, unique_phases: np.ndarray
) -> np.ndarray:
    """
    Convert discontinuous phase indices to continuous sequential numbers.

    This function remaps phase indices to a continuous range starting from 0,
    which is useful for consistent colorbar mapping in phase diagrams.

    Example:
        [1, 4, 6, 15] → [0, 1, 2, 3]

    Args:
        array: Array of phase indices that may have gaps in numbering.
        unique_phases: Array of unique phase numbers present in the data.

    Returns
    -------
        Array with phase indices converted to continuous sequential numbers
        starting from 0.
    """
    array_copy = array.copy()
    num_phases = len(unique_phases)
    converted_unique_phase = np.linspace(0, num_phases - 1, num_phases)
    for i in range(len(unique_phases)):
        for j in range(len(array)):
            if array[j] == unique_phases[i]:
                array_copy[j] = converted_unique_phase[i]
    return array_copy


def get_ticks_and_levels(
    converted_phase: np.ndarray,
) -> tuple[list, np.ndarray]:
    """
    Calculate tick positions and color levels for phase diagram colorbar.

    This function determines the appropriate positions for colorbar ticks and
    creates discrete levels to ensure each phase has a unique color in the
    phase diagram visualization.

    Args:
        converted_phase: Array of converted (continuous) phase indices.

    Returns
    -------
        tuple: A tuple containing:
            - position: List of tick positions offset by -0.5 for centering
            - unique_levels: Array of discrete levels for colorbar boundaries,
              spanning from (min-1) to max with (n_phases+1) levels
    """
    converted_unique_phase = np.unique(converted_phase)

    position = [x - 0.5 for x in converted_unique_phase]

    min_num = min(converted_unique_phase)
    max_num = max(converted_unique_phase)
    unique_levels = np.linspace(
        min_num - 1, max_num, len(converted_unique_phase) + 1
    )
    return position, unique_levels


def get_compositions(
    dataframe: pd.DataFrame, num_files: int, species: str, ticks: list
) -> list[str]:
    """
    Calculate composition percentages for colorbar labels.

    This function computes the percentage composition of a specified species
    at each tick position for labeling the phase diagram colorbar.

    Args:
        dataframe: DataFrame containing composition data with species columns.
        num_files: Number of data files processed (affects composition calculation).
        species: Chemical species name (e.g., "Li", "O", "Ni") for which to
            calculate compositions.
        ticks: List of tick positions on the colorbar.

    Returns
    -------
        List of formatted strings showing percentage composition of the species,
        e.g., ["25.0%Li", "50.0%Li", "75.0%Li"].
    """
    total_relaxed = (
        max(dataframe.iloc[:]["O"]) - min(dataframe.iloc[:]["O"])
    ) * num_files
    boundary = int(len(dataframe) / num_files)
    labels = []
    for tick in ticks:
        if tick < boundary:
            comp = (
                dataframe.iloc[tick][species]
                - min(dataframe.iloc[0:boundary][species])
            ) / total_relaxed
            if num_files != 1:
                comp += 0.5
            labels.append(str(round(comp * 100, 1)) + "%" + str(species))
        else:
            comp = (
                dataframe.iloc[tick][species]
                - min(dataframe.iloc[boundary:][species])
            ) / total_relaxed
            labels.append(str(round(comp * 100, 1)) + "%" + str(species))
    return labels


def get_labels(
    dataframe: pd.DataFrame, num_files: int, species: list, ticks: list
) -> list[str]:
    """
    Generate multi-species composition labels for phase diagram colorbar.

    This function creates formatted labels showing the composition of multiple
    species at each tick position, combining individual species percentages
    into single label strings.

    Args:
        dataframe: DataFrame containing composition data for all species.
        num_files: Number of data files processed (affects composition calculation).
        species: List of chemical species names (e.g., ["Li", "O"]) to include
            in the labels.
        ticks: List of tick positions on the colorbar.

    Returns
    -------
        List of formatted strings combining all species compositions,
        e.g., ["25.0%Li 50.0%O", "50.0%Li 75.0%O"].
    """
    labels = collections.defaultdict(list)
    for s in species:
        species_label = get_compositions(dataframe, num_files, s, ticks)
        labels[s] = species_label

    display = labels[str(species[0])]
    for i in range(1, len(labels)):
        for j in range(len(display)):
            display[j] = str(display[j]) + " " + str(labels[species[i]][j])
    return display
