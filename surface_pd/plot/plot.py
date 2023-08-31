import collections

import numpy as np
import pandas as pd


def find_stable_phases(G: list,
                       dataframe: pd.DataFrame):
    """
    1. Split the list into sub-arrays first, i.e.
        [1, 2, 3, 4, 5, ..., 12]
        -->
        [[1, 2, 3, 4],
         [5, 6, 7, 8],
         [9, 10, 11, 12]]
    2. Stack arrays in sequence horizontally,
        [[1, 2, 3, 4],
         [5, 6, 7, 8],
         [9, 10, 11, 12]],
         -->
        [[1, 5, 9],
         [2, 6, 10],
         [3, 7, 11],
         [4, 8, 12]]
    3. Find the indices of the minimum values along an axis

    Args:
        G: list of phases
        dataframe:

    Returns:
        Indices of the most stable phase of each composition
    """

    # Split an array into multiple sub-arrays
    G_split = np.split(G, len(dataframe))

    # Stack arrays in sequence horizontally
    G_hstack = np.column_stack(G_split)

    # Get the indices of the minimum values along an axis
    stable_phases_index = np.argmin(G_hstack, axis=1)
    return stable_phases_index


def convert_numbers(array, unique_phases):
    """
    Convert unique phases into a continuous phase number
    i.e. [1, 4, 6, 15] --> [0, 1, 2, 3]

    Args:
        array:
        unique_phases:

    Returns:

    """
    array_copy = array.copy()
    num_phases = len(unique_phases)
    converted_unique_phase = np.linspace(0, num_phases - 1, num_phases)
    for i in range(len(unique_phases)):
        for j in range(len(array)):
            if array[j] == unique_phases[i]:
                array_copy[j] = converted_unique_phase[i]
    return array_copy


def get_ticks_and_levels(converted_phase):
    """
    Get the positions of ticks which will be labeled beside the colorbar.
    Get unique level of each phase which is used to make sure that each
    phase corresponds to only one color in the color bar.

    Args:
        converted_phase:

    Returns:
        Positions of the ticks and unique levels
    """
    converted_unique_phase = np.unique(converted_phase)

    position = [x - 0.5 for x in converted_unique_phase]

    min_num = min(converted_unique_phase)
    max_num = max(converted_unique_phase)
    unique_levels = np.linspace(min_num - 1, max_num,
                                len(converted_unique_phase) + 1)
    return position, unique_levels


def get_compositions(dataframe, num_files, species, ticks):
    """

    Args:
        dataframe:
        num_files:
        species:
        ticks:

    Returns:

    """
    total_relaxed = (max(dataframe.iloc[:]["O"]) -
                     min(dataframe.iloc[:]["O"])) * num_files
    boundary = int(len(dataframe) / num_files)
    labels = []
    for tick in ticks:
        if tick < boundary:
            comp = ((dataframe.iloc[tick][species] -
                     min(dataframe.iloc[0:boundary][species])) /
                    total_relaxed)
            if num_files != 1:
                comp += 0.5
            labels.append(
                str(round(comp * 100, 1)) + "%" + str(species))
        else:
            comp = ((dataframe.iloc[tick][species] -
                     min(dataframe.iloc[boundary:][species])) /
                    total_relaxed)
            labels.append(
                str(round(comp * 100, 1)) + "%" + str(species))
    return labels


def get_labels(dataframe: pd.DataFrame,
               num_files,
               species: list,
               ticks):
    """

    Args:
        dataframe:
        num_files:
        species:
        ticks:

    Returns:

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
