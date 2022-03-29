"""
Functions used to plot surface pd
"""

import numpy as np
import collections


__author__ = "Xinhao Li"
__email__ = "xl2778@columbia.edu"
__date__ = "2022-03-15"


E_O2_by_funtional = {"PBE+U": -9.86018 + 1.36,
                       "SCAN+rVV10+U": -12.00701}

def g_oxygen(T, functional="PBE+U"):
    """
    The equation was obtained from J. Osorio-Guillen, S. Lany,
    S. V. Barabash and A. Zunger, Phys. Rev. Lett., 2006, 96, 107203.
    :param T: Temperature
    :param functional: functional used to perform the DFT calculations
    :return: chemical potential of oxygen atom
    """
    E_O2 = E_O2_by_funtional[functional]
    T0 = 298  # K
    kB = 0.008314463  # kJ/(mol K)
    H0 = 8.683  # kJ/mol
    S0 = 205.147 * 1e-3  # kJ/(mol K)
    Cp = 3.5 * kB  # kJ/(mol K)
    delta_H = Cp * (T - T0)  # kJ/mol
    delta_S = Cp * np.log(T / T0)  # kJ/(mol K)
    delta_mu = (H0 + delta_H) - T * (S0 + delta_S)  # kJ/mol
    # print(delta_mu/96.487)
    mu_oxygen = (1 / 2) * (E_O2 + (delta_mu / 96.487))  # eV
    return mu_oxygen

E_Li_by_funtional = {"PBE+U": -1.89965,
                       "SCAN+rVV10+U": -2.33333}

E_bulk_by_funtional = {"PBE+U": -19.92283375,
                       "SCAN+rVV10+U": -36.8133525}

def gibbs_free_energy(V, T,
                      nLi, nNi, nO,
                      dft_energy,
                      a, b, gamma,
                      functional="PBE+U"):
    E_bulk = E_bulk_by_funtional[functional]
    E_Li = E_Li_by_funtional[functional]
    A = np.sin(gamma * np.pi / 180) * a * b
    G = (dft_energy + (nNi - nLi) * (E_Li - V)
         + (2 * nNi - nO) * g_oxygen(T, functional=functional)
         - nNi * E_bulk
         ) / (2 * A)
    return G


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
    converted_unique_phase = np.linspace(0, num_phases-1, num_phases)
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
        converted_uniqe_phase: unique and continuous phase number

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

def get_surface_energy(dataframe,
                       V, T,
                       functional="PBE+U"):
    E = []
    for i in range(len(dataframe)):
        E = np.append(
            E, gibbs_free_energy(
                V, T, *dataframe.iloc[i].to_list()[0:],
                functional=functional
            )
        )
    return E

def find_stable_phases(G, dataframe):

    # Split an array into multiple sub-arrays
    G_splitted = np.split(G, len(dataframe))

    # Stack arrays in sequence horizontally
    G_hstack = np.column_stack(G_splitted)

    # Get the indices of the minimum values along an axis
    stable_phases_index = np.argmin(G_hstack, axis=1)
    return stable_phases_index

def get_compositions(dataframe, species, ticks):
        full = max(dataframe.iloc[:][species])
        empty = min(dataframe.iloc[:][species])
        labels = []
        for i in (ticks):
            comp = ((dataframe.iloc[i][species] - empty) /
                    (full - empty))
            labels.append(str(comp * 100) + "%" + str(species))
        return labels

def get_labels(dataframe, species: list, ticks):
    labels = collections.defaultdict(list)
    for s in species:
        species_label = get_compositions(dataframe, s, ticks)
        labels[s] = species_label

    display = labels[str(species[0])]
    for i in range(1, len(labels)):
        for j in range(len(display)):

            display[j] = str(display[j]) + str(labels[species[i]][j])
    return display


