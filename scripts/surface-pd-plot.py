#!/usr/bin/env python

"""
This code will be used to construct the surface phase diagram.
"""

__author__ = "Xinhao Li"
__email__ = "xl2778@columbia.edu"
__date__ = "2022-03-22"

import argparse

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
rcparams = {"font.size": 16,
            "legend.frameon": False,
            "xtick.top": True,
            "xtick.direction": "in",
            "xtick.minor.visible": True,
            "xtick.major.size": 16,
            "xtick.minor.size": 8,
            "ytick.right": True,
            "ytick.direction": "in",
            "ytick.minor.visible": True,
            "ytick.major.size": 16,
            "ytick.minor.size": 8}

from surface_pd.surface_plot import (convert_numbers, get_ticks_and_levels,
                                     get_surface_energy, find_stable_phases,
                                     get_labels)

def surface_pd_plot(data,
                    color_Li=False,
                    color_O=False,
                    functional="PBE+U"):
    # Read data file
    df = pd.read_table(data, sep="\s+", index_col=0)

    # Create x and y axis
    V = np.linspace(0.0, 5.0, 1000)
    T = np.linspace(1.0, 2000.0, 400)
    V_mesh, T_mesh = np.meshgrid(V, T)

    # Get surface energy
    G = get_surface_energy(df, V_mesh, T_mesh, functional=functional)

    # Split an array into multiple sub-arrays
    stable_phases_index = find_stable_phases(G, df)
    ticks = np.unique([stable_phases_index])
    converted_stable_phases_index = convert_numbers(stable_phases_index, ticks)
    #
    stable_phases_reshaped_global = np.reshape(converted_stable_phases_index,
                                        (T.size, V.size))

    if color_Li:
        stable_phases_index = stable_phases_index % 5

        ticks = np.unique([stable_phases_index])
        converted_stable_phases_index = convert_numbers(stable_phases_index,
                                                        ticks)

        stable_phases_reshaped = np.reshape(converted_stable_phases_index,
                                            (T.size, V.size))

        # plot used
        ticky, levels = get_ticks_and_levels(converted_stable_phases_index)
        labels = get_labels(df, species=["Li"], ticks=ticks)
        cmap = "Greens_r"
        fig = plt.figure(figsize=(9.5, 6))
    elif color_O:
        stable_phases_index = stable_phases_index // 5 * 5

        ticks = np.unique([stable_phases_index])
        converted_stable_phases_index = convert_numbers(stable_phases_index,
                                                        ticks)

        stable_phases_reshaped = np.reshape(converted_stable_phases_index,
                                            (T.size, V.size))

        # plot used
        ticky, levels = get_ticks_and_levels(converted_stable_phases_index)
        labels = get_labels(df, species=["O"], ticks=ticks)
        cmap = "Reds_r"
        fig = plt.figure(figsize=(9.5, 6))
    else:
        stable_phases_reshaped = stable_phases_reshaped_global
        # plot used
        ticky, levels = get_ticks_and_levels(converted_stable_phases_index)
        labels = get_labels(df, species=["Li", "O"], ticks=ticks)
        # print(labels)
        cmap = 'tab20c'
        fig = plt.figure(figsize=(10.5, 6))

    plt.rcParams.update(rcparams)
    ax = fig.add_subplot(111)
    ax.contour(V, T, stable_phases_reshaped_global,
               colors="black", linewidths=1.5)
    PD = ax.contourf(V, T, stable_phases_reshaped,
                     levels=levels, cmap=cmap)
    ax.set_xlabel('Potential vs. Li/Li$^+$ (V)')
    ax.set_ylabel('Temperature (K)')
    colorbar = fig.colorbar(PD, ticks=ticky, pad=0.05)
    colorbar.ax.set_yticklabels(labels)
    colorbar.ax.tick_params(size=0)
    colorbar.ax.minorticks_off()
    plt.tight_layout()
    plt.show()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__ + "\n{}{}".format(__date__, __author__),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "surface_pd_data",
        help="Path to the surface pd data."
    )

    parser.add_argument(
        "--color-by-Li-content", "-cl",
        help="Whether color the surface pd by the Li content.",
        action="store_true"
    )

    parser.add_argument(
        "--color-by-O-content", "-co",
        help="Whether color the surface pd by the O content.",
        action="store_true"
    )

    parser.add_argument(
        "--functional", "-f",
        help="Functional used to perform the calculations.",
        type=str,
        default="PBE+U"
    )

    args = parser.parse_args()

    surface_pd_plot(args.surface_pd_data,
                    args.color_by_Li_content,
                    args.color_by_O_content,
                    args.functional
                    )