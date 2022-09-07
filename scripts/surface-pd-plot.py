#!/usr/bin/env python

"""
This code will be used to construct the surface phase diagram.
"""

__author__ = "Xinhao Li"
__email__ = "xl2778@columbia.edu"
__date__ = "2022-08-05"

import argparse

from matplotlib import pyplot as plt
from surface_pd.plot.pd_data import PdData
from surface_pd.plot.plot import *

rcparams = {"font.size": 20,
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


def surface_pd_plot(data_files,
                    lithium_like_species,
                    oxygen_like_species,
                    low_T, high_T,
                    color_Li=False,
                    color_O=False,
                    functional="PBE+U",
                    discharge=False,
                    save=False):
    # Create x and y axes
    V = np.linspace(0.0, 5.0, 500)
    T = np.linspace(low_T, high_T, 500)
    V_mesh, T_mesh = np.meshgrid(V, T)

    # Get the number of files
    num_files = len(data_files)

    if num_files != 1:
        checked_phases = []
        df = []
        for data in data_files:
            temp_df = pd.read_table(data, sep="\s+", index_col=0)

            # Initialize the PdData class and do standardization
            phase_data = PdData(dataframe=temp_df,
                                lithium_like_species=lithium_like_species,
                                oxygen_like_species=oxygen_like_species,
                                functional=functional)

            # Standardize the surface pd data (with the same number of TM
            # species)
            temp_df = phase_data.standardize_pd_data().dataframe

            # Get the phases that will be used to calculate the shift energy
            checked_phases.append(phase_data.get_check_phases())
            df.append(temp_df)

        # Combine all dataframes
        df = pd.concat(df)

        # Get surface energy
        phase_data = PdData(dataframe=df,
                            lithium_like_species=lithium_like_species,
                            oxygen_like_species=oxygen_like_species,
                            functional=functional)

        # Calculate the shift energy for the polar surface only
        shift_energy = phase_data.get_the_shift_energy(
            checked_phases=checked_phases)

        temp_G = phase_data.get_surface_energy(V=V_mesh, T=T_mesh)

        shifted_G = temp_G[int(len(temp_G) / num_files):] + shift_energy
        G = np.append(temp_G[0:int(len(temp_G) / num_files)], shifted_G)
    else:
        # Read data file
        df = pd.read_table(data_files[0], sep="\s+", index_col=0)
        phase_data = PdData(dataframe=df,
                            lithium_like_species=lithium_like_species,
                            oxygen_like_species=oxygen_like_species,
                            functional=functional)
        phase_data.standardize_pd_data()

        # Get surface energy
        G = phase_data.get_surface_energy(V=V_mesh, T=T_mesh)

    # Split an array into multiple sub-arrays
    stable_phases_index = find_stable_phases(G, df)

    # Get unique tick labels
    ticks = np.unique([stable_phases_index])

    # Convert unique phases into a continuous phase number
    converted_stable_phases_index = convert_numbers(stable_phases_index, ticks)

    # Colorbar related parameters
    ticky, global_levels = get_ticks_and_levels(converted_stable_phases_index)

    # Reshape the array
    stable_phases_reshaped_global = np.reshape(converted_stable_phases_index,
                                               (T.size, V.size))

    if color_Li:
        if len(data_files) != 1:
            boundary = int(len(df) / len(data_files))
            for i in range(len(stable_phases_index)):
                if stable_phases_index[i] < boundary:
                    stable_phases_index[i] %= 5
                else:
                    stable_phases_index[i] = stable_phases_index[i] % 5 + 25
        else:
            stable_phases_index = stable_phases_index % 5

        ticks = np.unique([stable_phases_index])

        converted_stable_phases_index = convert_numbers(stable_phases_index,
                                                        ticks)
        stable_phases_reshaped = np.reshape(converted_stable_phases_index,
                                            (T.size, V.size))

        # temp_plot used
        ticky, levels = get_ticks_and_levels(converted_stable_phases_index)
        labels = get_labels(dataframe=df, num_files=num_files,
                            species=[lithium_like_species],
                            ticks=ticks)
        cmap = "Greens_r"
        fig = plt.figure(figsize=(9.5, 6))

    elif color_O:
        # if len(data_files) != 1:
        #     boundary = int(len(df) / len(data_files))
        #     for i in range(len(stable_phases_index)):
        #         if stable_phases_index[i] < boundary:
        #             stable_phases_index[i] = stable_phases_index[i] // 5 * 5
        #         else:
        #             stable_phases_index[i] = (stable_phases_index[i] // 5 *
        #                                       5 + 25)
        # else:
        #     stable_phases_index = stable_phases_index % 5

        stable_phases_index = stable_phases_index // 5 * 5
        # print(np.unique(stable_phases_index))
        stable_phases_index = [x - 5 if x == 25 else x for x in
                               stable_phases_index]

        ticks = np.unique([stable_phases_index])
        # print(ticks)
        converted_stable_phases_index = convert_numbers(stable_phases_index,
                                                        ticks)

        stable_phases_reshaped = np.reshape(converted_stable_phases_index,
                                            (T.size, V.size))

        # temp_plot used
        ticky, levels = get_ticks_and_levels(converted_stable_phases_index)
        labels = get_labels(dataframe=df, num_files=num_files,
                            species=[oxygen_like_species], ticks=ticks)
        cmap = "Reds_r"
        fig = plt.figure(figsize=(9.5, 6))
    else:
        stable_phases_reshaped = stable_phases_reshaped_global
        # temp_plot used
        ticky, levels = get_ticks_and_levels(converted_stable_phases_index)
        labels = get_labels(dataframe=df, num_files=num_files,
                            species=[lithium_like_species,
                                     oxygen_like_species],
                            ticks=ticks)
        cmap = 'tab20c'
        fig = plt.figure(figsize=(10.5, 6))

    plt.rcParams.update(rcparams)
    ax = fig.add_subplot(111)
    # ax = plt.axes(projection='3d')
    ax.contour(V, T, stable_phases_reshaped_global, levels=global_levels,
               colors="black", linewidths=2)
    PD = ax.contourf(V, T, stable_phases_reshaped,
                     levels=levels, cmap=cmap)
    ax.set_xlabel('Potential vs. Li/Li$^+$ (V)', fontsize=20)
    ax.set_ylabel('Temperature (K)', fontsize=20)
    if discharge:
        ax.invert_xaxis()
    else:
        pass
    colorbar = fig.colorbar(PD, ticks=ticky, pad=0.05)
    colorbar.ax.set_yticklabels(labels)
    colorbar.ax.tick_params(size=0)
    colorbar.ax.minorticks_off()
    plt.tight_layout()
    if save:
        if discharge:
            plt.savefig('discharge-surface-pd.pdf', dpi=500)
        else:
            plt.savefig('charge-surface-pd.pdf', dpi=500)
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__ + "\n{} {}".format(__date__, __author__),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "surface_pd_data",
        help="Path to the surface pd data.",
        nargs="+",
        type=str
    )

    parser.add_argument(
        "--lithium-like-species", "-L",
        help="Define the lithium like species",
        type=str,
        default="Li"
    )

    parser.add_argument(
        "--oxygen-like-species", "-O",
        help="Define the oxygen like species",
        type=str,
        default="O"
    )

    parser.add_argument(
        "--low-T", "-lt",
        help="Low temperature boundary",
        type=float,
        default=1
    )

    parser.add_argument(
        "--high-T", "-ht",
        help="High temperature boundary",
        type=float,
        default=1500
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

    parser.add_argument(
        '--discharge', '-d',
        help='If the surface pd represents charge or discharge',
        action='store_true'
    )

    parser.add_argument(
        '--save', '-s',
        help='Whether to save the charge/discharge surface pd.',
        action='store_true'
    )

    args = parser.parse_args()

    surface_pd_plot(args.surface_pd_data,
                    args.lithium_like_species,
                    args.oxygen_like_species,
                    args.low_T,
                    args.high_T,
                    args.color_by_Li_content,
                    args.color_by_O_content,
                    args.functional,
                    args.discharge,
                    args.save
                    )
