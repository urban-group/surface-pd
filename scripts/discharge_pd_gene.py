#!/usr/bin/env python

"""
This code will be used to generate the discharge surface pd data.
"""

__author__ = "Xinhao Li"
__email__ = "xl2778@columbia.edu"
__date__ = "2022-11-07"

import argparse

import pandas as pd

from surface_pd.plot.pd_data import PdData


def find_max_min_composition(data: pd.DataFrame,
                             column_name: str):
    a_max = data[column_name].max()
    b_min = data[column_name].min()
    return a_max - b_min


def add_composition2df(data1: str,
                       data2: str = None,
                       lithium_like_species: str = None,
                       oxygen_like_species: str = None
                       ):
    df1 = pd.read_csv(data1, sep='\s+', index_col=0)
    df1_c = PdData(df1,
                   lithium_like_species=lithium_like_species,
                   oxygen_like_species=oxygen_like_species,
                   functional=None)
    df1_c.standardize_pd_data()
    df1 = df1_c.dataframe

    if data2 is None:
        df2_li = 0
        df2_o = 0
    else:
        df2 = pd.read_csv(data2, sep='\s+', index_col=0)
        df2_c = PdData(df2,
                       lithium_like_species=lithium_like_species,
                       oxygen_like_species=oxygen_like_species,
                       functional=None)
        df2_c.standardize_pd_data()
        df2 = df2_c.dataframe

        df2_li = find_max_min_composition(df2,
                                          column_name=lithium_like_species)
        df2_o = find_max_min_composition(df2,
                                         column_name=oxygen_like_species)

    df1_li = find_max_min_composition(df1,
                                      column_name=lithium_like_species)
    df1_o = find_max_min_composition(df1,
                                     column_name=oxygen_like_species)

    overall_li = df1_li + df2_li
    overall_o = df1_o + df2_o

    df1['li_composition'] = (df1['Li'] - df1['Li'].min() + df2_li) / overall_li
    df1['o_composition'] = (df1['O'] - df1['O'].min() + df2_o) / overall_o

    if data2 is None:
        return df1
    else:
        df2['li_composition'] = (df2['Li'] - df2['Li'].min()) / overall_li
        df2['o_composition'] = (df2['O'] - df2['O'].min()) / overall_o
        return df1, df2


def create_discharge_pd(data_files,
                        charge_pd_end_composition: float,
                        lithium_like_species: str = None,
                        oxygen_like_species: str = None,
                        save=False):
    if len(data_files) == 1:
        df1 = add_composition2df(data_files[0],
                                 lithium_like_species=lithium_like_species,
                                 oxygen_like_species=oxygen_like_species)
    else:
        df1, df2 = add_composition2df(
            data_files[0], data_files[1],
            lithium_like_species=lithium_like_species,
            oxygen_like_species=oxygen_like_species)

    for i in df1.index:
        if df1['o_composition'][i] > charge_pd_end_composition:
            df1['E'][i] = 0
    df1.drop(['li_composition', 'o_composition'], axis=1, inplace=True)
    if save:
        df1.to_csv('discharge-data1.csv', sep='\t')
    if len(data_files) == 1:
        pass
    else:
        for i in df2.index:
            if df2['o_composition'][i] > charge_pd_end_composition:
                df2['E'][i] = 0
        df2.drop(['li_composition', 'o_composition'], axis=1,
                 inplace=True)
        if save:
            df2.to_csv('discharge-data2.csv', sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__ + "\n{} {}".format(__date__, __author__),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        'surface_pd_data',
        help='Path to the charge surface pd data',
        nargs='+',
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
        '--charge_pd_end_composition', '-x',
        help='Oxygen composition at the end of charge.',
        type=float
    )

    parser.add_argument(
        '--save', '-s',
        help='Whether to save the discharge surfce pd data.',
        action='store_true'
    )

    args = parser.parse_args()

    create_discharge_pd(
        args.surface_pd_data,
        args.charge_pd_end_composition,
        args.lithium_like_species,
        args.oxygen_like_species,
        args.save
    )
