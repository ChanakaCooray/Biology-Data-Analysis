import os
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import sys
import csv
import pandas as pd


def verify_headers_BG(df):
    header_list = list(df)

    for header in header_list[1:]:
        number = header[6:]
        if not number.isdigit():
            print("BG: header: {} is not in the proper format.".format(header))


def verify_headers_rhizo(df):
    header_list = list(df)

    for header in header_list[1:]:
        number = header[5:]
        if not number.isdigit():
            print("rhizo: header: {} is not in the proper format.".format(header))


def verify_headers_roots(df):
    header_list = list(df)

    for header in header_list[1:]:
        number = header[4:]
        if not number.isdigit():
            print("roots: header: {} is not in the proper format.".format(header))


def summarize(df, prefix_size):
    header_list = list(df)

    df_summarize = df[['OTU ID']].copy()

    group_name_prev = ""
    for header in header_list[1:]:
        group_name_new = header[:prefix_size]
        if group_name_new==group_name_prev:
            continue
        else:
            group_name_prev = group_name_new

        # spike_cols = [col for col in df.columns if group_name_new in col]

        group_columns = df.filter(regex=group_name_new)
        df_summarize["{}_sum".format(group_name_new)] = group_columns.sum(axis=1)
        df_summarize["{}_mean".format(group_name_new)] = group_columns.mean(axis=1)
        df_summarize["{}_non_zero".format(group_name_new)] = group_columns.astype(bool).sum(axis=1)

        # print(list(group_columns))

    return df_summarize


def main():
    data_dir = "data"
    input_file1 = os.path.join(data_dir, "BG.csv")
    input_file2 = os.path.join(data_dir, "rhizo.csv")
    input_file3 = os.path.join(data_dir, "roots.csv")

    df_BG = pd.read_csv(input_file1, sep=",", header=0)
    df_rhizo = pd.read_csv(input_file2, sep=",", header=0)
    df_roots = pd.read_csv(input_file3, sep=",", header=0)

    verify_headers_BG(df_BG)
    verify_headers_rhizo(df_rhizo)
    verify_headers_roots(df_roots)

    df_BG_summarize = summarize(df_BG, 6)
    df_rhizo_summarize = summarize(df_rhizo, 5)
    df_roots_summarize = summarize(df_roots, 4)

    output_dir = "output"

    # create output directory if not exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    df_BG_summarize.to_csv(os.path.join(output_dir,"BG_summarize.csv"), index=None, sep=',', mode='w')
    df_rhizo_summarize.to_csv(os.path.join(output_dir,"rhizo_summarize.csv"), index=None, sep=',', mode='w')
    df_roots_summarize.to_csv(os.path.join(output_dir,"roots_summarize.csv"), index=None, sep=',', mode='w')

if __name__ == '__main__':
    main()
