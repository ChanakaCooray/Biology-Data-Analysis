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


if __name__ == '__main__':
    main()
