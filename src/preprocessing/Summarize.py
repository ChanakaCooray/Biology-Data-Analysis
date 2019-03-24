import os
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import sys
import csv
import pandas as pd
import operator

import networkx as nx
from networkx.algorithms import bipartite


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
        if group_name_new == group_name_prev:
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


def convert_unweighted(df, prefix_size):
    header_list = list(df)

    df_unweighted_matrix = df[['OTU ID']].copy()

    group_name_prev = ""
    for header in header_list[1:]:
        group_name_new = header[:prefix_size]
        if group_name_new == group_name_prev:
            continue
        else:
            group_name_prev = group_name_new

        group_columns = df.filter(regex=group_name_new).copy()
        group_columns["{}_mean".format(group_name_new)] = group_columns.mean(axis=1)

        df_unweighted_matrix[group_name_new] = group_columns["{}_mean".format(group_name_new)].apply(
            lambda x: 1 if x > 0 else 0)

    return df_unweighted_matrix


def analyze_matrix(df):
    header_list = list(df)

    df_output = pd.DataFrame(columns=['Name', 'Dot product', 'probability 1', 'probability 2'])

    # dot_product_map = {}
    i = 0
    for header in header_list[1:]:
        if "CA" in header:
            plant = header[3:5]
            location = header[1:3]
        else:
            continue

        for header2 in header_list[1:]:
            if "CA" in header2 or location not in header2:
                continue

            count1 = df[header].sum()
            count2 = df[header2].sum()

            dot_product = np.dot(df[header], df[header2])
            plant2 = header2[3:5]

            # dot_product_map["{}_{}_{}_probability1".format(plant, plant2, location)] = dot_product / count1
            # dot_product_map["{}_{}_{}_probability2".format(plant, plant2, location)] = dot_product / count2

            df_output = df_output.append(
                {'Name': "{}_{}_{}".format(plant, plant2, location), 'Dot product': dot_product,
                 'probability 1': dot_product / count1, 'probability 2': dot_product / count2}, ignore_index=True)

            i += 1

    return df_output


def draw_graph(edge_list):

    B = nx.Graph()

    for edge in edge_list:
        B.add_edge(edge[0], edge[1])

    X, Y = bipartite.sets(B)

    print(len(X))
    print(len(Y))

    pos = dict()
    pos.update((n, (1, i)) for i, n in enumerate(X))  # put nodes from X at x=1
    pos.update((n, (2, i)) for i, n in enumerate(Y))  # put nodes from Y at x=2
    nx.draw(B, pos=pos)
    plt.show()


def convert_edgelist(df, output_dir):
    header_list = list(df)

    edge_list = []
    for header in header_list[1:]:
        df_temp = df.loc[df[header] == 1]['OTU ID']

        for entry in df_temp.values:
            edge_list.append((header, entry))

    output_file = os.path.join(output_dir, "edge_list.txt")
    out = open(output_file, "w")
    for entry in edge_list:
        out.write("{} {} {}\n".format(entry[0], entry[1], 1))
    out.close()
    
    draw_graph(edge_list)
    
    return edge_list


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

    # df_BG_summarize = summarize(df_BG, 6)
    # df_rhizo_summarize = summarize(df_rhizo, 5)
    # df_roots_summarize = summarize(df_roots, 4)

    output_dir = "output"

    # create output directory if not exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # df_BG_summarize.to_csv(os.path.join(output_dir, "BG_summarize.csv"), index=None, sep=',', mode='w')
    # df_rhizo_summarize.to_csv(os.path.join(output_dir, "rhizo_summarize.csv"), index=None, sep=',', mode='w')
    # df_roots_summarize.to_csv(os.path.join(output_dir, "roots_summarize.csv"), index=None, sep=',', mode='w')

    df_unweighted_matrix = convert_unweighted(df_rhizo, 5)
    # df_unweighted_matrix.to_csv(os.path.join(output_dir, "rhizo_unweighted_matrix.csv"), index=None, sep=',', mode='w')

    convert_edgelist(df_unweighted_matrix, output_dir)
    
    df_output = analyze_matrix(df_unweighted_matrix)

    df_output.to_csv(os.path.join(output_dir, "analyzed_output.csv"), index=None, sep=',', mode='w')


if __name__ == '__main__':
    main()
