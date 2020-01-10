#!/usr/bin/python3

import argparse

import numpy as np
import pandas as pd
import pygraphviz


CLUSTER_COLORMAP = ['#ff8c00', '#5ca904', '#ffa500', '#1e90ff', '#0652ff', '#6b8ba4', '#006400', '#c875c4', '#ff0000', '#0000ff', '#014d4e', '#15b01a', '#be0119', '#cf6275', '#580f41', '#7f2b0a', '#a83c09', '#7f5e00', '#030aa7', '#800080', '#029386']

def gen_graph(df, penwidth=15, k="1.0", r="0.5"):
    A = pygraphviz.AGraph(strict=False)
    A.graph_attr["overlap"] = "prism1000"
    A.graph_attr["splines"] = "true"
    A.graph_attr["bgcolor"] = "#ffffff"
    A.graph_attr["outputorder"] = "edgesfirst"
    A.graph_attr["sep"] = "+20,20"
    A.graph_attr["esep"] = "+15,15"
    A.graph_attr["K"] = k
    A.graph_attr["start"] = "12345"
    A.graph_attr["repulsiveforce"] = r
    A.node_attr["style"] = "filled"
    A.node_attr["fillcolor"] = "#ffffc2"
    A.node_attr["fontname"] = "Arial"
    A.node_attr["fontsize"] = "16"
    ew = dict()

    for g in df.index:
        z = df.loc[g]
        for c in (z[z > 0]).index:
            ekey = tuple(sorted([g, c]))
            ew[ekey] = df.loc[g, c]
            A.add_edge(g, c, weight=df.loc[g, c], style="filled")

    v = ew.values()
    mn = min(v)
    mx = max(v)
    d = mx - mn

    if d == 0:
        d = 1

    for k, v in ew.items():
        e = A.get_edge(*k)
        pw = ((float(v) - mn) / d)
        e.attr["penwidth"] = pw * penwidth + 2.0

    return A

def scale_cluster_sizes(cstable):
    scale_lo = 1
    scale_hi = 2

    mn = cstable.Size.min()
    mx = cstable.Size.max()

    c = scale_hi - scale_lo
    cstable["nsize"] = c * (cstable.Size - mn) / (mx - mn) + 1.0

def filter_adjacency_table(df, cutoff_score=2.0):
    df = df.where(df > cutoff_score, 0.0)
    centers = df.T[df.any()].T.columns
    return centers

def set_cluster_node_size(A, cstable):
    for i in cstable.index:
        if i in A:
            n = A.get_node(i)
            c = cstable.loc[i, "nsize"]
            w = 1.25 * c
            h = 1.25 * c
            n.attr["width"] = str(w)
            n.attr["height"] = str(h)
            n.attr["fontsize"] = "48"

def set_cluster_node_edge_color(A, cstable):
    i = 1
    for j in sorted(cstable.index[:-1], key=int):
        if j in A:
            for e in A.edges(j):
                e.attr["color"] = CLUSTER_COLORMAP[i]
            n = A.get_node(j)
            n.attr["fillcolor"] = "{}90".format(CLUSTER_COLORMAP[i])
            i += 1

    for e in A.edges("Ubiquitous"):
        e.attr["color"] = CLUSTER_COLORMAP[0]
    n = A.get_node("Ubiquitous")
    n.attr["fillcolor"] = "{}90".format(CLUSTER_COLORMAP[0])

def annotate_shared_nodes(A, centers):
    leaves = [n for n in A.nodes() if n not in centers]
    degree = A.degree(leaves, with_labels=True)

    for k, v in degree.items():
        if v > 1:
            n = A.get_node(k)
            n.attr["style"] = "filled"
            n.attr["fillcolor"] = "#e0d0d0"
            n.attr["penwidth"] = "3"

def adjust_node_labels(A):
    n = A.get_node("Ubiquitous")
    n.attr["label"] = "Ubiqui-\ntous"
    n.attr["fontsize"] = 32

    n = A.get_node("ENSMUSG00000044690")
    n.attr["label"] = "ENSMUSG-\n00000044690"
    n.attr["width"] = 1.25

    n = A.get_node("9430076C15Rik")
    n.attr["label"] = "9430076-\nC15Rik"
    n.attr["width"] = 1.25

    n = A.get_node("IRC900814")
    n.attr["label"] = "IRC90-\n0814"
    n.attr["width"] = 1.25

def main(args):
    df = pd.read_table(args.adjacency_table, index_col=0)
    cstable = pd.read_table(args.cluster_size_table, index_col=0)

    scale_cluster_sizes(cstable)

    centers = filter_adjacency_table(df)

    A = gen_graph(df[centers], penwidth=10, k="0.75", r="1.0")

    set_cluster_node_size(A, cstable)
    set_cluster_node_edge_color(A, cstable)
    annotate_shared_nodes(A, centers)
    adjust_node_labels(A)

    A.layout(prog="sfdp")
    A.draw("{}/{}".format(args.output_dir, args.output_file))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate motif blossom plot")
    parser.add_argument("adjacency_table", help="Motif adjacency table")
    parser.add_argument("cluster_size_table", help="Cluster size table")
    parser.add_argument("output_file", help="Blossom plot output")
    parser.add_argument("--output_dir", default=".", help="Directory to put output into.  default: current directory")
    args = parser.parse_args()
    main(args)

