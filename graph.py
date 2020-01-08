# Graph Structure for SDBL
# Author: Henry Amrhein
# Date: 19 OCT 2019

import io

import pygraphviz
import imageio

import edge_engine
import sql

import matplotlib.pyplot as plt


class SdblGraphException(Exception):
    pass


class SdblGraph:
    def __init__(self, dbfile, name=None):
        self.gattr = {
                "overlap": "false",
                "splines": "false",
                "bgcolor": "white",
                "outputorder": "edgesfirst",
                "sep": "+20,20",
                "esep": "+12,12",
                "start": "12345"
                }

        self.nattr = {
                "style": "filled",
                "fillcolor": "white",
                "fontname": "sans serif bold"
                }

        self.dbfile = dbfile
        self.gobj = pygraphviz.AGraph(name=name, directed=True, strict=False)
        self.looping = list()
        self.disconnected = None
        self.color_dict = None

    def __del__(self):
        self.gobj.close()

    def __str__(self):
        return str(self.gobj)

    def __contains__(self, node):
        if self.gobj is None:
            return False

        return node in self.gobj

    def __len__(self):
        if self.gobj is None:
            return 0
        else:
            return len(self.gobj)

    def initialize_graph(self, name, label, labelpos, labelfont,
            labelfontsize, graphattr, edgeattr, nodeattr):
        """Initialize the Graphviz AGraph structure"""

        self.gobj.graph_attr.update(self.gattr)
        self.gobj.node_attr.update(self.nattr)

        if label is not None:
            self.gobj.graph_attr["label"] = label
            self.gobj.graph_attr["labelloc"] = labelpos
            self.gobj.graph_attr["fontname"] = labelfont
            self.gobj.graph_attr["fontsize"] = labelfontsize

        if graphattr is not None:
            self.gobj.graph_attr.update(graphattr)

        if edgeattr is not None:
            self.gobj.edge_attr.update(edgeattr)

        if nodeattr is not None:
            self.gobj.node_attr.update(nodeattr)
        
    def build_graph(self, gene_list, cutoff, modes, schema="action",
            connected=True, looping=False, penwidth_multiplier=2, name="",
            graphattr=None, edgeattr=None, nodeattr=None, label=None,
            labelpos='t', labelfont="sans serif bold", labelfontsize=48):
        """Build the graph from a list of gene names, a cutoff score, and
        a set of edge modes."""

        self.initialize_graph(name, label, labelpos, labelfont,
                labelfontsize, graphattr, edgeattr, nodeattr)

        dbh = sql.SdblSql(self.dbfile)

        if schema == "action":
            res = dbh.actions_query_multiple_genes(gene_list,
                    cutoff_score=cutoff)
        else:
            res = dbh.evidence_query_multiple_genes(gene_list,
                    cutoff_score=cutoff)

        dbh.close()

        ee = edge_engine.SdblEdgeEngine(res)
        ee.generate_edges()

        for k, v in ee:
            if k[0] == k[1]:
                self.looping.append(k)
                if looping == False:
                    continue

            if k[2] not in modes:
                continue

            self.gobj.add_edge(k[0], k[1], weight=v.score, key=k,
                    dir=v.direction, style="solid", arrowhead=v.arrowtype,
                    arrowtail=v.arrowtype, color=self.color_dict[k[2]],
                    penwidth=v.penwidth * penwidth_multiplier)

            self.disconnected = [g for g in gene_list if g not in self.gobj]

    def layout(self, prog="sfdp"):
        """Arrange the nodes using a specified layout program"""

        self.gobj.layout(prog=prog)

        A = pygraphviz.AGraph(str(self.gobj))

        self.gattr = dict(A.graph_attr)
        self.eattr = dict(A.edge_attr)
        self.nattr = dict(A.node_attr)

    def draw(self, filename, format=None):
        """write a graphic file of the current graph.
        Use 'dot -T:' to list the available output formats"""

        self.gobj.draw(filename, format=format)

    def write(self, filename):
        """write a DOT file of the current graph"""
        self.gobj.write(filename)

    def set_node_fill_color(self, node, color):
        if node not in self.gobj:
            estr = "Node does not exist"
            raise SdblGraphException(estr)
        attr = self.gobj.get_node(node).attr
        attr["style"] = "filled"
        attr["fillcolor"] = color

    def uncolorize(self):
        for n in self.gobj:
            self.set_node_fill_color(n, "white")

    def load(self, dotfile):
        self.gobj = pygraphviz.AGraph(filename=dotfile)
        
