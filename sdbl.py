# Graph Structure for SDBL
# Author: Henry Amrhein
# Date: 3 JAN 2020

import os.path
import matplotlib.patches
import io
import imageio

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import graph


ACTION_COLOR_DICT = {"activation": "#008000a0", "binding": "#0000a080",
                     "catalysis": "#800080a0", "expression": "#ff8c00a0",
                     "inhibition": "#ff0000a0", "ptmod": "#ff00ff80",
                     "reaction": "#000000a0"}

EVIDENCE_COLOR_DICT = {"neighborhood": "#00a000a0", "fusion": "#ff0000c0",
                       "cooccurence": "#0000a080", "coexpression": "#000000a0",
                       "experimental": "#800080a0", "database": "#00808080",
                       "textmining": "#00d000a0"}


class SdblException(Exception):
    pass


class Sdbl:
    def __init__(self, dbfile):
        self.G = graph.SdblGraph(dbfile)

    def __len__(self):
        return len(self.G)

    def __contains__(self, name):
        return name in self.G.gobj

    def __str__(self):
        return str(self.G.gobj)

    def nodes(self):
        if len(self.G) > 0:
            return self.G.gobj.nodes()

    def edges(self):
        if len(self.G) > 0:
            return self.G.gobj.edges()

    def build_action_graph(self, gl, cutoff, modes, **kwargs):
        self.G.color_dict = ACTION_COLOR_DICT
        self.G.build_graph(gl, cutoff, modes, schema="action", **kwargs)

    def build_evidence_graph(self, gl, cutoff, modes, **kwargs):
        self.G.color_dict = EVIDENCE_COLOR_DICT
        self.G.build_graph(gl, cutoff, modes, schema="evidence", **kwargs)

    def reset(self):
        dbfile = self.G.dbfile
        self.G = graph.SdblGraph(dbfile)

    def draw(self, filename, format=None):
        self.G.draw(filename, format)

    def layout(self, prog="sfdp"):
        self.G.layout(prog)

    def write(self, filename):
        self.G.write(filename)

    def colorize_by_column(self, frame, column, warm_cm=plt.cm.Oranges,
            cool_cm=plt.cm.Blues, center=0.0, cm_lower_stop=0.333,
            cm_upper_stop=0.8, n_quant=3):
        cool = cool_cm(np.linspace(cm_upper_stop, cm_lower_stop, n_quant))
        warm = warm_cm(np.linspace(cm_lower_stop, cm_upper_stop, n_quant))
        cmap = colors.ListedColormap(np.vstack((cool, warm)))

        data = frame[column]

        qspace = np.linspace(0, 1, n_quant + 1)
        dnq = data[data < center].quantile(qspace[:n_quant])
        upq = data[data > center].quantile(qspace[1:])

        bins = np.concatenate([dnq, [center], upq])

        cmap_idx = np.digitize(data, bins) - 1

        for n, c in zip(data.index, cmap_idx):
            self.G.set_node_fill_color(n, colors.to_hex(cmap(c)))

    def to_matplotlib_figure(self, ax=None):
        """Returns a Figure object or draws directly to Axes"""

        if not self.gobj.has_layout:
            errstr = "Attempt to draw SdblGraph with no layout"
            raise SdblGraphException(errstr)

        bio = io.BytesIO()
        self.G.draw(bio, format="png")
        bio.seek(0)

        img = imageio.imread(bio)

        if ax is None:
            figsize = (img.shape[0] / 100.0, img.shape[1] / 100.0)
            fig = plt.figure(figsize=figsize, dpi=100)
            ax = fig.add_subplot(1, 1, 1)

        ax.imshow(img, interpolation="nearest")
        ax.set_axis_off()
        bio.close()

        if ax is None:
            fig.tight_layout()
            return fig

    def add_disconnected_right(self):
        if not self.G.gobj.has_layout:
            return

        bb = [float(c) for c in self.G.gobj.graph_attr["bb"].split(",")]

        x = bb[2] * 1.25
        y = bb[3] - 40

        for n in sorted(self.G.disconnected):
            self.G.gobj.add_node(n, pos="{},{}".format(x, y))
            y -= 80

            if y < 40:
                x += 120
                y = bb[3] - 40

    def to_adjacency_matrix(self):
        adict = dict()

        for e in sorted(self.edges()):
            if "weight" in e.attr:
                w = e.attr["weight"]
            else:
                w = 1.0

            if e[0] not in adict:
                adict[e[0]] = dict()

            adict[e[0]][e[1]] = w

            if e[1] not in adict:
                adict[e[1]] = dict()

            adict[e[1]][e[0]] = w

        adjmat = pd.DataFrame.from_dict(adict, dtype=float)
        adjmat.sort_index(inplace=True)
        adjmat.sort_index(axis=1, inplace=True)

        return adjmat

    def add_legend(self, imgfile):
        if not self.G.gobj.has_layout:
            return

        img = imageio.imread(imgfile)
        imsize = img.shape
        xpos = -20 - pixels_to_points(imsize[1], 100) / 2
        bb = [float(c) for c in self.graph_attr["bb"].split(",")]

        if imsize[0] > bb[3]:
            ypos = bb[3] / 2
        else:
            ypos = bb[3] - (((imsize[0] / 100) * 72) / 2)

        spos = "{},{}".format(xpos, ypos)

        self.G.gobj.add_node("legend", shape="box", image=image,
                pos=spos, label=" ")

    def bloom_centroids(self, b):
        ew = dict()

        for g in b.amat.index:
            z = b.amat.loc[g]

            for c in (z[z > 0]).index:
                ekey = tuple(sorted([g, c]))
                ew[ekey] = b.amat.loc[g, c]

                if g in self:
                    self.G.gobj.add_edge(g, c, weight=ew[ekey], style="invis")

        cstable = b.amat.count()
        mx = cstable.max()
        mn = cstable.min()
        cs = (cstable - mn) / (mx - mn) + 2

        for c in b.amat:
            s = cs.loc[c] * 1.25
            t = b.gosig.loc[int(c), "Term"]
            n = self.G.gobj.get_node(c)
            n.attr["label"] = t.replace(" ", "\n")
            n.attr["shape"] = "rect"
            n.attr["style"] = "filled"
            n.attr["fillcolor"] = "#ffffc2"
            n.attr["fontsize"] = 20
            n.attr["width"] = s
            n.attr["height"] = s

