#!/usr/bin/python3

import argparse

import sdbl
import colormap

def build_linear_cmap(hot, cold, ncolors=256, nbins=16):
    lcm = colors.LinearSegmentedColormap.from_list(name="lcm",
            colors=[cold, hot], N=ncolors)
    cm = colormap.SdblLinearColormap(lcm)
    cm.component.lower_stop = 0.0
    cm.component.upper_stop = 1.0
    cm.component.n_quantiles = nbins
    return cm

def colorize_graph(gobj, cm, data, bins, fontcolor1="black",
        fontcolor2="white"):
    color_idx = np.digitize(data, bins) - 1
    G = gobj.G

    for n, c in zip(data.index, color_idx):
        if n in G:
            nodecolor = colors.to_hex(cm(c))
            if cm.luminance_list[c] < 0.55:
                textcolor = fontcolor2
            else:
                textcolor = fontcolor1

            attr = G.gobj.get_node(n).attr
            attr["style"] = "filled"
            attr["fillcolor"] = nodecolor
            attr["fontcolor"] = fc

def build_graph(Sobj, gl):
    emodes = ["database", "experimental"]
    cutoff = 200
    gattr = {"splines": "true",
             "mode": "maxent",
             "K": "0.15",
             "repulsiveforce": "5.0"}
    eattr = {"len": "0.15"}
    nattr = {"fontname": "Arial",
             "fontsize": "26",
             "height": "0.80"}
    Sobj.build_evidence_graph(gl, cutoff=cutoff, modes=emodes,
            graphattr=gattr, edgeattr=eattr, nodeattr=nattr)
    Sobj.layout(prog="sfdp")
    Sobj.add_disconnected_right()

def build_colorbar(cm, data):
    norm = colors.Normalize(data.min(), data.max())
    fig = plt.figure(figsize=(1.28, 5.12), dpi=100)
    ax = fig.add_subplot(1, 1, 1)
    cb = matplotlib.colorbar.ColorbarBase(ax, cm.component.cmap, norm=norm,
            orientation="vertical", format="%.2f")
    for t in cb.ax.get_yticklabels():
        t.set_fontname("Arial")
    fig.tight_layout()
    return fig

def main(args):
    S = sdbl.Sdbl(args.sdblfile)

    meta_df = pd.read_hdf(args.datafile, "/metadata")
    marker_counts = pd.read_hdf(args.datafile, "/counts_by_cluster_normalized")

    fn_template = "{}/{}".format(args.output_dir, args.img_tmpl)
    cb_template = "{}/{}".format(args.output_dir, args.cb_tmpl)

    cmap = dict(zip(meta_df["unified_label"], meta_df["color"]))
    gmap = dict(zip(meta_df["unified_label"], meta_df["gene_series"]))

    qspace = np.linspace(0, 1, 16)

    for lbl in marker_counts:
        cm = build_linear_cmap(cmap[lbl])
        gl = pd.read_hdf(args.datafile, "/marker_genes/{}".format(gmap[lbl])).tolist()
        if len(gl) == 0:
            continue

        build_graph(S, gl)

        data = marker_counts.reindex(gl)[lbl].dropna()

        colorize_graph(S, cm, data, data.quantile(qspace))

        fig = build_colorbar(cm, data)

        clustername = gmap[lbl]

        for ext in ("png", "pdf", "svg"):
            S.draw(fn_template.format(label=clustername, ext=ext))
            fig.savefig(cb_template.format(label=clustername, ext=ext))

        S.write(fn_template.format(label=clustername, ext="dot"))
        S.reset()

        plt.close(fig)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("datafile", help="HDF5 file containing counts, markers, and metadata tables")
    parser.add_argument("sdblfile", help="sdbl database file")
    parser.add_argument("img_tmpl", help="filename template for network graphs - format: {cluster_name}_some_text.{ext}")
    parser.add_argument("cb_tmpl", help="filename template for colorbars - format: {cluster_name}_some_text.{ext}")
    parser.add_argument("output_dir", default=".", help="directory to put output into.  default: current directory")
    args = parser.parse_args()
    main(args)
