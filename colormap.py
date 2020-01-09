# Diverging and Linear Colormaps for use by SDBL
# Author: Henry Amrhein
# Date: 16 OCT 2019

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np

"""Diverging and linear colormaps for use by SDBL"""

__version__ = 1.0


class SdblColormapException(Exception):
    pass


class SdblColormapComponent:
    def __init__(self, cmap, parent=None):
        self.cmap = cmap
        self.parent = parent
        self.nbins = 3
        self.lstop = 0.3
        self.ustop = 0.8
        self.color_list = self.cmap(np.linspace(self.lstop, self.ustop, self.nbins))
        self.lum_list = np.sum(self.color_list * (0.299, 0.587, 0.114, 0.0), axis=1)

    def _refresh(self):
        self.color_list = self.cmap(np.linspace(self.lstop, self.ustop, self.nbins))
        self.lum_list = np.sum(self.color_list * (0.299, 0.587, 0.114, 0.0), axis=1)

        if self.parent is not None:
            self.parent._refresh()

    def __len__(self):
        return self.nbins

    @property
    def lower_stop(self):
        return self.lstop

    @lower_stop.setter
    def lower_stop(self, value):
        self.lstop = value
        self._refresh()

    @property
    def upper_stop(self):
        return self.ustop

    @upper_stop.setter
    def upper_stop(self, value):
        self.ustop = value
        self._refresh()

    @property
    def n_quantiles(self):
        return self.nbins

    @n_quantiles.setter
    def n_quantiles(self, value):
        self.nbins = value
        self._refresh()

    def reverse(self):
        lower = self.lstop
        self.lstop = self.ustop
        self.ustop = lower
        self._refresh()


class SdblDivergingColormap:
    def __init__(self, colormap0, colormap1, parent=None):
        self.component0 = SdblColormapComponent(colormap0, self)
        self.component1 = SdblColormapComponent(colormap1, self)
        self.parent = parent
        self.component0.reverse()
        self.cmap = colors.ListedColormap(np.vstack((
            self.component0.color_list,
            self.component1.color_list
            ))
        )

    def _refresh(self):
        self.cmap = colors.ListedColormap(np.vstack((
            self.component0.color_list,
            self.component1.color_list
            ))
        )

        if self.parent is not None:
            self.parent._refresh()

    @property
    def N(self):
        return self.cmap.N

    @property
    def color_list(self):
        return self.cmap.colors

    @property
    def luminance_list(self):
        return np.sum(self.cmap.colors * (0.299, 0.587, 0.114, 0.0), axis=1)

    def __call__(self, value):
        return self.cmap(value)


class SdblLinearColormap:
    def __init__(self, colormap0, parent=None):
        self.component = SdblColormapComponent(colormap0, self)
        self.parent = parent
        self.cmap = colors.ListedColormap(self.component.color_list)

    def _refresh(self):
        self.cmap = colors.ListedColormap(self.component.color_list)
        
        if self.parent is not None:
            self.parent._refresh()

    @property
    def N(self):
        return self.cmap.N

    @property
    def color_list(self):
        return self.cmap.colors

    @property
    def luminance_list(self):
        return np.sum(self.cmap.colors * (0.299, 0.587, 0.114, 0.0), axis=1)

    def __call__(self, value):
        return self.cmap(value)


def hex_colorlist(cmap):
    return [colors.to_hex(c) for c in cmap.color_list]


def sdbl_quantile_diverging_colorize_by_numeric_series(O, series,
        cool_cm=plt.cm.Blues, warm_cm=plt.cm.Oranges, center=0.0,
        cm_lower_stop=0.3, cm_upper_stop=0.8, n_quantiles=3):

    cmap = SdblDivergingColormap(cool_cm, warm_cm)

    cmap.component0.upper_stop = cm_lower_stop
    cmap.component0.lower_stop = cm_upper_stop
    cmap.component0.n_quantiles = n_quantiles

    cmap.component1.lower_stop = cm_lower_stop
    cmap.component1.upper_stop = cm_upper_stop
    cmap.component1.n_quantiles = n_quantiles

    qspace = np.linspace(0, 1, n_quantiles + 1)
    dnq = series[series < center].quantile(qspace[:n_quantiles])
    upq = series[series > center].quantile(qspace[1:])

    bins = np.concatenate([dnq, [center], upq])

    color_idx = np.digitize(series, bins) - 1

    for n, c in zip(series.index, color_idx):
        if n in O.A:
            attr = O.A.get_node(n).attr
            attr["style"] = "filled"
            attr["fillcolor"] = colors.to_hex(cmap(c))


def sdbl_quantile_linear_colorize_by_numeric_series(O, series,
        cm=plt.cm.Blues, cm_lower_stop=0.3, cm_upper_stop=0.8,
        n_quantiles=3):

    cmap = SdblLinearColormap(cm)

    cmap.component.lower_stop = cm_lower_stop
    cmap.component.upper_stop = cm_upper_stop
    cmap.component.n_quantiles = n_quantiles

    qspace = np.linspace(0, 1, n_quantiles)
    color_idx = np.digitize(series, series.quantile(qspace)) - 1

    for n, c in zip(series.index, color_idx):
        if n in O.A:
            attr = O.A.get_node(n).attr
            attr["style"] = "filled"
            attr["fillcolor"] = colors.to_hex(cmap(c))
