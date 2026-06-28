"""Matplotlib renderers shared by the pages (figures + saved PNGs).

Each function returns a Matplotlib ``Figure`` for ``st.pyplot`` and can be saved
to PNG for the run dir. Categorical pyroconvection panels reuse the colour/label
maps exported by :mod:`pyflam.atmosphere`.
"""

from __future__ import annotations

import numpy as np


def _fig(figsize=(7, 4)):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt.subplots(figsize=figsize)


def meteo_timeseries(report, variables=None):
    """Line charts of selected meteo variables over the run window.

    ``report`` is a :class:`pyflam.meteo_report.MeteoReport`. Returns a Figure with
    one subplot per variable.
    """
    from pyflam.meteo_report import _UNITS
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    variables = variables or ["temperature", "relative_humidity", "wind_speed",
                              "m_1h", "cape", "boundary_layer_height", "plume_factor"]
    variables = [v for v in variables if v in report.series
                 and any(x is not None for x in report.series[v])]
    n = len(variables) or 1
    fig, axes = plt.subplots(n, 1, figsize=(8, 1.7 * n), sharex=True, squeeze=False)
    t = report.times
    for ax, v in zip(axes[:, 0], variables):
        y = [np.nan if x is None else x for x in report.series[v]]
        ax.plot(t, y, marker="o", ms=3, lw=1.2)
        ax.set_ylabel(f"{v}\n[{_UNITS.get(v, '')}]", fontsize=8)
        ax.grid(alpha=0.3)
    axes[-1, 0].set_xlabel("time")
    fig.autofmt_xdate()
    fig.tight_layout()
    return fig


def pyroconv_panels(cats, lonlat_bbox, hours, *, title="", flip_for_display=False):
    """Categorical pyroconvection-type panels (one per forecast hour).

    ``cats`` is ``(nhours, nlat, nlon)`` of class levels 0..4; ``lonlat_bbox`` is
    ``(lon_min, lon_max, lat_min, lat_max)``. Colours/labels come from pyflam.
    """
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from pyflam.atmosphere import (
        PYROCONVECTION_TYPES, PYROCONVECTION_TYPE_LEVEL,
        PYROCONVECTION_TYPE_COLOR, PYROCONVECTION_TYPE_LABEL)

    cmap = ListedColormap([PYROCONVECTION_TYPE_COLOR[t] for t in PYROCONVECTION_TYPES])
    norm = BoundaryNorm(np.arange(-0.5, 5.5, 1), cmap.N)
    ext = list(lonlat_bbox)
    nh = len(hours)
    fig, ax = plt.subplots(1, nh, figsize=(2.1 * nh, 3.0), constrained_layout=True,
                           squeeze=False)
    ax = ax[0]
    for hi, hour in enumerate(hours):
        a = cats[hi][::-1] if flip_for_display else cats[hi]
        ax[hi].imshow(a, origin="lower", extent=ext, cmap=cmap, norm=norm,
                      aspect="auto", interpolation="nearest")
        ax[hi].set_title(f"{hour:02d}Z", fontsize=8)
        ax[hi].set_xticks([]); ax[hi].set_yticks([])
    if title:
        fig.suptitle(title, fontsize=11)
    leg = [Patch(facecolor=PYROCONVECTION_TYPE_COLOR[t], edgecolor="0.4",
                 label=f"{PYROCONVECTION_TYPE_LEVEL[t]}  {PYROCONVECTION_TYPE_LABEL[t]}")
           for t in PYROCONVECTION_TYPES]
    fig.legend(handles=leg, loc="lower center", ncol=5, fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, -0.12))
    return fig


def raster(array, *, title="", cmap="magma", lonlat_bbox=None, label=""):
    """A single continuous raster (burn prob, flame length, fireline intensity)."""
    fig, ax = _fig((6, 5))
    ext = list(lonlat_bbox) if lonlat_bbox is not None else None
    a = np.asarray(array, dtype=float)
    a = np.where(np.isfinite(a), a, np.nan)        # inf (unburned) -> blank
    im = ax.imshow(a, origin="upper", extent=ext, cmap=cmap,
                   interpolation="nearest")
    fig.colorbar(im, ax=ax, label=label, shrink=0.8)
    ax.set_title(title)
    fig.tight_layout()
    return fig


def flp_bars(result):
    """Flame-length probability distribution (mean over burned cells) as bars.

    ``result`` is a :class:`pyflam.mtt.BurnProbabilityResult`.
    """
    fig, ax = _fig((6, 3.5))
    centers = result.flame_length_class_centers()
    burned = result.burn_prob > 0
    if burned.any():
        mean_flp = np.nanmean(result.flp[:, burned], axis=1)
    else:
        mean_flp = np.zeros(len(centers))
    ax.bar([f"{c:.1f}" for c in centers], mean_flp, color="tab:orange")
    ax.set_xlabel("flame length class (ft)")
    ax.set_ylabel("P(class | burned)")
    ax.set_title("Flame-length probability (FLP)")
    fig.tight_layout()
    return fig


def fire_size_hist(fire_sizes, *, to_ha=1.0):
    """Histogram of per-fire burned area. ``fire_sizes`` in area units; ``to_ha``
    multiplies into hectares (e.g. 1e-4 when sizes are m²)."""
    fig, ax = _fig((6, 3.5))
    sizes_ha = np.asarray(fire_sizes, dtype=float) * to_ha
    sizes_ha = sizes_ha[np.isfinite(sizes_ha) & (sizes_ha > 0)]
    if sizes_ha.size:
        ax.hist(sizes_ha, bins=min(30, max(5, sizes_ha.size // 5)), color="tab:red")
    ax.set_xlabel("fire size (ha)")
    ax.set_ylabel("number of fires")
    ax.set_title("Fire-size distribution")
    fig.tight_layout()
    return fig


def operative_quiver(report):
    """Per-sector driving-force arrows (reuses OperativeReport.quiver)."""
    fig, ax = _fig((7, 4))
    report.quiver(ax=ax)
    fig.tight_layout()
    return fig


def save_png(fig, path: str, dpi: int = 140) -> str:
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    return path
