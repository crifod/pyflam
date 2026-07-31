"""3-day Tuscany pyroconvection forecast (hybrid) from a single 00Z run, as one report.

Runs the tested daily pipeline (tests/pyroconv_daily.py) for three valid days -- the run day
and the next two forecast days -- then stitches the three days into combined 3-row figures and
a markdown/PDF report. The maps are ordered to answer three questions, in escalating order:

  1. **POTENTIAL** -- *what is the most intense pyroconvection this column could sustain?*
     The atmospheric upper bound: the ladder run with no fire-side condition at all.
  2. **FUEL-GATED** -- *is there enough fire here to make a plume at all?*  The same ladder
     with each cell forced to a surface plume below 10 MW/m of Byram intensity on the .lcp
     fuels (Tedim et al. 2018).
  3. **PFT MARGIN** -- *is there enough fire here to make a pyroCb in this specific column?*
     Per-cell total firepower over that column's own PyroCb Firepower Threshold (Tory &
     Kepert 2021). A far harder test than (2), and a continuous one.

A fourth, supporting map rides along:

  * **dry-pyrocloud decoupling** -- the continuous fireABL/ABL ratio for a reference intense
    fire, the DRY counterpart to the moist class ladder (diagnostic, no class label).

Questions 2 and 3 share an input -- both descend from the same Byram intensity field -- but
they are different tests and they disagree: (2) compares power *per metre of front* against a
fixed constant, (3) converts it to a *total* power through an assumed head-fire length and
compares it against a locally computed threshold. Cells can clear (2) and fail (3) by orders
of magnitude, which is the point of carrying both.

Every day is a **complete 24 h cycle, 3-hourly, always starting at 00Z of the run day**
(00, 03, 06, 09, 12, 15, 18, 21Z), on the Tuscany domain with ISTAT province borders.

Uses the day+1/day+2 forecast steps available in the same ICON-EU + ICON-2I 00Z run (ICON-EU to
+120 h, ICON-2I to +72 h). Output goes to its own folder.

Usage:  PYTHONPATH=src python scripts/forecast_3day_tuscany.py [YYYY-MM-DD run-day] [run]
Env:    PYROCONV_OUT (folder, default docs/forecast_<rundate>_3day), plus the usual
        PYROCONV_* knobs honoured by tests/pyroconv_daily.py. PYROCONV_HOURS is *not*
        honoured here -- the product is defined on the full 00-21Z cycle.
        SKIP_BASE_RUN=1 re-stitches the figures and report from rasters already on disk,
        without re-running the three daily pipelines -- the only way to revise a published
        report after DWD has dropped the run.
"""
import os
import subprocess
import sys
from datetime import datetime, timedelta, timezone

import numpy as np
import rasterio

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from pyflam.atmosphere import (PYROCONVECTION_TYPES, PYROCONVECTION_TYPE_LEVEL,
                               PYROCONVECTION_TYPE_LABEL, pyroconvection_colors)
from pyflam_gui.core.pyroconv import ABL_MIN_M, _REFERENCE_THETA_EXCESS_K
sys.path.insert(0, HERE)
import pyroconv_i18n as I18N        # report prose + figure labels, en/it

RUNDATE = sys.argv[1] if len(sys.argv) > 1 else datetime.now(timezone.utc).strftime("%Y-%m-%d")
RUN = int(sys.argv[2]) if len(sys.argv) > 2 else 0
RUNDT = datetime.strptime(RUNDATE, "%Y-%m-%d")
DAYS = [(RUNDT + timedelta(days=k)).strftime("%Y-%m-%d") for k in (0, 1, 2)]
# One complete 24 h cycle per valid day, 3-hourly, always starting at 00Z of the run day.
# Fixed on purpose: the three day-rows of every figure must be directly comparable, and
# day 0 must open at the run's own initial time, so PYROCONV_HOURS is not honoured here
# (run_daily overrides it for the daily pipeline too).
HOURS = [0, 3, 6, 9, 12, 15, 18, 21]
DAYTIME = (9, 12, 15, 18)                                 # peak-of-day window for the tables
OUT = os.environ.get("PYROCONV_OUT") or os.path.join(REPO, "docs", f"forecast_{RUNDATE}_3day")
LON0, LON1, LAT0, LAT1 = 9.6, 12.5, 42.2, 44.6            # Tuscany bbox (matches the daily runner)
CACHE = os.environ.get("PYROCONV_CACHE") or f"/tmp/pyflam_icon2i/{RUNDT:%Y%m%d}{RUN:02d}"
PROV_GEOJSON = (os.environ.get("PYROCONV_PROVINCES")
                or "/Users/cristianofoderi/DATI/boundaries/limits_IT_provinces.geojson")
SEA_COLOR = "#cfe4ef"
# "pyflam" (default) or "graf" -- see atmosphere.pyroconvection_colors. Passed through to the
# daily runner so the per-day panels and the stitched ones cannot disagree.
PALETTE = os.environ.get("PYROCONV_PALETTE", "pyflam")
COLORS = pyroconvection_colors(PALETTE)
NODATA_COLOR = "0.92"                                     # land with no classifiable column
# Colour range of the dry-pyrocloud decoupling ratio (fireABL / ABL). 1 = the reference fire
# grows no boundary layer of its own above the ambient one; the daily product uses the same
# 1..8 range, so the panels are readable side by side.
DECOUP_VMIN, DECOUP_VMAX = 1.0, 8.0
# Decoupling levels the summary table counts (ratio >= x over the classified land).
DECOUP_LEVELS = (2.0, 3.0)
# PFT margin = firepower / PyroCb Firepower Threshold. 1.0 is the criterion, and it is the
# colour pivot -- everything below is blue, everything at or above is red, with no ambiguity
# about which side of the threshold a cell is on. The scale is logarithmic because the field
# spans four decades: the median burnable cell sits near 5e-3, the domain max near 7.
PFT_MARGIN_VMIN, PFT_MARGIN_VMAX = 1.0e-3, 10.0
# Margin levels the summary table counts. 1.0 is the criterion; 0.1 is the within-one-decade
# band, which is where a bigger head fire than the assumed 700 m would start to matter.
PFT_MARGIN_LEVELS = (1.0, 0.1)
# Cells with no burnable fuel have exactly zero firepower, so their margin is exactly zero --
# a categorical "no", not a small number. Drawn white rather than at the bottom of the log
# scale, which would read as "very weak fire" when the truth is "no fire".
NOFUEL_COLOR = "#ffffff"


def provenance():
    """Code version + generation time, for the report footer.

    Products cannot be regenerated after DWD drops the run (~24 h), so a published raster set
    is frozen with whatever physics produced it. Without a version stamp a reader has no way to
    tell which corrections a given product predates -- and this pipeline has had several that
    changed magnitudes. See docs/graf_vs_pyflam_2026-07-26.md sec. 20.
    """
    import subprocess as _sp
    from datetime import datetime as _dt, timezone as _tz
    try:
        h = _sp.run(["git", "rev-parse", "--short", "HEAD"], cwd=REPO, capture_output=True,
                    text=True, timeout=10, check=True).stdout.strip()
        # -uno: untracked files are ignored. The run writes its own output folder into the
        # repo, so without this every report stamps itself "+local-changes" and the flag
        # loses all meaning. What matters for reproducibility is whether *tracked* code
        # differs from HEAD.
        dirty = _sp.run(["git", "status", "--porcelain", "-uno"], cwd=REPO, capture_output=True,
                        text=True, timeout=10).stdout.strip()
        ver = h + ("+local-changes" if dirty else "")
    except Exception:
        ver = "unknown"
    return ver, _dt.now(_tz.utc).strftime("%Y-%m-%d %H:%M UTC")


def sea_mask(lat, lon):
    """ICON-2I land fraction on the raster grid -> sea (frland < 0.5) as True, or None."""
    import warnings
    warnings.simplefilter("ignore")
    import xarray as xr
    path = os.path.join(CACHE, "FRLAND.grib")
    if not os.path.exists(path):
        return None
    ds = xr.open_dataset(path, engine="cfgrib", backend_kwargs={"indexpath": ""})
    v = list(ds.data_vars)[0]
    la = ds["latitude"].values
    ds = (ds.sel(latitude=slice(LAT1, LAT0), longitude=slice(LON0, LON1)) if la[0] > la[-1]
          else ds.sel(latitude=slice(LAT0, LAT1), longitude=slice(LON0, LON1)))
    fr = np.asarray(ds[v].values, float)
    if ds["latitude"].values[0] < ds["latitude"].values[-1]:
        fr = fr[::-1]                                    # north-up, to match the rasters
    return fr < 0.5


def provinces():
    """Tuscany province boundaries clipped to the bbox, or None if unavailable."""
    try:
        import geopandas as gpd
    except Exception:
        return None
    if not os.path.exists(PROV_GEOJSON):
        return None
    try:
        g = gpd.read_file(PROV_GEOJSON)
        if "reg_name" in g.columns:
            g = g[g["reg_name"] == "Toscana"]
        return g.clip((LON0, LAT0, LON1, LAT1)).boundary
    except Exception:
        return None


def run_daily(valid):
    """Invoke the daily pipeline for one valid day into OUT (hybrid default).

    ``PYROCONV_HOURS`` is forced to the full 00-21Z cycle: an inherited subset would
    leave a day-row short and break the raster loaders below.
    """
    env = dict(os.environ, PYROCONV_OUT=OUT, PYROCONV_VALID=valid,
               PYROCONV_HOURS=",".join(str(h) for h in HOURS), PYROCONV_PALETTE=PALETTE)
    print(f"\n=== {valid} (forecast day +{(datetime.strptime(valid,'%Y-%m-%d')-RUNDT).days}) ===")
    subprocess.run([sys.executable, os.path.join(REPO, "tests", "pyroconv_daily.py"),
                    RUNDATE, str(RUN), valid], env=env, check=True)


def load_stack(valid, kind):
    """Class stack (nhours, ny, nx) + (lat, lon) 1-D arrays, from the daily rasters."""
    rdir = os.path.join(OUT, f"rasters_hybrid_{valid}")
    arrs, lat, lon = [], None, None
    for h in HOURS:
        with rasterio.open(os.path.join(rdir, f"pyroconv_{kind}_{h:02d}Z.tif")) as d:
            arrs.append(d.read(1))
            if lat is None:
                t = d.transform
                lon = t.c + t.a * (np.arange(d.width) + 0.5)
                lat = t.f + t.e * (np.arange(d.height) + 0.5)      # north-up (t.e < 0)
    return np.stack(arrs), lat, lon


def load_diag_stack(valid, name):
    """One diagnostic stacked over the day's hours (nhours, ny, nx), north-up float32.

    The daily pipeline exports these next to the class rasters (``export_diagnostics``),
    so the 3-day product re-reads them instead of recomputing the columns.
    """
    rdir = os.path.join(OUT, f"rasters_hybrid_{valid}")
    arrs = []
    for h in HOURS:
        with rasterio.open(os.path.join(rdir, f"diag_{name}_{h:02d}Z.tif")) as d:
            arrs.append(d.read(1))
    return np.stack(arrs)


def valid_stack(valid):
    """Classifiable-column mask (nhours, ny, nx), rebuilt from the diagnostic rasters.

    The daily pipeline does not export its boolean ``valid`` field, so it is reconstructed
    here exactly as ``pyroconv.regrid_diagnostics`` defines it: a finite ABL, LCL/ABL ratio
    and mixed-layer gradient, with the ABL above the mixing floor. Needed because the
    decoupling ratio is a bare fireABL/ABL quotient -- it stays finite in columns the
    classifier rejects, and those must not be painted.
    """
    abl = load_diag_stack(valid, "abl")
    return (np.isfinite(abl) & np.isfinite(load_diag_stack(valid, "lcl_ratio"))
            & np.isfinite(load_diag_stack(valid, "ml_grad")) & (abl >= ABL_MIN_M))


def combined_figure(kind, path, lang):
    """WRF-style day x hour grid: geographic panels, sea masked, province borders, legend."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch

    # Class colours preceded by two sentinels: -2 = sea, -1 = PYROCONV_NODATA (a land column
    # the ladder could not run on). The sea sentinel moved from -1 to -2 when the classifier
    # started emitting -1 for unclassifiable cells -- sharing the value would have painted
    # every collapsed evening column as sea.
    colors = [SEA_COLOR, NODATA_COLOR] + [COLORS[t] for t in PYROCONVECTION_TYPES]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(-2.5, 5.5, 1), cmap.N)
    _, lat0, lon0 = load_stack(DAYS[0], kind)
    ext = [lon0.min(), lon0.max(), lat0.min(), lat0.max()]
    sea = sea_mask(lat0, lon0)
    prov = provinces()

    fig, ax = plt.subplots(len(DAYS), len(HOURS), figsize=(2.05 * len(HOURS), 2.5 * len(DAYS)),
                           constrained_layout=True, squeeze=False)
    for r, valid in enumerate(DAYS):
        stack, _, _ = load_stack(valid, kind)
        for c, hour in enumerate(HOURS):
            a = stack[c].astype(float)
            if sea is not None and sea.shape == a.shape:
                a = np.where(sea, -2, a)                 # paint sea cells with the sea colour
            axc = ax[r][c]
            axc.imshow(a, origin="upper", extent=ext, cmap=cmap, norm=norm,
                       aspect="auto", interpolation="nearest")
            if prov is not None:
                prov.plot(ax=axc, color="0.25", linewidth=0.4)
                axc.set_xlim(ext[0], ext[1]); axc.set_ylim(ext[2], ext[3])
            axc.set_xticks([]); axc.set_yticks([])
            if r == 0:
                axc.set_title(f"{hour:02d}Z", fontsize=9)
            if c == 0:
                axc.set_ylabel(f"{valid}\n(+{r}d)", fontsize=9)
    title = I18N.f(lang, "gated_title" if kind == "gated" else "potential_title")
    head = ("Piroconvezione in Toscana -- previsione 3 giorni" if lang == "it"
            else "Tuscany pyroconvection -- 3-day forecast")
    fig.suptitle(f"{head} -- {title}\n"
                 f"{I18N.f(lang, 'src_hybrid')} -- {I18N.f(lang, 'run')} {RUNDATE} {RUN:02d}Z",
                 fontsize=12)
    leg = [Patch(facecolor=COLORS[t], edgecolor="0.4",
                 label=f"{PYROCONVECTION_TYPE_LEVEL[t]}  {I18N.CLASS_LABEL[lang][t]}")
           for t in PYROCONVECTION_TYPES]
    leg.append(Patch(facecolor=NODATA_COLOR, edgecolor="0.4", label=I18N.f(lang, "nc")))
    fig.legend(handles=leg, loc="lower center", ncol=6, fontsize=9, frameon=False,
               title=I18N.f(lang, "class_legend_title"), bbox_to_anchor=(0.5, -0.09))
    fig.savefig(path, dpi=140, bbox_inches="tight"); plt.close(fig)


def decoupling_figure(path, lang):
    """Day x hour grid of the dry-pyrocloud decoupling ratio (fireABL / ABL).

    The DRY counterpart to the class maps, on the same geographic panels: how high a
    reference intense fire would grow its own boundary layer by sensible heat alone,
    divided by the ambient ABL. Continuous, diagnostic, no class label -- columns the
    classifier rejects stay blank rather than being painted with a bare quotient.
    """
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    _, lat0, lon0 = load_stack(DAYS[0], "potential")
    ext = [lon0.min(), lon0.max(), lat0.min(), lat0.max()]
    sea = sea_mask(lat0, lon0)
    prov = provinces()
    cmap = plt.get_cmap("magma").copy()
    cmap.set_bad(alpha=0.0)                      # unclassifiable land -> panel facecolor
    sea_cmap = ListedColormap([SEA_COLOR])

    fig, ax = plt.subplots(len(DAYS), len(HOURS), figsize=(2.05 * len(HOURS), 2.5 * len(DAYS)),
                           constrained_layout=True, squeeze=False)
    im = None
    for r, valid in enumerate(DAYS):
        dec, ok = load_diag_stack(valid, "decoupling"), valid_stack(valid)
        for c, hour in enumerate(HOURS):
            axc = ax[r][c]
            axc.set_facecolor(NODATA_COLOR)
            im = axc.imshow(np.ma.masked_invalid(np.where(ok[c], dec[c], np.nan)),
                            origin="upper", extent=ext, cmap=cmap,
                            vmin=DECOUP_VMIN, vmax=DECOUP_VMAX,
                            aspect="auto", interpolation="nearest")
            if sea is not None and sea.shape == dec[c].shape:
                axc.imshow(np.ma.masked_where(~sea, np.zeros_like(dec[c])), origin="upper",
                           extent=ext, cmap=sea_cmap, aspect="auto", interpolation="nearest")
            if prov is not None:
                prov.plot(ax=axc, color="0.25", linewidth=0.4)
                axc.set_xlim(ext[0], ext[1]); axc.set_ylim(ext[2], ext[3])
            axc.set_xticks([]); axc.set_yticks([])
            if r == 0:
                axc.set_title(f"{hour:02d}Z", fontsize=9)
            if c == 0:
                axc.set_ylabel(f"{valid}\n(+{r}d)", fontsize=9)
    fig.suptitle(f"{I18N.f(lang, 'decoup_title')} "
                 f"({I18N.f(lang, 'decoup_ref', k=_REFERENCE_THETA_EXCESS_K)})\n"
                 f"{I18N.f(lang, 'src_hybrid')} -- {I18N.f(lang, 'run')} {RUNDATE} {RUN:02d}Z",
                 fontsize=12)
    cb = fig.colorbar(im, ax=ax, shrink=0.6, aspect=34, pad=0.01)
    cb.set_label(I18N.f(lang, "decoup_cb"), fontsize=9)
    fig.savefig(path, dpi=140, bbox_inches="tight"); plt.close(fig)


def fire_bridge_constants():
    """``(head-fire length m, convective fraction)`` -- the Byram-to-PFT bridge, from the
    daily runner that actually computed the margin.

    Read from ``tests/pyroconv_daily.py`` rather than restated here so the figure cannot end
    up labelled with a head-fire length the field was not computed at; that script owns them
    (and honours ``PYROCONV_HEADFIRE_M``). It parses ``sys.argv`` at import, so argv is
    neutralised for the duration -- this script's own arguments are not its arguments.
    """
    import importlib.util
    src = os.path.join(REPO, "tests", "pyroconv_daily.py")
    argv = sys.argv
    try:
        sys.argv = [src]
        spec = importlib.util.spec_from_file_location("_pyroconv_daily_consts", src)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return float(mod._HEADFIRE_LENGTH_M), float(mod._CONVECTIVE_FRACTION)
    except Exception as e:
        sys.stderr.write(f"[3day] head-fire constants unavailable ({e}); figure will not state them\n")
        return None, None
    finally:
        sys.argv = argv


def pft_margin_figure(path, lang):
    """Day x hour grid of the PFT margin: firepower / PyroCb Firepower Threshold.

    The third question the report answers -- *is there enough fire here to make a pyroCb in
    this specific column?* -- and the only one of the three whose threshold is computed per
    column rather than fixed. Logarithmic, diverging about the criterion at 1.0, so that the
    single fact a reader must not misread -- which side of 1 a cell is on -- is carried by
    the colour pivot rather than by reading a number off a bar.

    Three exclusions, each drawn differently on purpose: columns the classifier rejects are
    grey (no atmosphere to test against), cells with no burnable fuel are white (zero
    firepower -- a categorical no, not a small margin), and sea is masked.
    """
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, TwoSlopeNorm
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    _, lat0, lon0 = load_stack(DAYS[0], "potential")
    ext = [lon0.min(), lon0.max(), lat0.min(), lat0.max()]
    sea = sea_mask(lat0, lon0)
    prov = provinces()
    cmap = plt.get_cmap("RdYlBu_r").copy()
    cmap.set_bad(alpha=0.0)                      # unclassifiable land -> panel facecolor
    sea_cmap, nofuel_cmap = ListedColormap([SEA_COLOR]), ListedColormap([NOFUEL_COLOR])
    lo, hi = np.log10(PFT_MARGIN_VMIN), np.log10(PFT_MARGIN_VMAX)
    norm = TwoSlopeNorm(vcenter=0.0, vmin=lo, vmax=hi)          # log10(margin), pivot at 1.0

    fig, ax = plt.subplots(len(DAYS), len(HOURS), figsize=(2.05 * len(HOURS), 2.5 * len(DAYS)),
                           constrained_layout=True, squeeze=False)
    im = None
    for r, valid in enumerate(DAYS):
        mg, ok = load_diag_stack(valid, "pft_margin"), valid_stack(valid)
        gate, _, _ = load_stack(valid, "gated")
        for c, hour in enumerate(HOURS):
            axc = ax[r][c]
            axc.set_facecolor(NODATA_COLOR)
            m = np.where(ok[c], mg[c].astype(float), np.nan)
            nofuel = np.isfinite(m) & (m <= 0)
            # Clip into the drawn range before the log: below the floor the exact value is
            # not resolvable on this scale anyway, and above the ceiling the pivot is what
            # matters, not how far past it the cell went (the table carries the true max).
            lm = np.log10(np.clip(np.where(nofuel, np.nan, m), PFT_MARGIN_VMIN, PFT_MARGIN_VMAX))
            im = axc.imshow(np.ma.masked_invalid(lm), origin="upper", extent=ext, cmap=cmap,
                            norm=norm, aspect="auto", interpolation="nearest")
            axc.imshow(np.ma.masked_where(~nofuel, np.zeros_like(m)), origin="upper", extent=ext,
                       cmap=nofuel_cmap, aspect="auto", interpolation="nearest")
            # Ring every cell that meets the criterion. Without this the map fails at the one
            # job it has: the field spans three decades below 1 and the cells at or above it
            # are a handful of 2 km pixels, so on a printed panel the answer to question 3 is
            # invisible inside the colour ramp that is supposed to convey it.
            # ...but split by whether the cell also clears question 2. The two criteria are
            # not nested: a column whose PFT has collapsed to a few GW can be "passed" by a
            # fire far too weak to raise any pyroCu at all, and marking those the same way
            # as a genuine candidate would be the most misleading thing this figure could do.
            hit = np.isfinite(m) & (m >= 1.0)
            gpass = hit & (gate[c] > 0)
            gfail = hit & ~gpass
            if gpass.any():
                yy, xx = np.nonzero(gpass)
                axc.scatter(lon0[xx], lat0[yy], s=15, facecolors="none", edgecolors="black",
                            linewidths=0.7, zorder=6)
            if gfail.any():
                yy, xx = np.nonzero(gfail)
                axc.scatter(lon0[xx], lat0[yy], s=11, marker="x", color="0.35",
                            linewidths=0.6, zorder=5)
            if sea is not None and sea.shape == m.shape:
                axc.imshow(np.ma.masked_where(~sea, np.zeros_like(m)), origin="upper",
                           extent=ext, cmap=sea_cmap, aspect="auto", interpolation="nearest")
            if prov is not None:
                prov.plot(ax=axc, color="0.25", linewidth=0.4)
                axc.set_xlim(ext[0], ext[1]); axc.set_ylim(ext[2], ext[3])
            axc.set_xticks([]); axc.set_yticks([])
            if r == 0:
                axc.set_title(f"{hour:02d}Z", fontsize=9)
            if c == 0:
                axc.set_ylabel(f"{valid}\n(+{r}d)", fontsize=9)
    hf, cf = fire_bridge_constants()
    bridge = "" if hf is None else f" ({I18N.f(lang, 'margin_bridge', m=hf, c=cf)})"
    fig.suptitle(f"{I18N.f(lang, 'margin_title')}{bridge}\n"
                 f"{I18N.f(lang, 'src_hybrid')} -- {I18N.f(lang, 'run')} {RUNDATE} {RUN:02d}Z",
                 fontsize=12)
    ticks = [t for t in range(int(lo), int(hi) + 1)]
    cb = fig.colorbar(im, ax=ax, shrink=0.6, aspect=34, pad=0.01, ticks=ticks)
    cb.ax.set_yticklabels([("1" if t == 0 else f"$10^{{{t}}}$") for t in ticks])
    cb.set_label(I18N.f(lang, "margin_cb"), fontsize=9)
    cb.ax.axhline(0.0, color="0.1", linewidth=1.6)                # the criterion, drawn on the bar
    fig.legend(handles=[Line2D([], [], marker="o", linestyle="none", markersize=7,
                               markerfacecolor="none", markeredgecolor="black",
                               label=I18N.f(lang, "margin_hit")),
                        Line2D([], [], marker="x", linestyle="none", markersize=6, color="0.35",
                               label=I18N.f(lang, "margin_miss")),
                        Patch(facecolor=NOFUEL_COLOR, edgecolor="0.4",
                              label=I18N.f(lang, "nofuel")),
                        Patch(facecolor=NODATA_COLOR, edgecolor="0.4",
                              label=I18N.f(lang, "nc"))],
               loc="lower center", ncol=4, fontsize=8.5, frameon=False, bbox_to_anchor=(0.5, -0.05))
    fig.savefig(path, dpi=140, bbox_inches="tight"); plt.close(fig)


def pft_margin_table():
    """Per-day daytime PFT-margin summary over the classifiable, burnable land.

    Rows are (day, hour, % grid with firepower, % of that at margin >= 1, >= 0.1, domain max).
    The denominator is deliberately the *burnable and classifiable* share rather than the
    domain: a percentage of Tuscany would be dominated by the cells that carry no fuel and can
    never contribute, and would fall towards zero for reasons that have nothing to do with the
    forecast. The carried share makes that denominator visible.
    """
    rows = []
    for valid in DAYS:
        mg, ok = load_diag_stack(valid, "pft_margin"), valid_stack(valid)
        gate, _, _ = load_stack(valid, "gated")
        for hi, h in enumerate(HOURS):
            if h not in DAYTIME:
                continue
            m = np.where(ok[hi], mg[hi].astype(float), np.nan)
            fire = np.isfinite(m) & (m > 0)          # a column to test, and fuel to test it with
            n = max(int(fire.sum()), 1)
            pct = [round(100.0 * int((fire & (m >= lv)).sum()) / n, 1) for lv in PFT_MARGIN_LEVELS]
            hit = fire & (m >= 1.0)
            rows.append((valid, h, round(100.0 * int(fire.sum()) / fire.size, 1), *pct,
                         round(float(m[fire].max()), 2) if fire.any() else float("nan"),
                         int(hit.sum()), int((hit & (gate[hi] > 0)).sum())))
    return rows


def decoupling_table():
    """Per-day daytime dry-decoupling summary over the classifiable land.

    Rows are (day, hour, %grid classifiable, %area >= 2, %area >= 3, domain max). The
    ratio percentages are of the cells classifiable at that hour, so they read as "of the
    columns the classifier accepted, this share would decouple by at least 2x" -- and the
    classifiable share is carried with them precisely because that denominator collapses
    around sunset: an evening row can read 100% off a handful of surviving columns, which
    is a statement about those columns, not about Tuscany.
    """
    rows = []
    for valid in DAYS:
        dec, ok = load_diag_stack(valid, "decoupling"), valid_stack(valid)
        for hi, h in enumerate(HOURS):
            if h not in DAYTIME:
                continue
            d = np.where(ok[hi], dec[hi], np.nan)
            fin = np.isfinite(d)
            n = max(int(fin.sum()), 1)
            pct = [round(100.0 * int((fin & (d >= lv)).sum()) / n, 1) for lv in DECOUP_LEVELS]
            rows.append((valid, h, round(100.0 * int(fin.sum()) / fin.size, 1), *pct,
                         round(float(d[fin].max()), 1) if fin.any() else float("nan")))
    return rows


def dist_table(kind):
    """Per-day, per-hour class % over land (land = cells ever classed > 0 across the 3 days).

    The trailing column is the share that could **not** be classified, which used to be
    silently folded into class 0 and is what made the collapsed evening hours read as a
    forecast of surface plumes. Each row therefore sums to 100 across classes 0-4 plus n/c.
    """
    stacks = {v: load_stack(v, "potential")[0] for v in DAYS}
    land = np.zeros(next(iter(stacks.values()))[0].shape, bool)
    for s in stacks.values():
        land |= (s > 0).any(0)
    n = max(int(land.sum()), 1)
    rows = []
    for valid in DAYS:
        st = load_stack(valid, kind)[0]
        for hi, h in enumerate(HOURS):
            pct = [round(100 * ((st[hi] == c) & land).sum() / n, 1) for c in range(5)]
            pct.append(round(100 * ((st[hi] < 0) & land).sum() / n, 1))
            rows.append((valid, h, pct))
    return rows, n


def build_report(png_pot, png_gate, png_margin, png_decoup, lang):
    """Write the .md and build the .pdf for one language, from the same rasters.

    The two editions are generated in the same call chain from one set of computed tables, so
    they cannot drift apart on a number: only the prose and the column headers come from the
    language table. Numbers are formatted once, here.
    """
    sfx = I18N.SUFFIX[lang]
    md = os.path.join(OUT, f"forecast_3day_{RUNDATE}{sfx}.md")
    pdf = os.path.join(OUT, f"forecast_3day_{RUNDATE}{sfx}.pdf")
    # The parameter is named _key, not k: the prose blocks take a kwarg called k (the
    # reference plume excess), and a shorter name here would shadow it.
    def T(_key, **kw):
        return I18N.t(lang, _key, **kw)
    rows, n = dist_table("gated")
    hf, cf = fire_bridge_constants()
    mlines = [T("mtbl_head"), "|:--|:--|--:|--:|--:|--:|--:|--:|"]
    for valid, h, fire, p1, p01, mx, n1, nboth in pft_margin_table():
        mlines.append(f"| {valid} | {h:02d}Z | {fire} | {p1} | {p01} | {mx} | {n1} | **{nboth}** |")
    mtbl = "\n".join(mlines)
    lines = [T("ctbl_head"), "|:--|:--|--:|--:|--:|--:|--:|--:|"]
    for valid, h, pct in rows:
        if h in DAYTIME:                # daytime rows only, to keep the table compact
            lines.append(f"| {valid} | {h:02d}Z | " + " | ".join(str(x) for x in pct) + " |")
    tbl = "\n".join(lines)
    dlines = [T("dtbl_head"), "|:--|:--|--:|--:|--:|--:|"]
    for valid, h, cls, p2, p3, mx in decoupling_table():
        dlines.append(f"| {valid} | {h:02d}Z | {cls} | {p2} | {p3} | {mx} |")
    dtbl = "\n".join(dlines)
    hours_txt = ", ".join(f"{h:02d}Z" for h in HOURS)
    code_ver, gen_at = provenance()
    abl = int(ABL_MIN_M)
    with open(md, "w") as f:
        f.write(f"""---
title: "{T('title_3day', rundate=RUNDATE, run=RUN)}"
subtitle: "{T('subtitle_3day', d0=DAYS[0], d2=DAYS[2])}"
geometry: a4paper, landscape, margin=1.1cm
fontsize: 9pt
lang: {lang}
header-includes: |
  \\usepackage{{float}}
  \\floatplacement{{figure}}{{H}}
---

## {T('h_intro_3day')}

{T('p_intro_3day', rundate=RUNDATE, run=RUN, hours=hours_txt)}

### {T('h_questions')}

{T('questions_intro')}

{T('questions_table')}

{T('questions_note')}

{T('questions_fourth')}

## {T('h1')}

![{T('cap1')}]({os.path.basename(png_pot)}){{width=100%}}

{T('p1')}

## {T('h2')}

![{T('cap2')}]({os.path.basename(png_gate)}){{width=100%}}

{T('p2')}

## {T('h3')}

![{T('cap3')}]({os.path.basename(png_margin)}){{width=100%}}

{T('p3a')}

{T('p3b', m=hf or 0.0, c=cf or 0.0)}

{T('p3c')}

{T('p3d')}

### {T('h_nested')}

{T('p_nested')}

## {T('h_margin_table')}

{mtbl}

{T('p_mtbl')}

## {T('h4', n=4)}

![{T('cap4')}]({os.path.basename(png_decoup)}){{width=100%}}

{T('p_decoup', k=_REFERENCE_THETA_EXCESS_K, abl=abl)}

## {T('h_dtbl')}

{dtbl}

{T('p_dtbl', abl=abl)}

## {T('h_ctbl', n=n)}

{tbl}

{T('p_ctbl', abl=abl)}

## {T('h_evening')}

{T('p_evening')}

{T('p_footer', ver=code_ver, at=gen_at)}
""")
    try:
        # Run from OUT: the image links are basenames, and pandoc resolves relative paths
        # against the working directory, not the .md's location -- from anywhere else the
        # figures silently drop out and the PDF builds anyway, just without them.
        subprocess.run(["pandoc", os.path.basename(md), "-o", os.path.basename(pdf),
                        "--pdf-engine=tectonic"], check=True, capture_output=True,
                       timeout=300, cwd=OUT)
        return pdf
    except Exception as e:
        sys.stderr.write(f"PDF build skipped for {lang} ({e})\n"); return None


def main():
    if RUN != 0:
        sys.exit(f"run {RUN:02d}Z: this product starts at 00Z of the run day, so day 0 would be "
                 f"missing its first {RUN} h. Use the 00Z run.")
    os.makedirs(OUT, exist_ok=True)
    if os.environ.get("SKIP_BASE_RUN") == "1":
        print("SKIP_BASE_RUN=1: re-stitching from the rasters already in "
              f"{OUT} -- the daily pipelines are not re-run.")
    else:
        for valid in DAYS:
            run_daily(valid)
    # One figure set and one report per language, from the identical rasters. The English
    # filenames keep their historical form so existing links still resolve.
    made = []
    for lang in I18N.LANGS:
        sfx = I18N.SUFFIX[lang]
        png_pot = os.path.join(OUT, f"forecast_3day_potential_{RUNDATE}{sfx}.png")
        png_gate = os.path.join(OUT, f"forecast_3day_gated_{RUNDATE}{sfx}.png")
        png_margin = os.path.join(OUT, f"forecast_3day_pft_margin_{RUNDATE}{sfx}.png")
        png_decoup = os.path.join(OUT, f"forecast_3day_decoupling_{RUNDATE}{sfx}.png")
        combined_figure("potential", png_pot, lang)
        combined_figure("gated", png_gate, lang)
        pft_margin_figure(png_margin, lang)
        decoupling_figure(png_decoup, lang)
        pdf = build_report(png_pot, png_gate, png_margin, png_decoup, lang)
        made += [png_pot, png_gate, png_margin, png_decoup] + ([pdf] if pdf else [])
    print(f"\nOK 3-day forecast -> {OUT}")
    for m in made:
        print(f"  {m}")


if __name__ == "__main__":
    main()
