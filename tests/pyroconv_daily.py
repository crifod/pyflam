"""Daily Tuscany pyroconvection-type maps from ICON-2I 2.2 km (MISTRAL open data).

Self-contained, date-parametric runner for a scheduled job. Downloads today's
ICON-2I GRIB via pyflam.fetch_icon2i_mistral, classifies every native cell
(potential + fuel-gated), and writes PNG panels + 2 km GeoTIFFs + a dated PDF.

Classification runs the *profile* path: real per-cell heights from the model's
geopotential, a bulk-Richardson ABL, an exact (Bolton) LCL and moisture at the ABL
top. The archive publishes only 5-6 pressure levels, which cannot locate a
shear-maximum height, so the ladder used is the four-diagnostic one -- the fifth
(shear) diagnostic is dropped rather than fabricated from an interpolant. See
pyflam.atmosphere.shear_height_adaptive.

Config via env (all optional):
  PYROCONV_DATE   YYYY-MM-DD   (default: today, UTC)
  PYROCONV_RUN    0 | 12       (default: 0)
  PYROCONV_OUT    output dir   (default: <repo>/docs/daily)
  PYROCONV_CACHE  GRIB cache   (default: /tmp/pyflam_icon2i/<stamp>)
  PYROCONV_LCP    .lcp path    (default: the Tuscany canopy .lcp)
  PYROCONV_LADDER adaptive | noshear | shear | castellnou  (default: adaptive)
  PYROCONV_ML_METHOD  surface_to_parcel | surface_to_abl | mid_layer  (default: surface_to_parcel)
  PYROCONV_SOURCE     hybrid | icon2i   (default: hybrid, auto-falls back to icon2i)
  PYROCONV_HOURS      subset of valid hours, e.g. "12" or "9,12,15" (default: all 8)
  PYROCONV_PALETTE    pyflam | graf  (default: pyflam; "graf" = the Catalan Bombers key)
Usage:  PYTHONPATH=src python tests/pyroconv_daily.py [YYYY-MM-DD] [run]

SOURCE=hybrid (the default) takes the atmospheric profile from ICON-EU model levels -- the
mixed-layer dtheta/dz is a genuine measurement there (~10 levels inside the mixed layer),
against a mixing-depth proxy on ICON-2I's 5 pressure levels -- and keeps ICON-2I's 2.2 km
surface fields for the fuel gate. It downloads one ICON-EU step per valid hour (~200 MB each).
If the ICON-EU run is not yet published (or any step fails), the whole run **automatically
falls back to icon2i**, so a missing ICON-EU run never breaks the product; the filenames and
PDF then reflect the source actually used. Set SOURCE=icon2i to force the 5-pressure-level
product. Use PYROCONV_HOURS to test the (bandwidth-heavy) hybrid path on one hour.
"""
from __future__ import annotations

import os, sys, subprocess, warnings
from datetime import datetime, timezone
import numpy as np
import xarray as xr
from pyproj import Transformer

import pyflam
from pyflam import units, fuel_models
from pyflam.atmosphere import (
    critical_growth_rate_grid,
    equilibrium_moisture_content, relative_humidity_from_dewpoint,
    fetch_icon2i_mistral, fetch_icon_eu, ICON2I_PROFILE_FIELDS, ICON_EU_MODEL_LEVELS,
    PYROCONVECTION_TYPES, PYROCONVECTION_TYPE_LEVEL,
    PYROCONVECTION_TYPE_LABEL, pyroconvection_colors,
)

# Shared compute core (also used by the Streamlit GUI). Add the repo root to the
# path so ``pyflam_gui`` is importable when running with ``PYTHONPATH=src``.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pyflam_gui.core.pyroconv import (
    read_icon2i_profile, profile_diagnostics, classify_profile, fli_grid,
    read_icon_eu, iconeu_diagnostics, regrid_diagnostics,
    _REFERENCE_THETA_EXCESS_K, ML_FIT_MIN_PTS, PYROCONV_NODATA, ABL_MIN_M,
    lcp_fields as _core_lcp_fields)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                "scripts"))
import pyroconv_i18n as I18N        # report prose + figure labels, en/it

warnings.simplefilter("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
HOURS = [0, 3, 6, 9, 12, 15, 18, 21]
NODATA_COLOR = "#d9d9d9"        # land with no classifiable column (PYROCONV_NODATA)
LON0, LON1, LAT0, LAT1 = 9.6, 12.5, 42.2, 44.6
FLI_GATE_KW = 1.0e4
# Short titles per variant; the potential map is an unconditional atmospheric upper
# bound (assumes a pyroCu-capable fire everywhere), the gated map requires real fire
# power from the .lcp fuels -- the latter is comparable to the Catalan product.
TAG_TITLE = {
    "potential": "POTENTIAL -- atmosphere only (upper bound: assumes a pyroCu-capable fire in every cell)",
    "gated": "FUEL-GATED -- expected (only where fire power >= 10 MW/m on the Tuscany .lcp fuels)",
}
_TO3035 = Transformer.from_crs("EPSG:4326", "EPSG:3035", always_xy=True)

RUNDATE = (sys.argv[1] if len(sys.argv) > 1 else os.environ.get("PYROCONV_DATE")
           or datetime.now(timezone.utc).strftime("%Y-%m-%d"))
RUN = int(sys.argv[2] if len(sys.argv) > 2 else os.environ.get("PYROCONV_RUN", 0))
# Optional 3rd arg / PYROCONV_VALID: the VALID forecast day to render (default = run
# day). Lets a single run produce later forecast days, e.g. day+1, day+2 (<=72 h).
VALID = (sys.argv[3] if len(sys.argv) > 3 else os.environ.get("PYROCONV_VALID") or RUNDATE)
RUNDT = datetime.strptime(RUNDATE, "%Y-%m-%d").replace(hour=RUN)
VALIDDT = datetime.strptime(VALID, "%Y-%m-%d")
DATE = VALID                                   # labels / filenames use the valid day
DT = datetime.strptime(RUNDATE, "%Y-%m-%d")    # run day, for fetch
STAMP = f"{DT:%Y%m%d}{RUN:02d}"
OUTDIR = os.environ.get("PYROCONV_OUT") or os.path.join(REPO, "docs", "daily")
CACHE = os.environ.get("PYROCONV_CACHE") or f"/tmp/pyflam_icon2i/{STAMP}"
LCP = (os.environ.get("PYROCONV_LCP")
       or "/Users/cristianofoderi/DATI/FUEL_TOS/pyflam_canopy_tuscany/canopy_tuscany.lcp")
# Data source for the ATMOSPHERE (defined here because filenames below depend on it).
# "hybrid" (default) takes the profile diagnostics from ICON-EU model levels (~10 levels
# inside the mixed layer, so dtheta/dz is measurable via ml_method=fit_in_ml), regrids them
# onto the 2.2 km grid, and keeps ICON-2I's surface fields for the fuel gate. "icon2i" uses
# ICON-2I 2.2 km's 5 pressure levels throughout. SOURCE is what was *requested*;
# EFFECTIVE_SOURCE is what actually ran (main() downgrades hybrid -> icon2i if ICON-EU is
# unavailable). Filenames and PDF track EFFECTIVE_SOURCE.
SOURCE = os.environ.get("PYROCONV_SOURCE", "hybrid")
EFFECTIVE_SOURCE = SOURCE
MODEL_TAG = "hybrid" if EFFECTIVE_SOURCE == "hybrid" else "icon2i"
RASTERDIR = os.path.join(OUTDIR, f"rasters_{MODEL_TAG}_{DATE}")


def _apply_source(src):
    """Set the effective source and the filename/raster paths that depend on it."""
    global EFFECTIVE_SOURCE, MODEL_TAG, RASTERDIR
    EFFECTIVE_SOURCE = src
    MODEL_TAG = "hybrid" if src == "hybrid" else "icon2i"
    RASTERDIR = os.path.join(OUTDIR, f"rasters_{MODEL_TAG}_{DATE}")
# Which decision ladder to run. "adaptive" (default) uses the richest one the data
# support: on ICON-2I open data that is the 4-diagnostic ladder, because 5-6 pressure
# levels cannot locate a shear maximum. Set "castellnou" to reproduce the old
# 3-diagnostic product, or "shear" to force the full 5 (needs a richer profile).
LADDER = os.environ.get("PYROCONV_LADDER", "adaptive")
# How the mixed-layer dtheta/dz is measured. "surface_to_parcel" (default) takes it
# across the well-mixed layer, below the entrainment jump -- the quantity the
# Castellnou thresholds are defined on. "surface_to_abl" and "mid_layer" both fold
# the jump into the gradient and are kept only for comparison.
ML_METHOD = os.environ.get("PYROCONV_ML_METHOD", "surface_to_parcel")
# Optional comma-separated subset of valid hours, e.g. "12" or "9,12,15" -- handy for a
# quick single-hour product or for testing the (bandwidth-heavy) hybrid path.
if os.environ.get("PYROCONV_HOURS"):
    HOURS = [int(h) for h in os.environ["PYROCONV_HOURS"].split(",")]
# ICON-EU GRIB cache (per run; one subdir per forecast step).
EU_CACHE = os.environ.get("PYROCONV_EU_CACHE") or f"/tmp/pyflam_iconeu/{STAMP}"
# Tuscany province borders (ISTAT-derived, openpolis geojson-italy, EPSG:4326).
# Class colour palette: "pyflam" (default, ColorBrewer RdYlGn) or "graf" (the Catalan
# Bombers operational key, for side-by-side reading of the two products).
PALETTE = os.environ.get("PYROCONV_PALETTE", "pyflam")
COLORS = pyroconvection_colors(PALETTE)
PROV_GEOJSON = (os.environ.get("PYROCONV_PROVINCES")
                or "/Users/cristianofoderi/DATI/boundaries/limits_IT_provinces.geojson")


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


def _tuscany_provinces():
    """Tuscany province polygons clipped to the map bbox, or None if unavailable."""
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
    except Exception as e:
        sys.stderr.write(f"province borders skipped ({e})\n")
        return None


def lcp_fields(lat, lon):
    """Sample the Tuscany .lcp onto the weather grid (delegates to the shared core)."""
    if not os.path.exists(LCP):
        return None, None
    ls = pyflam.Landscape.from_lcp(LCP)
    return _core_lcp_fields(ls, lat, lon, _TO3035)


def render(cats, lat, lon, tag, lang="en", write_rasters=True):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch
    import rasterio; from rasterio.transform import from_origin
    os.makedirs(RASTERDIR, exist_ok=True)
    flip = lat[0] > lat[-1]; ext = [lon.min(), lon.max(), lat.min(), lat.max()]
    # Index -1 (PYROCONV_NODATA) gets its own grey: a column the ladder could not run on is
    # not a quiet one, and must not borrow class 0's colour.
    cmap = ListedColormap([NODATA_COLOR] + [COLORS[t] for t in PYROCONVECTION_TYPES])
    norm = BoundaryNorm(np.arange(-1.5, 5.5, 1), cmap.N)
    prov = _tuscany_provinces()
    fig, ax = plt.subplots(1, len(HOURS), figsize=(2.1*len(HOURS), 3.0),
                           constrained_layout=True, squeeze=False)
    ax = ax[0]                         # squeeze=False -> always a 1xN row; take the row
    for hi, hour in enumerate(HOURS):
        a = cats[hi][::-1] if flip else cats[hi]
        ax[hi].imshow(a, origin="lower", extent=ext, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
        if prov is not None:
            prov.plot(ax=ax[hi], color="0.15", linewidth=0.4)
            ax[hi].set_xlim(ext[0], ext[1]); ax[hi].set_ylim(ext[2], ext[3])
        ax[hi].set_title(f"{DATE} {hour:02d}Z", fontsize=8); ax[hi].set_xticks([]); ax[hi].set_yticks([])
    src_label = I18N.f(lang, "src_hybrid" if EFFECTIVE_SOURCE == "hybrid" else "src_icon2i")
    fig.suptitle(f"{I18N.f(lang, 'type_prefix')} -- {I18N.f(lang, tag + '_title_daily')}\n"
                 f"{src_label} -- {I18N.f(lang, 'tuscany')} -- {I18N.f(lang, 'valid')} {DATE} "
                 f"({I18N.f(lang, 'run')} {RUNDATE} {RUN:02d}Z)", fontsize=11)
    leg = [Patch(facecolor=COLORS[t], edgecolor="0.4",
                 label=f"{PYROCONVECTION_TYPE_LEVEL[t]}  {I18N.CLASS_LABEL[lang][t]}")
           for t in PYROCONVECTION_TYPES]
    leg.append(Patch(facecolor=NODATA_COLOR, edgecolor="0.4", label=I18N.f(lang, "nc_short")))
    fig.legend(handles=leg, loc="lower center", ncol=6, fontsize=8.5, frameon=False,
               title=I18N.f(lang, "class_legend_title"), bbox_to_anchor=(0.5, -0.12))
    png = os.path.join(OUTDIR,
                       f"pyroconv_tuscany_{MODEL_TAG}_{tag}_{DATE}{I18N.SUFFIX[lang]}.png")
    fig.savefig(png, dpi=140, bbox_inches="tight"); plt.close(fig)
    if not write_rasters:                # rasters carry no language; one pass writes them
        return png
    dlon = float(abs(lon[1]-lon[0])); dlat = float(abs(lat[1]-lat[0]))
    tr = from_origin(lon.min()-dlon/2, lat.max()+dlat/2, dlon, dlat)
    for hi, hour in enumerate(HOURS):
        arr = cats[hi] if flip else cats[hi][::-1]
        with rasterio.open(os.path.join(RASTERDIR, f"pyroconv_{tag}_{hour:02d}Z.tif"), "w",
                           driver="GTiff", height=arr.shape[0], width=arr.shape[1], count=1,
                           dtype="int16", crs="EPSG:4326", transform=tr, nodata=-1) as d:
            d.write(arr.astype("int16"), 1)
    return png


def render_decoupling(diags, lat, lon, lang="en"):
    """Render the dry-pyrocloud decoupling ratio (fireABL / ABL) as an 8-hour panel.

    A continuous heatmap, the DRY counterpart to the moist class map: how far a
    reference intense fire would grow its own boundary layer above the ambient ABL
    by sensible heat alone. Diagnostic only -- no class label.
    """
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(RASTERDIR, exist_ok=True)
    flip = lat[0] > lat[-1]; ext = [lon.min(), lon.max(), lat.min(), lat.max()]
    prov = _tuscany_provinces()
    fig, ax = plt.subplots(1, len(HOURS), figsize=(2.1*len(HOURS), 3.0),
                           constrained_layout=True, squeeze=False)
    ax = ax[0]
    cmap = plt.get_cmap("magma").copy(); cmap.set_bad("0.9")
    vmin, vmax = 1.0, 8.0
    im = None
    for hi, hour in enumerate(HOURS):
        r = np.asarray(diags[hi]["decoupling"], float)
        r = np.where(np.asarray(diags[hi]["valid"], bool), r, np.nan)
        a = r[::-1] if flip else r
        im = ax[hi].imshow(a, origin="lower", extent=ext, cmap=cmap, vmin=vmin, vmax=vmax,
                           aspect="auto", interpolation="nearest")
        if prov is not None:
            prov.plot(ax=ax[hi], color="0.15", linewidth=0.4)
            ax[hi].set_xlim(ext[0], ext[1]); ax[hi].set_ylim(ext[2], ext[3])
        ax[hi].set_title(f"{DATE} {hour:02d}Z", fontsize=8)
        ax[hi].set_xticks([]); ax[hi].set_yticks([])
    src_label = I18N.f(lang, "src_hybrid" if EFFECTIVE_SOURCE == "hybrid" else "src_icon2i")
    fig.suptitle(f"{I18N.f(lang, 'decoup_title_daily')}  "
                 f"({I18N.f(lang, 'decoup_ref', k=_REFERENCE_THETA_EXCESS_K)})\n"
                 f"{src_label} -- {I18N.f(lang, 'tuscany')} -- {I18N.f(lang, 'valid')} {DATE} "
                 f"({I18N.f(lang, 'run')} {RUNDATE} {RUN:02d}Z)", fontsize=11)
    cb = fig.colorbar(im, ax=ax, shrink=0.72, aspect=30, pad=0.01)
    cb.set_label(I18N.f(lang, "decoup_cb"), fontsize=8)
    png = os.path.join(OUTDIR,
                       f"pyroconv_tuscany_{MODEL_TAG}_decoupling_{DATE}{I18N.SUFFIX[lang]}.png")
    fig.savefig(png, dpi=140, bbox_inches="tight"); plt.close(fig)
    return png


PFT_MARGIN_VMIN, PFT_MARGIN_VMAX = 1.0e-3, 10.0
PFT_MARGIN_LEVELS = (1.0, 0.1)      # the criterion, and the within-one-decade band
NOFUEL_COLOR = "#ffffff"            # zero firepower: a categorical no, not a small margin


def render_pft_margin(diags, gate, lat, lon, lang="en"):
    """Render firepower / PyroCb Firepower Threshold as an 8-hour panel, or ``None``.

    The third question the product answers -- *is there enough fire here to make a pyroCb in
    this specific column?* -- and the only one whose threshold is computed per column rather
    than fixed. Returns ``None`` when the run has no margin field (no fuel gate, or an
    ICON-2I-only run that cannot form a PFT).

    Logarithmic and diverging about the criterion at 1.0, because the field spans four
    decades and the one fact a reader must not misread is which side of 1 a cell is on. Cells
    meeting the criterion are marked individually: at 2 km they are a few pixels, and inside a
    four-decade ramp the answer would otherwise be invisible.
    """
    if not gate or not all("pft_margin" in d for d in diags):
        return None
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, TwoSlopeNorm
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    flip = lat[0] > lat[-1]; ext = [lon.min(), lon.max(), lat.min(), lat.max()]
    prov = _tuscany_provinces()
    cmap = plt.get_cmap("RdYlBu_r").copy(); cmap.set_bad(NODATA_COLOR)
    nofuel_cmap = ListedColormap([NOFUEL_COLOR])
    lo, hi = np.log10(PFT_MARGIN_VMIN), np.log10(PFT_MARGIN_VMAX)
    norm = TwoSlopeNorm(vcenter=0.0, vmin=lo, vmax=hi)          # log10(margin), pivot at 1.0
    fig, ax = plt.subplots(1, len(HOURS), figsize=(2.1*len(HOURS), 3.0),
                           constrained_layout=True, squeeze=False)
    ax = ax[0]
    im = None
    for hi_, hour in enumerate(HOURS):
        m = np.where(np.asarray(diags[hi_]["valid"], bool),
                     np.asarray(diags[hi_]["pft_margin"], float), np.nan)
        g = np.asarray(gate[hi_])
        nofuel = np.isfinite(m) & (m <= 0)
        lm = np.log10(np.clip(np.where(nofuel, np.nan, m), PFT_MARGIN_VMIN, PFT_MARGIN_VMAX))
        a = (lambda x: x[::-1] if flip else x)
        im = ax[hi_].imshow(np.ma.masked_invalid(a(lm)), origin="lower", extent=ext, cmap=cmap,
                            norm=norm, aspect="auto", interpolation="nearest")
        ax[hi_].imshow(np.ma.masked_where(~a(nofuel), np.zeros_like(a(m))), origin="lower",
                       extent=ext, cmap=nofuel_cmap, aspect="auto", interpolation="nearest")
        # Split by whether the cell also clears the fuel gate. The two criteria are not
        # nested: a column whose PFT has collapsed to a few GW is "passed" by a fire far too
        # weak to raise any pyroCu at all, and marking those as candidates would be the most
        # misleading thing this panel could do.
        hit = np.isfinite(m) & (m >= 1.0)
        for sel, kw in ((hit & (g > 0), dict(s=15, facecolors="none", edgecolors="black",
                                             linewidths=0.7, zorder=6)),
                        (hit & ~(g > 0), dict(s=11, marker="x", color="0.35",
                                              linewidths=0.6, zorder=5))):
            if sel.any():
                yy, xx = np.nonzero(sel)
                ax[hi_].scatter(lon[xx], lat[yy], **kw)
        if prov is not None:
            prov.plot(ax=ax[hi_], color="0.15", linewidth=0.4)
            ax[hi_].set_xlim(ext[0], ext[1]); ax[hi_].set_ylim(ext[2], ext[3])
        ax[hi_].set_title(f"{DATE} {hour:02d}Z", fontsize=8)
        ax[hi_].set_xticks([]); ax[hi_].set_yticks([])
    src_label = I18N.f(lang, "src_hybrid" if EFFECTIVE_SOURCE == "hybrid" else "src_icon2i")
    fig.suptitle(f"{I18N.f(lang, 'margin_title_daily')}  "
                 f"({I18N.f(lang, 'margin_bridge', m=_HEADFIRE_LENGTH_M, c=_CONVECTIVE_FRACTION)})\n"
                 f"{src_label} -- {I18N.f(lang, 'tuscany')} -- {I18N.f(lang, 'valid')} {DATE} "
                 f"({I18N.f(lang, 'run')} {RUNDATE} {RUN:02d}Z)", fontsize=11)
    ticks = list(range(int(lo), int(hi) + 1))
    cb = fig.colorbar(im, ax=ax, shrink=0.72, aspect=30, pad=0.01, ticks=ticks)
    cb.ax.set_yticklabels([("1" if t == 0 else f"$10^{{{t}}}$") for t in ticks])
    cb.set_label(I18N.f(lang, "margin_cb"), fontsize=8)
    cb.ax.axhline(0.0, color="0.1", linewidth=1.6)              # the criterion, drawn on the bar
    fig.legend(handles=[
        Line2D([], [], marker="o", linestyle="none", markersize=7, markerfacecolor="none",
               markeredgecolor="black", label=I18N.f(lang, "margin_hit")),
        Line2D([], [], marker="x", linestyle="none", markersize=6, color="0.35",
               label=I18N.f(lang, "margin_miss")),
        Patch(facecolor=NOFUEL_COLOR, edgecolor="0.4", label=I18N.f(lang, "nofuel")),
        Patch(facecolor=NODATA_COLOR, edgecolor="0.4", label=I18N.f(lang, "nc_short"))],
        loc="lower center", ncol=4, fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.12))
    png = os.path.join(OUTDIR,
                       f"pyroconv_tuscany_{MODEL_TAG}_pft_margin_{DATE}{I18N.SUFFIX[lang]}.png")
    fig.savefig(png, dpi=140, bbox_inches="tight"); plt.close(fig)
    return png


def pft_margin_rows(diags, gate):
    """``(hour, %grid with firepower, %>=1, %>=0.1, max, n>=1, n>=1 and gate-passing)`` rows.

    The denominator for the two percentages is the burnable *and* classifiable share, not the
    domain: a percentage of Tuscany would be dominated by cells that carry no fuel and can
    never contribute, and would move for reasons unrelated to the forecast. The last column
    applies the joint criterion and is the count to act on.
    """
    if not gate or not all("pft_margin" in d for d in diags):
        return []
    rows = []
    for hi, hour in enumerate(HOURS):
        m = np.where(np.asarray(diags[hi]["valid"], bool),
                     np.asarray(diags[hi]["pft_margin"], float), np.nan)
        fire = np.isfinite(m) & (m > 0)
        n = max(int(fire.sum()), 1)
        pct = [round(100.0 * int((fire & (m >= lv)).sum()) / n, 1) for lv in PFT_MARGIN_LEVELS]
        hit = fire & (m >= 1.0)
        rows.append((hour, round(100.0 * int(fire.sum()) / fire.size, 1), *pct,
                     round(float(m[fire].max()), 2) if fire.any() else float("nan"),
                     int(hit.sum()), int((hit & (np.asarray(gate[hi]) > 0)).sum())))
    return rows


def build_pdf(png_pot, png_gate, ladders, n_levels, png_decoup=None, ml_fit_support=None,
              png_margin=None, margin_rows=None, lang="en"):
    """Write the .md and build the .pdf for one language, from one set of computed values.

    Called once per language from :func:`main`; the numbers are formatted here and only the
    prose comes from the language table, so the two editions cannot disagree about a value.
    """
    sfx = I18N.SUFFIX[lang]
    md = os.path.join(OUTDIR, f"pyroconv_{MODEL_TAG}_{DATE}{sfx}.md")
    pdf = os.path.join(OUTDIR, f"pyroconv_tuscany_{MODEL_TAG}_{DATE}{sfx}.pdf")
    def T(_key, **kw):
        return I18N.t(lang, _key, **kw)
    ladder_txt = ", ".join(ladders) if ladders else T("ladder_none")
    code_ver, gen_at = provenance()
    # The fifth (shear) diagnostic is available only where the profile resolves a
    # shear-maximum height: the ICON-EU model levels do, ICON-2I's 5 pressure levels do not,
    # and classify_profile(ladder="adaptive") decides that per cell. So report which case this
    # run actually hit -- asserting the pressure-level answer would misdescribe every hybrid
    # run, where the full ladder does run and class 4 does carry its shear clause.
    has_shear = "shear" in ladders
    which = ("partial" if has_shear and len(ladders) > 1 else "full" if has_shear else "none")
    shear_para = T(f"shear_{which}", lad=ladder_txt)
    shear_row = T(f"shear_row_{which}")
    hybrid = EFFECTIVE_SOURCE == "hybrid"
    src_title = I18N.f(lang, "src_hybrid" if hybrid else "src_icon2i")
    heights_txt = T("heights_hybrid" if hybrid else "heights_icon2i")
    ml_rows = T("ml_rows_hybrid" if hybrid else "ml_rows_icon2i")
    forcing_txt = T("forcing_hybrid" if hybrid else "forcing_icon2i", n=n_levels)
    if hybrid:
        sup = (T("sup_txt", pct=100.0 * ml_fit_support, pts=ML_FIT_MIN_PTS)
               if ml_fit_support is not None else T("sup_none"))
        ml_para = T("ml_para_hybrid", n=n_levels, sup=sup)
    else:
        ml_para = T("ml_para_icon2i", n=n_levels)
    abl = int(ABL_MIN_M)
    margin_block = ""
    if png_margin:
        mlines = [T("mtbl_head_daily"), "|:--|--:|--:|--:|--:|--:|--:|"]
        for hour, fire, p1, p01, mx, n1, nboth in (margin_rows or []):
            mlines.append(f"| {hour:02d}Z | {fire} | {p1} | {p01} | {mx} | {n1} | **{nboth}** |")
        hf, cf = _HEADFIRE_LENGTH_M, _CONVECTIVE_FRACTION
        margin_block = (
            f"## {T('h3')}\n\n"
            f"![{T('cap3')}]({png_margin}){{width=100%}}\n\n"
            f"{T('p3a')}\n\n{T('p3b', m=hf, c=cf)}\n\n{T('p3c')}\n\n"
            f"### {T('h_nested')}\n\n{T('p_nested')}\n\n"
            + "\n".join(mlines) + "\n\n" + T("p_mtbl") + "\n\n")
    decoup_block = ""
    if png_decoup:
        decoup_block = (
            f"## {T('h4', n=4 if png_margin else 3)}\n\n"
            f"![{T('cap4')}]({png_decoup}){{width=100%}}\n\n"
            f"{T('p_decoup', k=_REFERENCE_THETA_EXCESS_K, abl=abl)}\n\n")
    with open(md, "w") as f:
        f.write(f"""---
title: "{T('title_daily', src=src_title, date=DATE)}"
subtitle: "{T('subtitle_daily', rundate=RUNDATE, run=RUN)}"
geometry: a4paper, landscape, margin=1.2cm
fontsize: 9pt
lang: {lang}
header-includes: |
  \\usepackage{{float}}
  \\floatplacement{{figure}}{{H}}
---

## {T('h_howto')}

{T('questions_intro')}

{T('questions_table')}

{T('questions_note')}

## {T('h_method')}

{T('p_method', heights=heights_txt)}

{ml_para}

{shear_para}

{T('p_predisposition')}

## {T('h1')}

![{T('cap1')}]({png_pot}){{width=100%}}

{T('p1')}

## {T('h2')}

![{T('cap2')}]({png_gate}){{width=100%}}

{T('p2')}

{margin_block}{decoup_block}## {T('h_classscale')}

{T('classscale')}

## {T('h_thresholds', lad=ladder_txt)}

{T('thr_head')}
{ml_rows}
{T('thr_rest')}
{shear_row}
{T('thr_tail', abl=abl)}

## {T('h_refcases')}

{T('p_refcases')}

{T('refcases')}

{forcing_txt}
{T('p_footer_daily', ver=code_ver, at=gen_at)}
""")
    try:
        subprocess.run(["pandoc", md, "-o", pdf, "--pdf-engine=tectonic"], check=True,
                       capture_output=True, timeout=300)
        return pdf
    except Exception as e:
        sys.stderr.write(f"PDF build skipped for {lang} ({e}); PNGs + GeoTIFFs still written\n")
        return None


PEAK_ML_HOURS = (12, 15)          # mature mixed layer -- the regime the fit claim is about


def ml_fit_resolution(diags):
    """Peak-of-day mixed-layer resolution: ``(median levels, support fraction | None)``.

    The report quotes this to justify the hybrid path ("the model levels put ~N inside the
    mixed layer, so dtheta/dz is measured, not proxied"), and both the hour and the
    statistic have to match that claim:

    * **Hours.** Reading a single hour -- 00Z, as this once did -- samples the collapsed
      nocturnal layer and returns ~1, understating the resolution in the very sentence meant
      to establish it, on hours nobody classifies. Averaging the whole daytime span is no
      better: the growth (09Z), mature (12-15Z) and collapse (18Z) phases are different
      regimes, and their median describes no hour that occurred -- at 18Z the count is
      bimodal (collapsed columns against still-mixed ones), so a central value there is a
      number about nothing. The mature window is the regime the gradient claim concerns.
    * **Statistic.** The median says what a typical column has; the claim rests on the
      *worst-supported* ones, because a column below :data:`ML_FIT_MIN_PTS` returns a nan
      slope and drops out of the class map. ``ml_fit_support`` -- the share of land columns
      clearing that floor -- is what licenses the word "measured", and it is also the
      mechanism behind the evening coverage collapse. Both are reported.

    ``ml_fit_support`` is ``None`` on the ICON-2I pressure-level path, which proxies the
    gradient rather than fitting it. Falls back to all rendered hours when PYROCONV_HOURS
    excludes the peak window.
    """
    idx = [i for i, h in enumerate(HOURS) if h in PEAK_ML_HOURS] or list(range(len(diags)))
    n = int(np.median([diags[i]["n_levels"] for i in idx]))
    sup = [diags[i].get("ml_fit_support") for i in idx]
    sup = [s for s in sup if s is not None]
    return n, (float(np.mean(sup)) if sup else None)


def export_diagnostics(diags, lat, lon):
    """Write the per-hour profile diagnostics (ABL, LCL, ML dtheta/dz, cap, RH-top).

    The intermediate fields the class is built from, so a forecaster can see *why* a
    cell got its class rather than only the class.
    """
    import rasterio
    from rasterio.transform import from_origin
    os.makedirs(RASTERDIR, exist_ok=True)
    flip = lat[0] > lat[-1]
    dlon = float(abs(lon[1] - lon[0])); dlat = float(abs(lat[1] - lat[0]))
    tr = from_origin(lon.min() - dlon / 2, lat.max() + dlat / 2, dlon, dlat)
    for hi, hour in enumerate(HOURS):
        # The last four exist only on the model-level (hybrid) path; the ICON-2I pressure-level
        # fallback does not produce them, so missing keys are skipped rather than fatal.
        for name in ("abl", "parcel_ml", "lcl", "lcl_ratio", "ml_grad", "gamma", "rh_top",
                     "fireabl", "decoupling",
                     "residual_ml", "delta_theta", "firecape", "penetration",
                     "pft_gw", "z_fc", "delta_theta_fc", "u_ml", "abl_rib",
                     "fuel_load", "burnable_fraction", "firepower_gw", "pft_margin",
                     "crit_growth_ha_h"):
            if name not in diags[hi]:
                continue
            arr = np.asarray(diags[hi][name], "float32")
            arr = arr if flip else arr[::-1]
            path = os.path.join(RASTERDIR, f"diag_{name}_{hour:02d}Z.tif")
            with rasterio.open(path, "w", driver="GTiff", height=arr.shape[0],
                               width=arr.shape[1], count=1, dtype="float32",
                               crs="EPSG:4326", transform=tr, nodata=np.nan) as d:
                d.write(arr, 1)


def eu_diag_for_hour(h, lat, lon):
    """Hybrid path: ICON-EU model-level diagnostics for valid hour ``h``, on the 2.2 km grid.

    Fetches one ICON-EU forecast step (cached per step), computes the profile diagnostics
    from the native model levels (real in-mixed-layer dtheta/dz), and regrids them onto the
    ICON-2I grid the product renders on.
    """
    step = int((VALIDDT.replace(hour=h) - RUNDT).total_seconds() // 3600)
    files = fetch_icon_eu(DT, run=RUN, step=step, cache_dir=os.path.join(EU_CACHE, f"{step:03d}"))
    eud = read_icon_eu(files, (LAT1, LON0, LAT0, LON1), ICON_EU_MODEL_LEVELS)
    ediag = iconeu_diagnostics(eud, ml_method="fit_in_ml")
    return regrid_diagnostics(ediag, eud["lat"], eud["lon"], lat, lon)


def hybrid_diags_or_none(lat, lon):
    """All-hours ICON-EU diagnostics, or ``None`` if the ICON-EU run is unavailable.

    Built atomically: if *any* hour's ICON-EU fetch/read fails (e.g. the run is not yet
    published, or a step 404s), the whole attempt is abandoned and the caller falls back to
    ICON-2I -- a run is never a mix of two sources.
    """
    try:
        out = []
        for h in HOURS:
            diag = eu_diag_for_hour(h, lat, lon)
            out.append(diag)
            sys.stderr.write(f"  {h:02d}Z  EU model levels, ~{diag['n_levels']} in ML\n")
        return out
    except Exception as e:
        sys.stderr.write(f"[pyroconv_daily] hybrid (ICON-EU) unavailable ({e}); "
                         f"falling back to icon2i\n")
        return None


# --- fire-side conditioning ------------------------------------------------------------
# The pyroconvection classes describe a *plume*, and whether one forms depends on fire power
# as well as atmosphere. The POTENTIAL map posits a pyroCu-capable fire in every cell, which
# is why it saturates; these give the per-cell fire power to condition on instead.
#
# Byram intensity is power per metre of front (W/m); the PyroCb Firepower Threshold is a
# TOTAL power (W). Bridging them needs an active head-fire length. Rather than assume one,
# it is taken from the fire behaviour the Wildfire Data Portal publishes: L = (burn ratio) /
# ROS over its 20 fires with both, giving a median of 703 m (p25 398, p75 1300). The largest,
# Varnavas at 6.1 km, matches Tory & Kepert's "a head fire of about 5 km" for an extreme case.
_HEADFIRE_LENGTH_M = float(os.environ.get("PYROCONV_HEADFIRE_M", 700.0))
# Fraction of released heat entering the plume; the rest is radiated (Tory & Kepert app. D).
_CONVECTIVE_FRACTION = 0.7
# Available fuel load + burnable fraction on the forecast grid, from the 10 m FBFM40 map
# (scripts/fuel_load_10m.py). Optional: absent, the fuel diagnostics are simply not exported.
FUEL_LOAD_TIF = os.environ.get("PYFLAM_FUEL_LOAD_TIF",
                               os.path.join(REPO, "docs", "fuel_load_tuscany.tif"))


def sample_fuel_grid(lat, lon):
    """Available load (kg/m2) and burnable fraction sampled onto the (lat, lon) forecast grid.

    Returns ``(load, burnable_fraction)`` or ``(None, None)`` when the raster is absent. The
    raster is EPSG:3035; the forecast grid is lat/lon, so this reprojects by point sampling --
    adequate because the raster is already block-averaged to about the forecast cell size.
    """
    if not os.path.exists(FUEL_LOAD_TIF):
        return None, None
    try:
        import rasterio
        from pyproj import Transformer
        with rasterio.open(FUEL_LOAD_TIF) as src:
            tr = Transformer.from_crs("EPSG:4326", src.crs, always_xy=True)
            lon2d, lat2d = np.meshgrid(lon, lat)
            x, y = tr.transform(lon2d.ravel(), lat2d.ravel())
            rows, cols = rasterio.transform.rowcol(src.transform, x, y)
            rows = np.clip(np.asarray(rows), 0, src.height - 1)
            cols = np.clip(np.asarray(cols), 0, src.width - 1)
            load = src.read(1)[rows, cols].reshape(lat2d.shape)
            frac = src.read(2)[rows, cols].reshape(lat2d.shape)
        return load, frac
    except Exception as e:
        sys.stderr.write(f"[pyroconv_daily] fuel grid unavailable ({e})\n")
        return None, None


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    sys.stderr.write(f"[pyroconv_daily] {DATE} {RUN:02d}Z -> {OUTDIR}  [requested={SOURCE}]\n")
    # ICON-2I is always fetched: it supplies the surface fields for the fuel gate and the
    # 2.2 km grid everything renders on, whether or not the atmosphere comes from ICON-EU.
    files = fetch_icon2i_mistral(DT, run=RUN, cache_dir=CACHE,
                                 fields=ICON2I_PROFILE_FIELDS)
    d = read_icon2i_profile(files, (LAT1, LON0, LAT0, LON1), RUNDT, VALIDDT, HOURS)
    lat, lon, idx = d["lat"], d["lon"], d["idx"]
    lf, burn = lcp_fields(lat, lon)

    T2m_c = d["T2m"] - 273.15
    RH_sfc = relative_humidity_from_dewpoint(T2m_c, d["Td2m"] - 273.15)

    # Resolve the atmosphere source: hybrid if requested and ICON-EU is available, else
    # icon2i. Everything downstream (filenames, PDF, labels) follows EFFECTIVE_SOURCE.
    hybrid_diags = hybrid_diags_or_none(lat, lon) if SOURCE == "hybrid" else None
    _apply_source("hybrid" if hybrid_diags is not None else "icon2i")
    sys.stderr.write(f"[pyroconv_daily] effective source: {EFFECTIVE_SOURCE}\n")

    pot, gate, diags, ladders = [], [], [], set()
    fli_by_hour = []
    for hi, si in enumerate(idx):
        if hybrid_diags is not None:
            diag = hybrid_diags[hi]
        else:
            diag = profile_diagnostics(d, si, ml_method=ML_METHOD)
        diags.append(diag)
        cls, used = classify_profile(diag, ladder=LADDER)
        pot.append(cls); ladders.update(used)
        if lf is not None:
            wsp = np.hypot(d["U10"][si], d["V10"][si])
            fli = fli_grid(T2m_c[si], RH_sfc[si], wsp, lf, burn)
            fli_by_hour.append(fli)
            g, _ = classify_profile(diag, fli=fli, ladder=LADDER)
            gate.append(g)

    # Fire-side conditioning: per-cell firepower against the column's own PFT.
    fuel_load, burn_frac = sample_fuel_grid(lat, lon)
    for hi, diag in enumerate(diags):
        if fuel_load is not None:
            diag["fuel_load"] = fuel_load
            diag["burnable_fraction"] = burn_frac
            # Recompute the critical growth rate on the mapped load, replacing the default
            # constant the compute core had to assume (it does not see the fuel map).
            if "pft_gw" in diag:
                diag["crit_growth_ha_h"] = critical_growth_rate_grid(
                    diag["pft_gw"], fuel_load_kg_m2=fuel_load)
        if gate and hi < len(gate) and "pft_gw" in diag:
            fli_w_m = fli_by_hour[hi] * 1.0e3                    # kW/m -> W/m
            fp_gw = _CONVECTIVE_FRACTION * fli_w_m * _HEADFIRE_LENGTH_M / 1.0e9
            diag["firepower_gw"] = fp_gw
            with np.errstate(invalid="ignore", divide="ignore"):
                diag["pft_margin"] = fp_gw / np.where(diag["pft_gw"] > 0, diag["pft_gw"], np.nan)

    export_diagnostics(diags, lat, lon)
    n_lev, support = ml_fit_resolution(diags)
    mrows = pft_margin_rows(diags, gate)
    # One figure set and one report per language, off the identical arrays. The GeoTIFFs are
    # written on the first pass only -- a raster has no language, and rewriting them per
    # edition would only risk the two passes disagreeing.
    png_pot = pdf = None
    for li, lang in enumerate(I18N.LANGS):
        p_pot = render(np.stack(pot), lat, lon, "potential", lang, write_rasters=(li == 0))
        p_gate = (render(np.stack(gate), lat, lon, "gated", lang, write_rasters=(li == 0))
                  if gate else p_pot)
        p_margin = render_pft_margin(diags, gate, lat, lon, lang)
        p_decoup = render_decoupling(diags, lat, lon, lang)
        pdf_l = build_pdf(p_pot, p_gate, sorted(ladders), n_lev, png_decoup=p_decoup,
                          ml_fit_support=support, png_margin=p_margin, margin_rows=mrows,
                          lang=lang)
        if li == 0:
            png_pot, pdf = p_pot, pdf_l
    print(f"OK {DATE} {RUN:02d}Z [ladder={','.join(sorted(ladders))}]: {png_pot}"
          + (f" | {pdf}" if pdf else ""))


if __name__ == "__main__":
    main()
