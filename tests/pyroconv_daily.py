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
  PYROCONV_SOURCE     icon2i | hybrid   (default: icon2i)
  PYROCONV_HOURS      subset of valid hours, e.g. "12" or "9,12,15" (default: all 8)
Usage:  PYTHONPATH=src python tests/pyroconv_daily.py [YYYY-MM-DD] [run]

SOURCE=hybrid takes the atmospheric profile from ICON-EU model levels (the mixed-layer
dtheta/dz is a genuine measurement there, ~10 levels inside the mixed layer, vs a proxy on
ICON-2I's 5 pressure levels) and keeps ICON-2I's 2.2 km surface fields for the fuel gate.
It downloads one ICON-EU step per valid hour (~200 MB each), so use PYROCONV_HOURS to test.
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
    equilibrium_moisture_content, relative_humidity_from_dewpoint,
    fetch_icon2i_mistral, fetch_icon_eu, ICON2I_PROFILE_FIELDS, ICON_EU_MODEL_LEVELS,
    PYROCONVECTION_TYPES, PYROCONVECTION_TYPE_LEVEL, PYROCONVECTION_TYPE_COLOR,
    PYROCONVECTION_TYPE_LABEL,
)

# Shared compute core (also used by the Streamlit GUI). Add the repo root to the
# path so ``pyflam_gui`` is importable when running with ``PYTHONPATH=src``.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pyflam_gui.core.pyroconv import (
    read_icon2i_profile, profile_diagnostics, classify_profile, fli_grid,
    read_icon_eu, iconeu_diagnostics, regrid_diagnostics,
    lcp_fields as _core_lcp_fields)

warnings.simplefilter("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
HOURS = [0, 3, 6, 9, 12, 15, 18, 21]
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
# "icon2i" (default) uses ICON-2I 2.2 km's 5 pressure levels throughout. "hybrid" takes
# the profile diagnostics from ICON-EU model levels (~10 levels inside the mixed layer,
# so dtheta/dz is measurable via ml_method=fit_in_ml), regrids them onto the 2.2 km grid,
# and keeps ICON-2I's surface fields for the fuel gate.
SOURCE = os.environ.get("PYROCONV_SOURCE", "icon2i")
# Model token in output filenames, so hybrid and icon2i products do not clobber each other.
MODEL_TAG = "hybrid" if SOURCE == "hybrid" else "icon2i"
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
PROV_GEOJSON = (os.environ.get("PYROCONV_PROVINCES")
                or "/Users/cristianofoderi/DATI/boundaries/limits_IT_provinces.geojson")


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


def render(cats, lat, lon, tag):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch
    import rasterio; from rasterio.transform import from_origin
    os.makedirs(RASTERDIR, exist_ok=True)
    flip = lat[0] > lat[-1]; ext = [lon.min(), lon.max(), lat.min(), lat.max()]
    cmap = ListedColormap([PYROCONVECTION_TYPE_COLOR[t] for t in PYROCONVECTION_TYPES])
    norm = BoundaryNorm(np.arange(-0.5, 5.5, 1), cmap.N)
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
    src_label = ("ICON-EU model levels + ICON-2I 2.2 km gate" if SOURCE == "hybrid"
                 else "ICON-2I 2.2 km")
    fig.suptitle(f"Pyroconvection type -- {TAG_TITLE.get(tag, tag)}\n{src_label} -- "
                 f"Tuscany -- VALID {DATE} (run {RUNDATE} {RUN:02d}Z)", fontsize=11)
    leg = [Patch(facecolor=PYROCONVECTION_TYPE_COLOR[t], edgecolor="0.4",
                 label=f"{PYROCONVECTION_TYPE_LEVEL[t]}  {PYROCONVECTION_TYPE_LABEL[t]}")
           for t in PYROCONVECTION_TYPES]
    fig.legend(handles=leg, loc="lower center", ncol=5, fontsize=8.5, frameon=False,
               title="Pyroconvection class (0 = lowest activity -> 4 = highest)",
               bbox_to_anchor=(0.5, -0.12))
    png = os.path.join(OUTDIR, f"pyroconv_tuscany_{MODEL_TAG}_{tag}_{DATE}.png")
    fig.savefig(png, dpi=140, bbox_inches="tight"); plt.close(fig)
    dlon = float(abs(lon[1]-lon[0])); dlat = float(abs(lat[1]-lat[0]))
    tr = from_origin(lon.min()-dlon/2, lat.max()+dlat/2, dlon, dlat)
    for hi, hour in enumerate(HOURS):
        arr = cats[hi] if flip else cats[hi][::-1]
        with rasterio.open(os.path.join(RASTERDIR, f"pyroconv_{tag}_{hour:02d}Z.tif"), "w",
                           driver="GTiff", height=arr.shape[0], width=arr.shape[1], count=1,
                           dtype="int16", crs="EPSG:4326", transform=tr, nodata=-1) as d:
            d.write(arr.astype("int16"), 1)
    return png


def build_pdf(png_pot, png_gate, ladders, n_levels):
    md = os.path.join(OUTDIR, f"pyroconv_{MODEL_TAG}_{DATE}.md")
    pdf = os.path.join(OUTDIR, f"pyroconv_tuscany_{MODEL_TAG}_{DATE}.pdf")
    ladder_txt = ", ".join(ladders) if ladders else "none (no classifiable cell)"
    if SOURCE == "hybrid":
        src_title = "ICON-EU model levels + ICON-2I 2.2 km gate"
        heights_txt = ("per-cell heights from the ICON-EU model-level heights (HHL)")
        ml_rows = (
            "| ML dtheta/dz (least-squares fit, in mixed layer) | > 1.1e-3 K/m (stable) | "
            "Convection plume only -- no pyroCu |\n"
            "| ML dtheta/dz (measured on model levels) | <= 1.1e-3 K/m | "
            "Column is pyroCu-capable (bias ~-0.4e-4 K/m vs radiosondes) |")
        forcing_txt = (
            "Forcing: atmosphere from ICON-EU 6.5 km native model levels (DWD open data, "
            "CC-BY), lowest ~24 levels; ~{n} inside the mixed layer here, regridded to the "
            "ICON-2I 2.2 km grid. Surface fields and fuel gate from ICON-2I 2.2 km "
            "(MISTRAL / AgenziaItaliaMeteo). Classifier: pyflam.pyroconvection_type "
            "(bulk-Richardson ABL, Bolton LCL, measured mixed-layer dtheta/dz, cap "
            "gamma-theta, ABL-top RH; no surface CAPE)."
        ).format(n=n_levels)
        ml_para = (
"""The **mixed-layer stability** here is a genuine measurement, not a proxy. The ICON-EU
native model levels put ~10 levels inside the mixed layer (this run: {n} median), so the
mixed-layer dtheta/dz is a real least-squares fit across those levels. Validated against
IGRA radiosondes (JJA 12Z, period of record), the model-level fit cuts the gradient bias to
~-0.4e-4 K/m (from ~-4.7e-4 on ICON-2I's 5 pressure levels) and the Rib ABL bias to ~-70 m
(from ~-700 m). This is why the hybrid product exists: the ICON-2I 2.2 km open data cannot
resolve the mixed layer, and ICON-EU can.

The trade is horizontal resolution: the atmosphere is 6.5 km (regridded to the 2.2 km grid),
while the **fuel gate keeps ICON-2I's 2.2 km surface fields**, where fine terrain matters.
The mixed-layer gradient is now measured, but the 1.1e-3 K/m threshold still sits inside the
validated error bar (+/- ~2.6e-4), so treat class *counts* as indicative, not exact.""")
    else:
        src_title = "ICON-2I 2.2 km"
        heights_txt = "per-cell heights from the model geopotential (FI)"
        ml_rows = (
            "| ML stability *proxy*: parcel mixing depth (see Method) | depth < 455 m "
            "(\"stable\") | Convection plume only -- no pyroCu |\n"
            "| ML stability proxy | depth >= 455 m (= dtheta/dz <= 1.1e-3 K/m) | "
            "Column is pyroCu-capable (weak filter: spec. 0.29 inland) |")
        forcing_txt = (
            "Forcing: ICON-2I 2.2 km full-Italy GRIB (MISTRAL / AgenziaItaliaMeteo, CC-BY), "
            "Tuscany subset: geopotential, temperature, RH and wind on {n} pressure levels "
            "plus the 2 m / 10 m state, surface pressure and orography. Classifier: "
            "pyflam.pyroconvection_type (bulk-Richardson ABL, Bolton LCL, mixed-layer "
            "dtheta/dz proxy, cap gamma-theta, ABL-top RH; no surface CAPE)."
        ).format(n=n_levels)
        ml_para = (
"""The **mixed-layer stability** diagnostic is the weakest link in this product, and is a
*proxy*, not a measurement. Validated against IGRA radiosondes (JJA 12Z, period of record,
n=2835), a 5-pressure-level column contains only 0-2 model levels inside the mixed layer --
one or none in 92% of coastal columns -- and the lowest sits in the superadiabatic surface
layer. The mixed-layer dtheta/dz is therefore **not measurable from this archive**: fitting
it across the in-mixed-layer levels has no skill (Youden J ~ 0.00), and measuring theta up to
the Rib ABL top folds in the entrainment jump (Delta-theta), which Castellnou et al. (2022,
sec.2.1.1) treat as a variable *separate* from the gradient the ladder conditions on.

What is used instead is the **parcel mixing depth** as a proxy: at the 1.1e-3 K/m threshold
the criterion reduces to "well-mixed layer >= 455 m deep". It is the only candidate with
skill (J = 0.29 inland / 0.55 coastal, r = +0.50 against the radiosonde truth). It
**over-flags**: inland specificity is 0.29, so of the columns that are truly stable it still
calls ~71% pyroCu-capable (sensitivity 0.99 -- it rarely misses a capable column).

**Consequence: read the classes as a screening flag, not as calibrated counts.** The class
*totals* on these maps are not quantitatively trustworthy. The hybrid product
(`PYROCONV_SOURCE=hybrid`, ICON-EU model levels) measures this gradient properly; on this
5-level source, set `PYROCONV_ML_METHOD=surface_to_abl` or `mid_layer` for the superseded
measurements.""")
    ml_para = ml_para.format(n=n_levels)
    with open(md, "w") as f:
        f.write(f"""---
title: "Tuscany Pyroconvection-Type Forecast -- {src_title} -- VALID {DATE}"
subtitle: "3-hourly, 24 h. Run {RUNDATE} {RUN:02d}Z. Method after Castellnou et al. (2022), JGR-Atmos."
geometry: a4paper, landscape, margin=1.2cm
fontsize: 9pt
---

## How to read this product

Two panels are produced. The **fuel-gated** map is the expected, operationally
comparable product (the equivalent of the Catalan "tipus de piroconveccio" map):
a pyroCu/pyroCb class is assigned **only where a fire could actually reach >= 10
MW/m** of fireline intensity on the real Tuscany fuels. The **potential** map is an
unconditional **upper bound** -- it assumes a pyroCu-capable fire in *every* cell.
Use the gated panel for situational awareness; use the potential panel only to see
the atmospheric ceiling.

## Method and its limits (read this before using the classes)

The column is classified from a **real vertical profile**, not a standard atmosphere:
{heights_txt}, the mixing depth from the **bulk Richardson number** (first Rib >= 0.33
above 200 m AGL, referenced to the 2 m / 10 m state), the LCL from the exact
**Bolton (1980)** formula, and the humidity at the ABL top from the model's RH. The cap
gamma-theta is taken over ABL+200 m to ABL+1200 m.

{ml_para}

**Ladder actually used for this run: `{ladder_txt}`.** The full method has a fifth
diagnostic -- the distance from the ABL/LCL to the height of maximum wind shear -- which
this product does not yet compute, so the four-diagnostic ladder runs. The shear clause is a
*necessary* condition for the top class, so omitting it makes class 4 somewhat **easier** to
reach here than in the full method. Treat class 4 as an alert to inspect the column, not as a
calibrated probability.

The classes express atmospheric predisposition **given a fire of sufficient power**;
in the gated panel that power is computed, not assumed. They are paper-informed
thresholds, not locally validated ones.

## Expected pyroconvection type -- FUEL-GATED

![gated]({png_gate}){{width=100%}}

## Atmospheric potential -- UPPER BOUND (assumes a pyroCu-capable fire in every cell)

![potential]({png_pot}){{width=100%}}

## Class scale (low -> high pyroconvective activity)

| Level | Colour | Class | Meaning |
|:--:|:--|:--|:--|
| 0 | white | Surface plume | Buoyant smoke plume; no significant cloud development. |
| 1 | green | Convection plume | Plume penetrates a stable mixed layer; condensation possible, no pyroCu. |
| 2 | yellow | Overshooting pyroCu | Brief pyrocumulus; cloud base above the mixing height (LCL/ABL > 1). |
| 3 | orange | Resilient pyroCu | Persistent pyrocumulus in an unstable column (LCL/ABL < 1). |
| 4 | dark red | Deep pyroCu / pyroCb | Deep pyroconvection / pyrocumulonimbus; weak upper cap lets the plume deepen. |

## Classification thresholds (ladder in use: `{ladder_txt}`)

| Diagnostic | Threshold | Effect |
|:--|:--|:--|
{ml_rows}
| LCL / ABL ratio | 1.0 -- 1.60 | Overshooting pyroCu (brief) |
| LCL / ABL ratio | < 1.0 | Resilient pyroCu (persistent) |
| LCL / ABL ratio | <= 1.10 (+ conditions below) | Admissible for deep pyroCu / pyroCb |
| Cap gamma-theta (ABL+200 m -> ABL+1200 m) | <= 4.2e-3 K/m (weak cap) | Permits deepening to pyroCb |
| Cap gamma-theta | >= 4.8e-3 K/m (strong cap) | Inhibits deepening (resilient at most) |
| RH at the ABL top (mean, ABL +/- 150 m) | >= 80% | Required for classes 3 and 4 |
| Shear-maximum distance / ABL | <= 0.30 | Required for class 4 **in the 5-diagnostic ladder only** (not resolvable here -- see Method) |
| Fireline intensity (fuel gate, gated panel) | >= 10 MW/m | Minimum fire power for any pyroCu (Tedim et al. 2018) |
| ABL depth | < 600 m | Held at surface plume (mixing too shallow) |
| Usable pressure levels | < 4 | Cell not classified |

## Reference cases (Castellnou et al. 2022, Table 1)

The paper's labelled events, which anchor the published 3-diagnostic ladder
(`pyflam.pyroconvection_type(ladder="castellnou")`, still the library default and
regression-tested against these rows). The profile ladders used for this map are
stricter at the top: they additionally require a moist ABL top and an LCL close to it.

| Case | Observed type | LCL/ABL | ML dtheta/dz | gamma-theta (700-500) |
|:--|:--|:--:|:--|:--:|
| T21 | Convection plume | -- | stable | -- |
| SCQ32 | Overshooting pyroCu | > 1 | neutral/unstable | -- |
| M11 | Resilient pyroCu | < 1 | unstable | 4.2e-3 (resilient cap) |
| SCQ41 | pyroCu, not pyroCb | < 1 | unstable | 5.1e-3 (strong cap) |
| SCQ51 | Deep pyroCu / pyroCb | < 1 | unstable | 3.9e-3 (weak cap) |

{forcing_txt}
Gate: Rothermel + Cruz-2005 crown on the .lcp fuels with forecast moisture/wind.
Per-hour diagnostic rasters (ABL, LCL, LCL/ABL, ML dtheta/dz, cap, RH-top) accompany the classes.
Province borders: ISTAT-derived (openpolis geojson-italy). Generated by tests/pyroconv_daily.py.
""")
    try:
        subprocess.run(["pandoc", md, "-o", pdf, "--pdf-engine=tectonic"], check=True,
                       capture_output=True, timeout=300)
        return pdf
    except Exception as e:
        sys.stderr.write(f"PDF build skipped ({e}); PNGs + GeoTIFFs still written\n")
        return None


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
        for name in ("abl", "parcel_ml", "lcl", "lcl_ratio", "ml_grad", "gamma", "rh_top"):
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


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    sys.stderr.write(f"[pyroconv_daily] {DATE} {RUN:02d}Z -> {OUTDIR}  [source={SOURCE}]\n")
    # ICON-2I is always fetched: it supplies the surface fields for the fuel gate and the
    # 2.2 km grid everything renders on, whether or not the atmosphere comes from ICON-EU.
    files = fetch_icon2i_mistral(DT, run=RUN, cache_dir=CACHE,
                                 fields=ICON2I_PROFILE_FIELDS)
    d = read_icon2i_profile(files, (LAT1, LON0, LAT0, LON1), RUNDT, VALIDDT, HOURS)
    lat, lon, idx = d["lat"], d["lon"], d["idx"]
    lf, burn = lcp_fields(lat, lon)

    T2m_c = d["T2m"] - 273.15
    RH_sfc = relative_humidity_from_dewpoint(T2m_c, d["Td2m"] - 273.15)

    pot, gate, diags, ladders = [], [], [], set()
    for hi, si in enumerate(idx):
        if SOURCE == "hybrid":
            diag = eu_diag_for_hour(HOURS[hi], lat, lon)
            tag = f"EU model levels, ~{diag['n_levels']} in ML"
        else:
            diag = profile_diagnostics(d, si, ml_method=ML_METHOD)
            tag = f"2I levels={diag['n_levels']}"
        diags.append(diag)
        cls, used = classify_profile(diag, ladder=LADDER)
        pot.append(cls); ladders.update(used)
        if lf is not None:
            wsp = np.hypot(d["U10"][si], d["V10"][si])
            fli = fli_grid(T2m_c[si], RH_sfc[si], wsp, lf, burn)
            g, _ = classify_profile(diag, fli=fli, ladder=LADDER)
            gate.append(g)
        sys.stderr.write(f"  {HOURS[hi]:02d}Z  {tag}  "
                         f"ladder={','.join(sorted(used)) or '-'}\n")

    png_pot = render(np.stack(pot), lat, lon, "potential")
    png_gate = render(np.stack(gate), lat, lon, "gated") if gate else png_pot
    export_diagnostics(diags, lat, lon)
    pdf = build_pdf(png_pot, png_gate, sorted(ladders), diags[0]["n_levels"])
    print(f"OK {DATE} {RUN:02d}Z [ladder={','.join(sorted(ladders))}]: {png_pot}"
          + (f" | {pdf}" if pdf else ""))


if __name__ == "__main__":
    main()
