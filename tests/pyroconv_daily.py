"""Daily Tuscany pyroconvection-type maps from ICON-2I 2.2 km (MISTRAL open data).

Self-contained, date-parametric runner for a scheduled job. Downloads today's
ICON-2I GRIB via pyflam.fetch_icon2i_mistral, classifies every native cell
(potential + fuel-gated), and writes PNG panels + 2 km GeoTIFFs + a dated PDF.

Config via env (all optional):
  PYROCONV_DATE   YYYY-MM-DD   (default: today, UTC)
  PYROCONV_RUN    0 | 12       (default: 0)
  PYROCONV_OUT    output dir   (default: <repo>/docs/daily)
  PYROCONV_CACHE  GRIB cache   (default: /tmp/pyflam_icon2i/<stamp>)
  PYROCONV_LCP    .lcp path    (default: the Tuscany canopy .lcp)
Usage:  PYTHONPATH=src python tests/pyroconv_daily.py [YYYY-MM-DD] [run]
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
    lcl_height_m, theta_kelvin, pyroconvection_type, equilibrium_moisture_content,
    relative_humidity_from_dewpoint, fetch_icon2i_mistral,
    PYROCONVECTION_TYPES, PYROCONVECTION_TYPE_LEVEL, PYROCONVECTION_TYPE_COLOR,
    PYROCONVECTION_TYPE_LABEL,
)

# Shared compute core (also used by the Streamlit GUI). Add the repo root to the
# path so ``pyflam_gui`` is importable when running with ``PYTHONPATH=src``.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pyflam_gui.core.pyroconv import (
    classify, fli_grid, lcp_fields as _core_lcp_fields)

warnings.simplefilter("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
HOURS = [0, 3, 6, 9, 12, 15, 18, 21]
LON0, LON1, LAT0, LAT1 = 9.6, 12.5, 42.2, 44.6
LEVELS = [850, 700, 500]
Z = {1000: 110.0, 850: 1457.0, 700: 3012.0, 500: 5574.0}
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
RASTERDIR = os.path.join(OUTDIR, f"rasters_{DATE}")
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


def _ds(path):
    ds = xr.open_dataset(path, engine="cfgrib", backend_kwargs={"indexpath": ""})
    la = ds["latitude"].values
    ds = (ds.sel(latitude=slice(LAT1, LAT0), longitude=slice(LON0, LON1)) if la[0] > la[-1]
          else ds.sel(latitude=slice(LAT0, LAT1), longitude=slice(LON0, LON1)))
    return ds, list(ds.data_vars)[0]


def load(files):
    T = {}
    for p in LEVELS:
        ds, v = _ds(files[f"T{p}"]); T[p] = ds[v].values - 273.15
    t2, v = _ds(files["T2M"]); T2m = t2[v].values - 273.15
    td, v = _ds(files["TD2M"]); Td = td[v].values - 273.15
    u, v = _ds(files["U10"]); U = u[v].values
    vv, vn = _ds(files["V10"]); V = vv[vn].values
    lat = t2["latitude"].values; lon = t2["longitude"].values
    sh = (t2["step"].values / np.timedelta64(1, "h")).astype(int)
    # forecast step (h from run) for each desired valid hour on the VALID day
    idx = []
    for h in HOURS:
        target = int((VALIDDT.replace(hour=h) - RUNDT).total_seconds() // 3600)
        m = np.where(sh == target)[0]
        if not m.size:
            raise ValueError(f"valid {VALID} {h:02d}Z needs forecast step +{target} h, "
                             f"beyond this run (max +{int(sh.max())} h)")
        idx.append(int(m[0]))
    return dict(lat=lat, lon=lon, idx=idx, T=T, T2m=T2m, Td=Td, U=U, V=V)


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
    fig, ax = plt.subplots(1, len(HOURS), figsize=(2.1*len(HOURS), 3.0), constrained_layout=True)
    for hi, hour in enumerate(HOURS):
        a = cats[hi][::-1] if flip else cats[hi]
        ax[hi].imshow(a, origin="lower", extent=ext, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
        if prov is not None:
            prov.plot(ax=ax[hi], color="0.15", linewidth=0.4)
            ax[hi].set_xlim(ext[0], ext[1]); ax[hi].set_ylim(ext[2], ext[3])
        ax[hi].set_title(f"{DATE} {hour:02d}Z", fontsize=8); ax[hi].set_xticks([]); ax[hi].set_yticks([])
    fig.suptitle(f"Pyroconvection type -- {TAG_TITLE.get(tag, tag)}\nICON-2I 2.2 km -- "
                 f"Tuscany -- VALID {DATE} (run {RUNDATE} {RUN:02d}Z)", fontsize=11)
    leg = [Patch(facecolor=PYROCONVECTION_TYPE_COLOR[t], edgecolor="0.4",
                 label=f"{PYROCONVECTION_TYPE_LEVEL[t]}  {PYROCONVECTION_TYPE_LABEL[t]}")
           for t in PYROCONVECTION_TYPES]
    fig.legend(handles=leg, loc="lower center", ncol=5, fontsize=8.5, frameon=False,
               title="Pyroconvection class (0 = lowest activity -> 4 = highest)",
               bbox_to_anchor=(0.5, -0.12))
    png = os.path.join(OUTDIR, f"pyroconv_tuscany_icon2i_{tag}_{DATE}.png")
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


def build_pdf(png_pot, png_gate):
    md = os.path.join(OUTDIR, f"pyroconv_{DATE}.md")
    pdf = os.path.join(OUTDIR, f"pyroconv_tuscany_icon2i_{DATE}.pdf")
    with open(md, "w") as f:
        f.write(f"""---
title: "Tuscany Pyroconvection-Type Forecast -- ICON-2I 2.2 km -- VALID {DATE}"
subtitle: "3-hourly, 24 h. Run {RUNDATE} {RUN:02d}Z. Method after Castellnou et al. (2022), JGR-Atmos."
geometry: a4paper, landscape, margin=1.2cm
fontsize: 9pt
---

## How to read this product

Two panels are produced. The **fuel-gated** map is the expected, operationally
comparable product (the equivalent of the Catalan "tipus de piroconveccio" map):
a pyroCu/pyroCb class is assigned **only where a fire could actually reach >= 10
MW/m** of fireline intensity on the real Tuscany fuels. The **potential** map is an
unconditional **upper bound** -- it assumes a pyroCu-capable fire in *every* cell,
so on a well-mixed summer afternoon it saturates to the high classes almost
everywhere. Use the gated panel for situational awareness; use the potential panel
only to see the atmospheric ceiling. (A standard-atmosphere free troposphere already
has gamma-theta ~= 3.3e-3 K/m, below the 4.0e-3 "deep" threshold -- hence the broad
red in the unconditional view; see the threshold table.)

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

## Classification thresholds (Castellnou et al. 2022)

| Diagnostic | Threshold | Effect |
|:--|:--|:--|
| Mixed-layer dtheta/dz (sfc -> 850 hPa) | > 1.1e-3 K/m (stable) | Convection plume only -- no pyroCu |
| Mixed-layer dtheta/dz | <= 1.1e-3 K/m (neutral/unstable) | Column is pyroCu-capable |
| LCL / ABL ratio | > 1 | Overshooting pyroCu (brief) |
| LCL / ABL ratio | < 1 | Resilient pyroCu (persistent) |
| Upper-cap gamma-theta (700 -> 500 hPa) | < 4.0e-3 K/m | Deepens to deep pyroCu / pyroCb |
| Fireline intensity (fuel gate, gated panel) | >= 10 MW/m | Minimum fire power for any pyroCu (Tedim et al. 2018) |
| ABL depth | < 600 m | Held at surface plume (mixing too shallow) |

## Reference cases (Castellnou et al. 2022, Table 1)

| Case | Observed type | LCL/ABL | ML dtheta/dz | gamma-theta (700-500) |
|:--|:--|:--:|:--|:--:|
| T21 | Convection plume | -- | stable | -- |
| SCQ32 | Overshooting pyroCu | > 1 | neutral/unstable | -- |
| M11 | Resilient pyroCu | < 1 | unstable | 4.2e-3 (resilient cap) |
| SCQ41 | pyroCu, not pyroCb | < 1 | unstable | 5.1e-3 (strong cap) |
| SCQ51 | Deep pyroCu / pyroCb | < 1 | unstable | 3.9e-3 (weak cap) |

Forcing: ICON-2I 2.2 km full-Italy GRIB (MISTRAL / AgenziaItaliaMeteo, CC-BY), Tuscany subset.
Classifier: pyflam.pyroconvection_type (LCL/ABL ratio, mixed-layer dtheta/dz, upper gamma-theta;
no surface CAPE). Gate: Rothermel + Cruz-2005 crown on the .lcp fuels with forecast moisture/wind.
Province borders: ISTAT-derived (openpolis geojson-italy). Generated by tests/pyroconv_daily.py.
""")
    try:
        subprocess.run(["pandoc", md, "-o", pdf, "--pdf-engine=tectonic"], check=True,
                       capture_output=True, timeout=300)
        return pdf
    except Exception as e:
        sys.stderr.write(f"PDF build skipped ({e}); PNGs + GeoTIFFs still written\n")
        return None


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    sys.stderr.write(f"[pyroconv_daily] {DATE} {RUN:02d}Z -> {OUTDIR}\n")
    files = fetch_icon2i_mistral(DT, run=RUN, cache_dir=CACHE)
    d = load(files); lat, lon, idx = d["lat"], d["lon"], d["idx"]
    lf, burn = lcp_fields(lat, lon)
    RH = relative_humidity_from_dewpoint(d["T2m"], d["Td"])
    pot, gate = [], []
    for si in idx:
        Tl = {p: d["T"][p][si] for p in LEVELS}
        pot.append(classify(d["T2m"][si], RH[si], Tl))
        if lf is not None:
            wsp = np.hypot(d["U"][si], d["V"][si])
            gate.append(classify(d["T2m"][si], RH[si], Tl, fli=fli_grid(d["T2m"][si], RH[si], wsp, lf, burn)))
    png_pot = render(np.stack(pot), lat, lon, "potential")
    png_gate = render(np.stack(gate), lat, lon, "gated") if gate else png_pot
    pdf = build_pdf(png_pot, png_gate)
    print(f"OK {DATE} {RUN:02d}Z: {png_pot}" + (f" | {pdf}" if pdf else ""))


if __name__ == "__main__":
    main()
