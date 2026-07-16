"""3-day Tuscany pyroconvection forecast (hybrid) from a single 00Z run, as one report.

Runs the tested daily pipeline (tests/pyroconv_daily.py) for three valid days -- the run day
and the next two forecast days -- then stitches the three days into a combined 3-row figure
(potential + gated) and a markdown/PDF report with per-day class distributions.

Uses the day+1/day+2 forecast steps available in the same ICON-EU + ICON-2I 00Z run (ICON-EU to
+120 h, ICON-2I to +72 h). Output goes to its own folder.

Usage:  PYTHONPATH=src python scripts/forecast_3day_tuscany.py [YYYY-MM-DD run-day] [run]
Env:    PYROCONV_OUT (folder, default docs/forecast_<rundate>_3day), plus the usual
        PYROCONV_* knobs honoured by tests/pyroconv_daily.py.
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
                               PYROCONVECTION_TYPE_COLOR, PYROCONVECTION_TYPE_LABEL)

RUNDATE = sys.argv[1] if len(sys.argv) > 1 else datetime.now(timezone.utc).strftime("%Y-%m-%d")
RUN = int(sys.argv[2]) if len(sys.argv) > 2 else 0
RUNDT = datetime.strptime(RUNDATE, "%Y-%m-%d")
DAYS = [(RUNDT + timedelta(days=k)).strftime("%Y-%m-%d") for k in (0, 1, 2)]
HOURS = [0, 3, 6, 9, 12, 15, 18, 21]
OUT = os.environ.get("PYROCONV_OUT") or os.path.join(REPO, "docs", f"forecast_{RUNDATE}_3day")
LON0, LON1, LAT0, LAT1 = 9.6, 12.5, 42.2, 44.6            # Tuscany bbox (matches the daily runner)
CACHE = os.environ.get("PYROCONV_CACHE") or f"/tmp/pyflam_icon2i/{RUNDT:%Y%m%d}{RUN:02d}"
PROV_GEOJSON = (os.environ.get("PYROCONV_PROVINCES")
                or "/Users/cristianofoderi/DATI/boundaries/limits_IT_provinces.geojson")
SEA_COLOR = "#cfe4ef"


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
    """Invoke the daily pipeline for one valid day into OUT (hybrid default)."""
    env = dict(os.environ, PYROCONV_OUT=OUT, PYROCONV_VALID=valid)
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


def combined_figure(kind, path):
    """WRF-style day x hour grid: geographic panels, sea masked, province borders, legend."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch

    # class colours preceded by a sea colour (index -1 = sea).
    colors = [SEA_COLOR] + [PYROCONVECTION_TYPE_COLOR[t] for t in PYROCONVECTION_TYPES]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(-1.5, 5.5, 1), cmap.N)
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
                a = np.where(sea, -1, a)                 # paint sea cells with the sea colour
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
    title = ("FUEL-GATED -- expected (fire power >= 10 MW/m on the Tuscany fuels)"
             if kind == "gated"
             else "POTENTIAL -- atmospheric upper bound (assumes a pyroCu-capable fire everywhere)")
    fig.suptitle(f"Tuscany pyroconvection -- 3-day forecast -- {title}\n"
                 f"ICON-EU model levels + ICON-2I 2.2 km gate -- run {RUNDATE} {RUN:02d}Z",
                 fontsize=12)
    leg = [Patch(facecolor=PYROCONVECTION_TYPE_COLOR[t], edgecolor="0.4",
                 label=f"{PYROCONVECTION_TYPE_LEVEL[t]}  {PYROCONVECTION_TYPE_LABEL[t]}")
           for t in PYROCONVECTION_TYPES]
    fig.legend(handles=leg, loc="lower center", ncol=5, fontsize=9, frameon=False,
               title="Pyroconvection class (0 = lowest -> 4 = highest)",
               bbox_to_anchor=(0.5, -0.09))
    fig.savefig(path, dpi=140, bbox_inches="tight"); plt.close(fig)


def dist_table(kind):
    """Per-day, per-hour class % over land (land = cells ever classed > 0 across the 3 days)."""
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
            rows.append((valid, h, pct))
    return rows, n


def build_report(png_pot, png_gate):
    md = os.path.join(OUT, f"forecast_3day_{RUNDATE}.md")
    pdf = os.path.join(OUT, f"forecast_3day_{RUNDATE}.pdf")
    rows, n = dist_table("gated")
    lines = ["| Day | Hour | surface | convect | overshoot | resilient | deep |",
             "|:--|:--|--:|--:|--:|--:|--:|"]
    for valid, h, pct in rows:
        if h in (9, 12, 15, 18):        # daytime rows only, to keep the table compact
            lines.append(f"| {valid} | {h:02d}Z | " + " | ".join(str(x) for x in pct) + " |")
    tbl = "\n".join(lines)
    with open(md, "w") as f:
        f.write(f"""---
title: "Tuscany Pyroconvection -- 3-Day Forecast -- run {RUNDATE} {RUN:02d}Z"
subtitle: "Valid {DAYS[0]} .. {DAYS[2]}. Hybrid: ICON-EU model levels + ICON-2I 2.2 km fuel gate. Method after Castellnou et al. (2022)."
geometry: a4paper, landscape, margin=1.1cm
fontsize: 9pt
---

## 3-day forecast

Three valid days from the single {RUNDATE} {RUN:02d}Z run (day 0, +1, +2), using the day+1/+2
forecast steps in the same ICON-EU + ICON-2I datasets. Atmosphere from ICON-EU native model
levels (bulk-Richardson ABL, Bolton LCL, measured mixed-layer dtheta/dz, cap gamma-theta,
ABL-top RH, shear); the 10 MW/m fuel gate uses ICON-2I's 2.2 km surface fields on the Tuscany
.lcp fuels. See scripts/validation/README.md for the diagnostic validation (radiosondes + ERA5).

The **fuel-gated** map is the expected product; the **potential** map is the atmospheric
upper bound (assumes a pyroCu-capable fire in every cell). Read class *counts* as indicative.

## Fuel-gated (expected)

![gated]({os.path.basename(png_gate)}){{width=100%}}

## Potential (atmospheric upper bound)

![potential]({os.path.basename(png_pot)}){{width=100%}}

## Fuel-gated class distribution -- daytime, % of land cells ({n} land cells)

{tbl}

Classes: 0 surface plume, 1 convection plume, 2 overshooting pyroCu, 3 resilient pyroCu,
4 deep pyroCu/pyroCb. Generated by scripts/forecast_3day_tuscany.py.
""")
    try:
        subprocess.run(["pandoc", md, "-o", pdf, "--pdf-engine=tectonic"], check=True,
                       capture_output=True, timeout=300)
        return pdf
    except Exception as e:
        sys.stderr.write(f"PDF build skipped ({e})\n"); return None


def main():
    os.makedirs(OUT, exist_ok=True)
    for valid in DAYS:
        run_daily(valid)
    png_pot = os.path.join(OUT, f"forecast_3day_potential_{RUNDATE}.png")
    png_gate = os.path.join(OUT, f"forecast_3day_gated_{RUNDATE}.png")
    combined_figure("potential", png_pot)
    combined_figure("gated", png_gate)
    pdf = build_report(png_pot, png_gate)
    print(f"\nOK 3-day forecast -> {OUT}")
    print(f"  {png_pot}\n  {png_gate}" + (f"\n  {pdf}" if pdf else ""))


if __name__ == "__main__":
    main()
