"""3-day Tuscany pyroconvection forecast (hybrid) from a single 00Z run, as one report.

Runs the tested daily pipeline (tests/pyroconv_daily.py) for three valid days -- the run day
and the next two forecast days -- then stitches the three days into three combined 3-row
figures and a markdown/PDF report:

  * **FUEL-GATED** -- the expected product (a class only where the .lcp fuels support
    >= 10 MW/m of fireline intensity);
  * **POTENTIAL** -- the atmospheric upper bound (assumes a pyroCu-capable fire everywhere);
  * **dry-pyrocloud decoupling** -- the continuous fireABL/ABL ratio for a reference intense
    fire, the DRY counterpart to the moist class ladder (diagnostic, no class label).

Every day is a **complete 24 h cycle, 3-hourly, always starting at 00Z of the run day**
(00, 03, 06, 09, 12, 15, 18, 21Z), on the Tuscany domain with ISTAT province borders.

Uses the day+1/day+2 forecast steps available in the same ICON-EU + ICON-2I 00Z run (ICON-EU to
+120 h, ICON-2I to +72 h). Output goes to its own folder.

Usage:  PYTHONPATH=src python scripts/forecast_3day_tuscany.py [YYYY-MM-DD run-day] [run]
Env:    PYROCONV_OUT (folder, default docs/forecast_<rundate>_3day), plus the usual
        PYROCONV_* knobs honoured by tests/pyroconv_daily.py. PYROCONV_HOURS is *not*
        honoured here -- the product is defined on the full 00-21Z cycle.
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
from pyflam_gui.core.pyroconv import ABL_MIN_M, _REFERENCE_FIRE_FLUX_W_M2

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
NODATA_COLOR = "0.92"                                     # land with no classifiable column
# Colour range of the dry-pyrocloud decoupling ratio (fireABL / ABL). 1 = the reference fire
# grows no boundary layer of its own above the ambient one; the daily product uses the same
# 1..8 range, so the panels are readable side by side.
DECOUP_VMIN, DECOUP_VMAX = 1.0, 8.0
# Decoupling levels the summary table counts (ratio >= x over the classified land).
DECOUP_LEVELS = (2.0, 3.0)


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
               PYROCONV_HOURS=",".join(str(h) for h in HOURS))
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


def combined_figure(kind, path):
    """WRF-style day x hour grid: geographic panels, sea masked, province borders, legend."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch

    # Class colours preceded by two sentinels: -2 = sea, -1 = PYROCONV_NODATA (a land column
    # the ladder could not run on). The sea sentinel moved from -1 to -2 when the classifier
    # started emitting -1 for unclassifiable cells -- sharing the value would have painted
    # every collapsed evening column as sea.
    colors = [SEA_COLOR, NODATA_COLOR] + [PYROCONVECTION_TYPE_COLOR[t]
                                          for t in PYROCONVECTION_TYPES]
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
    title = ("FUEL-GATED -- expected (fire power >= 10 MW/m on the Tuscany fuels)"
             if kind == "gated"
             else "POTENTIAL -- atmospheric upper bound (assumes a pyroCu-capable fire everywhere)")
    fig.suptitle(f"Tuscany pyroconvection -- 3-day forecast -- {title}\n"
                 f"ICON-EU model levels + ICON-2I 2.2 km gate -- run {RUNDATE} {RUN:02d}Z",
                 fontsize=12)
    leg = [Patch(facecolor=PYROCONVECTION_TYPE_COLOR[t], edgecolor="0.4",
                 label=f"{PYROCONVECTION_TYPE_LEVEL[t]}  {PYROCONVECTION_TYPE_LABEL[t]}")
           for t in PYROCONVECTION_TYPES]
    leg.append(Patch(facecolor=NODATA_COLOR, edgecolor="0.4",
                     label="n/c  not classifiable (no usable column)"))
    fig.legend(handles=leg, loc="lower center", ncol=6, fontsize=9, frameon=False,
               title="Pyroconvection class (0 = lowest -> 4 = highest)",
               bbox_to_anchor=(0.5, -0.09))
    fig.savefig(path, dpi=140, bbox_inches="tight"); plt.close(fig)


def decoupling_figure(path):
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
    fig.suptitle("Tuscany dry-pyrocloud decoupling -- 3-day forecast -- fireABL / ABL "
                 f"(reference {int(_REFERENCE_FIRE_FLUX_W_M2)} W/m^2 fire -- DIAGNOSTIC, no class)\n"
                 f"ICON-EU model levels + ICON-2I 2.2 km gate -- run {RUNDATE} {RUN:02d}Z",
                 fontsize=12)
    cb = fig.colorbar(im, ax=ax, shrink=0.6, aspect=34, pad=0.01)
    cb.set_label("fireABL / ABL   (1 = no decoupling; higher = deeper dry decoupling)", fontsize=9)
    fig.savefig(path, dpi=140, bbox_inches="tight"); plt.close(fig)


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


def build_report(png_pot, png_gate, png_decoup):
    md = os.path.join(OUT, f"forecast_3day_{RUNDATE}.md")
    pdf = os.path.join(OUT, f"forecast_3day_{RUNDATE}.pdf")
    rows, n = dist_table("gated")
    lines = ["| Day | Hour | surface | convect | overshoot | resilient | deep | n/c |",
             "|:--|:--|--:|--:|--:|--:|--:|--:|"]
    for valid, h, pct in rows:
        if h in DAYTIME:                # daytime rows only, to keep the table compact
            lines.append(f"| {valid} | {h:02d}Z | " + " | ".join(str(x) for x in pct) + " |")
    tbl = "\n".join(lines)
    dlines = ["| Day | Hour | % grid classifiable | % area ratio >= 2 | % area ratio >= 3 "
              "| max ratio |", "|:--|:--|--:|--:|--:|--:|"]
    for valid, h, cls, p2, p3, mx in decoupling_table():
        dlines.append(f"| {valid} | {h:02d}Z | {cls} | {p2} | {p3} | {mx} |")
    dtbl = "\n".join(dlines)
    hours_txt = ", ".join(f"{h:02d}Z" for h in HOURS)
    with open(md, "w") as f:
        f.write(f"""---
title: "Tuscany Pyroconvection -- 3-Day Forecast -- run {RUNDATE} {RUN:02d}Z"
subtitle: "Valid {DAYS[0]} .. {DAYS[2]}. Hybrid: ICON-EU model levels + ICON-2I 2.2 km fuel gate. Method after Castellnou et al. (2022)."
geometry: a4paper, landscape, margin=1.1cm
fontsize: 9pt
---

## 3-day forecast

Three valid days from the single {RUNDATE} {RUN:02d}Z run (day 0, +1, +2), using the day+1/+2
forecast steps in the same ICON-EU + ICON-2I datasets. Each day is a **complete 24 h cycle,
3-hourly, starting at 00Z of the run day** ({hours_txt}), so the three day-rows of every figure
are directly comparable and day 0 opens at the run's own initial time. Panels are the Tuscany
domain with ISTAT province borders; the sea is masked. Atmosphere from ICON-EU native model
levels (bulk-Richardson ABL, Bolton LCL, measured mixed-layer dtheta/dz, cap gamma-theta,
ABL-top RH, shear); the 10 MW/m fuel gate uses ICON-2I's 2.2 km surface fields on the Tuscany
.lcp fuels. See scripts/validation/README.md for the diagnostic validation (radiosondes + ERA5).

Three map series follow. The **fuel-gated** map is the expected product; the **potential** map
is the atmospheric upper bound (assumes a pyroCu-capable fire in every cell); the
**dry-pyrocloud decoupling** map is the dry counterpart to both, and carries no class label.
Read class *counts* as indicative.

## Fuel-gated (expected)

![gated]({os.path.basename(png_gate)}){{width=100%}}

## Potential (atmospheric upper bound)

![potential]({os.path.basename(png_pot)}){{width=100%}}

## Dry-pyrocloud decoupling -- DIAGNOSTIC (no class label)

![decoupling]({os.path.basename(png_decoup)}){{width=100%}}

The **decoupling ratio** fireABL / ABL is how high a reference intense fire
({int(_REFERENCE_FIRE_FLUX_W_M2)} W/m^2 convective flux) would grow its own boundary layer by
*sensible heat alone*, divided by the ambient ABL. It is the **dry** counterpart to the moist
class maps above (Castellnou et al. 2022; Castellnou Ribau et al. 2024): values well above 1
mark deep, hot, dry columns where a fire can punch through and decouple from the surface *even
where the moist ladder scores low*, which is exactly the situation the gated map under-reports.
Like the potential map it assumes a fire everywhere -- here a fixed reference flux rather than
a computed one -- so it is an upper bound, not an expectation, and no dry/moist LCL split is
applied (the +1 km literature offset is not supported by the GRAF prototype labels). fireABL
from `pyflam.atmosphere.fire_induced_abl_grid`, the parcel intersected with the real theta(z)
stack; cells the classifier rejects (no usable column, or ABL below {int(ABL_MIN_M)} m) are left
blank rather than painted with a bare quotient.

## Dry decoupling -- daytime, % of the classifiable land

{dtbl}

*% grid classifiable* is the share of the domain with a usable column at that hour (ABL above
{int(ABL_MIN_M)} m); the two ratio columns are percentages **of that share**, not of Tuscany.
The denominator collapses towards sunset as the mixed layer decays, so an 18Z row can read
100% off a small residual area -- read it together with the classifiable column.

## Fuel-gated class distribution -- daytime, % of land cells ({n} land cells)

{tbl}

Classes: 0 surface plume, 1 convection plume, 2 overshooting pyroCu, 3 resilient pyroCu,
4 deep pyroCu/pyroCb. **n/c** is land with no classifiable column (no usable profile, or an
ABL below {int(ABL_MIN_M)} m) -- distinct from class 0, which is the ladder's finding that the
column supports only a surface plume. Around sunset n/c takes most of the domain; those hours
carry no forecast, and must not be read as quiet ones. Per-hour class and diagnostic GeoTIFFs (ABL, LCL, LCL/ABL, ML dtheta/dz,
cap, RH-top, fireABL, decoupling) accompany this report in `rasters_hybrid_<day>/`.
Generated by scripts/forecast_3day_tuscany.py.
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
        sys.stderr.write(f"PDF build skipped ({e})\n"); return None


def main():
    if RUN != 0:
        sys.exit(f"run {RUN:02d}Z: this product starts at 00Z of the run day, so day 0 would be "
                 f"missing its first {RUN} h. Use the 00Z run.")
    os.makedirs(OUT, exist_ok=True)
    for valid in DAYS:
        run_daily(valid)
    png_pot = os.path.join(OUT, f"forecast_3day_potential_{RUNDATE}.png")
    png_gate = os.path.join(OUT, f"forecast_3day_gated_{RUNDATE}.png")
    png_decoup = os.path.join(OUT, f"forecast_3day_decoupling_{RUNDATE}.png")
    combined_figure("potential", png_pot)
    combined_figure("gated", png_gate)
    decoupling_figure(png_decoup)
    pdf = build_report(png_pot, png_gate, png_decoup)
    print(f"\nOK 3-day forecast -> {OUT}")
    print(f"  {png_pot}\n  {png_gate}\n  {png_decoup}" + (f"\n  {pdf}" if pdf else ""))


if __name__ == "__main__":
    main()
