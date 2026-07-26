"""3-day Tuscany forecast, enriched with per-province pyroconvection + FWI tables.

Wraps ``scripts/forecast_3day_tuscany.py`` (the tested hybrid pipeline: one 00Z run,
three valid days, ICON-EU model levels + ICON-2I 2.2 km gate) and adds, in a
**separate output folder**:

  * the three figure series it makes -- POTENTIAL, FUEL-GATED and the dry-pyrocloud
    DECOUPLING diagnostic -- 3 days x 8 hours, one complete 24 h cycle per day always
    starting at 00Z of the run day, on the Tuscany domain with ISTAT province borders;
  * a **per-province pyroconvection potential-metrics** table (zonal statistics of the
    class + diagnostic rasters over the ten Tuscany provinces, per valid day, including
    the peak daytime fireABL/ABL decoupling ratio); and
  * a **per-province FWI** table -- the full Canadian FWI System: the codes
    (FFMC / DMC / DC), the indices (ISI / BUI) and the final index (FWI) -- driven by an
    ERA5 spin-up of the slow drought codes (:func:`pyflam.fire_weather_history.fire_weather_spinup`,
    ~28 days ending the day before the run, SIR-rain-corrected, latency-bridged) carried
    forward through the three forecast days on the run's own noon (12Z) surface weather.

Everything (tables as CSV + Markdown, the two figures, and a combined PDF) lands in
``docs/forecast_<rundate>_3day_provinces/`` -- distinct from the base run folder.

Usage:  PYTHONPATH=src python scripts/forecast_3day_province_report.py [YYYY-MM-DD] [run]
Env:    PYROCONV_* knobs are passed through to the base pipeline; FWISPIN_DAYS (default 28),
        PYROCONV_PROVINCES (province geojson). Set SKIP_BASE_RUN=1 to reuse existing rasters.
"""
from __future__ import annotations

import csv
import os
import shutil
import subprocess
import sys
import warnings
from datetime import datetime, timedelta, timezone

import numpy as np
import rasterio
from rasterio.features import rasterize

warnings.simplefilter("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(REPO, "src"))
sys.path.insert(0, REPO)

from pyflam import fwi
from pyflam.atmosphere import (
    fetch_icon2i_mistral, ICON2I_PROFILE_FIELDS, relative_humidity_from_dewpoint,
    PYROCONVECTION_TYPE_LABEL)
from pyflam.fire_weather_history import (fire_weather_spinup, era5_daily_noon_history,
                                         blend_records)
from pyflam_gui.core.pyroconv import read_icon2i_profile, ABL_MIN_M, _REFERENCE_FIRE_FLUX_W_M2

RUNDATE = sys.argv[1] if len(sys.argv) > 1 else datetime.now(timezone.utc).strftime("%Y-%m-%d")
RUN = int(sys.argv[2]) if len(sys.argv) > 2 else 0
RUNDT = datetime.strptime(RUNDATE, "%Y-%m-%d").replace(hour=RUN)
DAYS = [(RUNDT + timedelta(days=k)).strftime("%Y-%m-%d") for k in (0, 1, 2)]
HOURS = [0, 3, 6, 9, 12, 15, 18, 21]
DAYTIME = [9, 12, 15, 18]                                  # peak-of-day window
NOON = 12
LON0, LON1, LAT0, LAT1 = 9.6, 12.5, 42.2, 44.6            # Tuscany bbox (matches the pipeline)
BASE_OUT = os.environ.get("PYROCONV_OUT") or os.path.join(REPO, "docs", f"forecast_{RUNDATE}_3day")
REPORT_OUT = os.path.join(REPO, "docs", f"forecast_{RUNDATE}_3day_provinces")
CACHE = os.environ.get("PYROCONV_CACHE") or f"/tmp/pyflam_icon2i/{RUNDT:%Y%m%d}{RUN:02d}"
PROV_GEOJSON = (os.environ.get("PYROCONV_PROVINCES")
                or "/Users/cristianofoderi/DATI/boundaries/limits_IT_provinces.geojson")
SPIN_DAYS = int(os.environ.get("FWISPIN_DAYS", 28))
# Local ERA5 daily-noon history, used as an offline drought-code anchor when the live CDS
# is unreachable: real reanalysis through its end date, then a persistence bridge to the run.
LOCAL_ERA5 = os.environ.get("FWI_LOCAL_ERA5") or os.path.join(REPO, "era5_hist_2026-06-05_2026-07-02.nc")
LOCAL_ERA5_START, LOCAL_ERA5_END = datetime(2026, 6, 5), datetime(2026, 7, 2)

# European Forest Fire danger classes on the final FWI index (EFFIS convention).
FWI_DANGER = [(5, "very low"), (11, "low"), (21, "moderate"),
              (38, "high"), (50, "very high"), (float("inf"), "extreme")]


def fwi_danger(v):
    for hi, label in FWI_DANGER:
        if v < hi:
            return label
    return "extreme"


# --------------------------------------------------------------------------- base run
def run_base_pipeline():
    """Produce the potential/gated/decoupling figures + class/diagnostic rasters."""
    if os.environ.get("SKIP_BASE_RUN") == "1":
        print(f"[report] SKIP_BASE_RUN=1 -> reusing {BASE_OUT}")
        missing = [p for p in (f"forecast_3day_potential_{RUNDATE}.png",
                               f"forecast_3day_gated_{RUNDATE}.png",
                               f"forecast_3day_decoupling_{RUNDATE}.png")
                   if not os.path.exists(os.path.join(BASE_OUT, p))]
        if missing:
            sys.exit(f"SKIP_BASE_RUN=1 but {BASE_OUT} is missing {', '.join(missing)} "
                     f"(a base run predating the decoupling series?). Re-run without it.")
        return
    print(f"[report] running base 3-day pipeline -> {BASE_OUT}")
    subprocess.run([sys.executable, os.path.join(HERE, "forecast_3day_tuscany.py"),
                    RUNDATE, str(RUN)], check=True, env=dict(os.environ, PYROCONV_OUT=BASE_OUT))


# --------------------------------------------------------------------------- provinces
def load_provinces():
    """Tuscany province GeoDataFrame clipped to the bbox (prov_name, geometry)."""
    import geopandas as gpd
    g = gpd.read_file(PROV_GEOJSON)
    if "reg_name" in g.columns:
        g = g[g["reg_name"] == "Toscana"]
    g = g.clip((LON0, LAT0, LON1, LAT1)).reset_index(drop=True)
    return g[["prov_name", "geometry"]]


def province_label_grid(sample_tif, gdf):
    """Rasterize provinces onto the class-raster grid: 0 = none, i+1 = gdf row i (north-up)."""
    with rasterio.open(sample_tif) as ds:
        transform, width, height = ds.transform, ds.width, ds.height
    shapes = [(geom, i + 1) for i, geom in enumerate(gdf.geometry)]
    return rasterize(shapes, out_shape=(height, width), transform=transform,
                     fill=0, dtype="int32", all_touched=True)


# --------------------------------------------------------------------------- rasters
def _raster_dir(valid):
    return os.path.join(BASE_OUT, f"rasters_hybrid_{valid}")


def load_class_stack(valid, kind):
    """(nhours, ny, nx) int16 class stack for 'potential'|'gated' (nodata -1)."""
    arrs = []
    for h in HOURS:
        with rasterio.open(os.path.join(_raster_dir(valid), f"pyroconv_{kind}_{h:02d}Z.tif")) as d:
            arrs.append(d.read(1))
    return np.stack(arrs)


def load_diag(valid, name, hour):
    """One diagnostic raster (float32, north-up) for a valid day + hour."""
    with rasterio.open(os.path.join(_raster_dir(valid), f"diag_{name}_{hour:02d}Z.tif")) as d:
        return d.read(1)


def classifiable(valid, hour):
    """Mask of columns the classifier accepted, rebuilt from the diagnostic rasters.

    Same definition as ``pyroconv.regrid_diagnostics``. The decoupling ratio is a bare
    fireABL/ABL quotient and stays finite in rejected columns, so it must be masked with
    this before any zonal statistic is taken.
    """
    abl = load_diag(valid, "abl", hour)
    return (np.isfinite(abl) & np.isfinite(load_diag(valid, "lcl_ratio", hour))
            & np.isfinite(load_diag(valid, "ml_grad", hour)) & (abl >= ABL_MIN_M))


# --------------------------------------------------------------------- pyroconv metrics
def pyroconv_metrics(gdf, lbl):
    """Per-province, per-day pyroconvection potential metrics from the class rasters.

    For each province the peak-of-day (09-18Z) statistics: the highest potential and
    gated class reached, the maximum daytime areal coverage of pyroCu (class >= 2) and of
    deep pyroCu/pyroCb (class == 4), the peak daytime dry-pyrocloud decoupling ratio
    (fireABL/ABL, diagnostic), and the noon column diagnostics (LCL/ABL, ABL depth,
    ABL-top RH) averaged over the province's classified land.
    """
    dt_idx = [HOURS.index(h) for h in DAYTIME]
    noon_i = HOURS.index(NOON)
    rows = []
    for valid in DAYS:
        pot = load_class_stack(valid, "potential")
        gate = load_class_stack(valid, "gated")
        lcl_ratio = load_diag(valid, "lcl_ratio", NOON)
        abl = load_diag(valid, "abl", NOON)
        rh_top = np.stack([load_diag(valid, "rh_top", h) for h in DAYTIME])
        decoup = np.stack([np.where(classifiable(valid, h), load_diag(valid, "decoupling", h),
                                    np.nan) for h in DAYTIME])
        for i, name in enumerate(gdf["prov_name"]):
            pm = (lbl == i + 1)
            land = pm & (pot[noon_i] >= 0)                  # classified land in the province
            n = max(int(land.sum()), 1)

            def cover(stack, thr):                          # max daytime %area with class >= thr
                return max(100.0 * ((stack[h] >= thr) & land).sum() / n for h in dt_idx)

            # -1 (nodata) never wins a max against a real class, so a peak of -1 means the
            # province had no classifiable column at any daytime hour -- report it as such
            # rather than letting it read as class 0.
            pot_peak = int(max((pot[h][land].max(initial=-1)) for h in dt_idx))
            gate_peak = int(max((gate[h][land].max(initial=-1)) for h in dt_idx))
            lr = float(np.nanmean(np.where(land, lcl_ratio, np.nan)))
            ab = float(np.nanmean(np.where(land, abl, np.nan)))
            rht = float(np.nanmin(np.where(land[None], rh_top, np.nan)))
            dc = np.where(land[None], decoup, np.nan)
            dcm = float(dc[np.isfinite(dc)].max()) if np.isfinite(dc).any() else float("nan")
            rows.append(dict(
                province=name, day=valid,
                pot_peak=pot_peak, pot_pyroCu=round(cover(pot, 2), 1),
                pot_deep=round(cover(pot, 4), 1),
                gate_peak=gate_peak, gate_pyroCu=round(cover(gate, 2), 1),
                decoup_max=round(dcm, 1),
                lcl_abl=round(lr, 2), abl_m=round(ab, 0), rh_top_min=round(rht, 0)))
    return rows


# ------------------------------------------------------------------------------- FWI
def read_noon_surface(valid):
    """ICON-2I 12Z surface weather (T2m, RH, wind km/h) + land mask for a valid day.

    Returns north-up arrays aligned with the class/diagnostic rasters, so the province
    label grid applies unchanged.
    """
    validdt = datetime.strptime(valid, "%Y-%m-%d")
    files = fetch_icon2i_mistral(RUNDT, run=RUN, cache_dir=CACHE, fields=ICON2I_PROFILE_FIELDS)
    d = read_icon2i_profile(files, (LAT1, LON0, LAT0, LON1), RUNDT, validdt, [NOON])
    si = d["idx"][0]
    flip = d["lat"][0] > d["lat"][-1]                       # north-first? then already north-up

    def north_up(a):
        return a if flip else a[::-1]

    t2m_c = north_up(d["T2m"][si] - 273.15)
    td2m_c = north_up(d["Td2m"][si] - 273.15)
    wind_kmh = north_up(np.hypot(d["U10"][si], d["V10"][si])) * 3.6
    rh = relative_humidity_from_dewpoint(t2m_c, td2m_c)
    land = north_up(d["frland"]) >= 0.5
    return t2m_c, rh, wind_kmh, land


def _persistence_bridge(state, last_rec, start_day, end_day):
    """Step FWI daily from ``start_day`` to ``end_day`` on the last ERA5 noon weather (rain=0).

    Bridges the gap between the local ERA5 anchor's end and the run day when the live CDS
    (which would supply that stretch) is down: a dry-spell persistence of the last observed
    noon temperature / RH / wind. Returns the carried :class:`FWIState` at ``end_day``.
    """
    day = start_day
    while day <= end_day:
        idx = fwi.FWISystem(state=state).step(
            temperature=last_rec["temperature"], relative_humidity=last_rec["relative_humidity"],
            wind_kmh=fwi.wind_ms_to_kmh(last_rec["wind_ms"]), rain_mm=0.0, month=day.month)
        state = fwi.FWIState(ffmc=idx.ffmc, dmc=idx.dmc, dc=idx.dc)
        day += timedelta(days=1)
    return state


def spinup_state_for_run(cen):
    """Antecedent FWIState for the run day, with the best drought source available.

    1. live CDS ERA5 spun up over ~SPIN_DAYS ending the day before the run (the requested
       path); on any failure ->
    2. the local ERA5 history file (real reanalysis through its end date) + a dry-spell
       persistence bridge to the run day; on failure ->
    3. spring-startup codes (indicative only).

    Returns ``(FWIState, method_text)`` where method_text documents which path ran.
    """
    end = RUNDT.date() - timedelta(days=1)
    try:
        st = fire_weather_spinup(RUNDT.date(), cen.y, cen.x, days=SPIN_DAYS,
                                 cache_path=os.path.join(REPORT_OUT, f"era5_spinup_{RUNDATE}.nc"))
        # fire_weather_spinup degrades to spring-startup internally on ERA5 failure; detect that
        if not (st.ffmc == fwi.STARTUP_FFMC and st.dmc == fwi.STARTUP_DMC and st.dc == fwi.STARTUP_DC):
            return st, (f"live ERA5 reanalysis, {SPIN_DAYS} d ending {end} (SIR-rain-corrected "
                        f"where the gauge window is current)")
        raise RuntimeError("live ERA5 returned spring-startup (CDS unavailable)")
    except Exception as e:
        sys.stderr.write(f"[report] live ERA5 spin-up unavailable ({e}); trying local ERA5\n")

    if os.path.exists(LOCAL_ERA5):
        recs = era5_daily_noon_history(cen.y, cen.x, LOCAL_ERA5_START.date(),
                                       LOCAL_ERA5_END.date(), cache_path=LOCAL_ERA5)
        anchor = fwi.spinup_state(blend_records(recs))
        bstart = LOCAL_ERA5_END.date() + timedelta(days=1)
        st = _persistence_bridge(anchor, recs[-1], bstart, end) if bstart <= end else anchor
        return st, (f"local ERA5 reanalysis {LOCAL_ERA5_START:%Y-%m-%d}..{LOCAL_ERA5_END:%Y-%m-%d} "
                    f"(anchor DMC {anchor.dmc:.0f} / DC {anchor.dc:.0f}) + dry-spell persistence "
                    f"bridge {bstart}..{end} -- **live CDS was unreachable at run time**")

    return fwi.FWIState(), ("spring-startup codes (FFMC 85 / DMC 6 / DC 15) -- **no ERA5 available; "
                            "DMC/DC/BUI/FWI are indicative only**")


def fwi_grids():
    """Per-day gridded FWI codes carried forward from the spin-up state, + land mask + method.

    The slow drought codes (DMC/DC) are spun up once at the Tuscany centroid (see
    :func:`spinup_state_for_run`) and broadcast as the day-0 antecedent state (near-uniform
    over the domain, per :func:`pyflam.fwi.gridded_fwi`); the full code grid is then stepped
    day by day on each valid day's 12Z weather (dry-forecast assumption, rain = 0).
    """
    cen = load_provinces().dissolve().geometry.iloc[0].centroid
    print(f"[report] FWI spin-up at Tuscany centroid {cen.y:.3f} N, {cen.x:.3f} E, "
          f"ending {DAYS[0]} - 1 d")
    spin, method = spinup_state_for_run(cen)
    print(f"[report] spin-up state: FFMC {spin.ffmc:.1f} / DMC {spin.dmc:.1f} / "
          f"DC {spin.dc:.1f}  [{method}]")
    prev = spin
    out, land0 = {}, None
    for valid in DAYS:
        t2m, rh, wind_kmh, land = read_noon_surface(valid)
        if land0 is None:
            land0 = land
        month = datetime.strptime(valid, "%Y-%m-%d").month
        codes = fwi.gridded_fwi(temperature=t2m, relative_humidity=rh, wind_kmh=wind_kmh,
                                month=month, previous=prev, rain_mm=0.0)
        out[valid] = codes
        prev = fwi.FWIState(ffmc=codes["ffmc"], dmc=codes["dmc"], dc=codes["dc"])
    return out, land0, spin, method


def fwi_metrics(gdf, lbl, grids, land):
    """Per-province, per-day province-mean of the six FWI quantities at 12Z."""
    rows = []
    for valid in DAYS:
        c = grids[valid]
        for i, name in enumerate(gdf["prov_name"]):
            m = (lbl == i + 1) & land
            if not m.any():
                m = (lbl == i + 1)
            def mean(k):
                return float(np.nanmean(c[k][m]))
            f = round(mean("fwi"), 1)
            rows.append(dict(
                province=name, day=valid,
                ffmc=round(mean("ffmc"), 1), dmc=round(mean("dmc"), 1), dc=round(mean("dc"), 1),
                isi=round(mean("isi"), 1), bui=round(mean("bui"), 1), fwi=f,
                danger=fwi_danger(f)))
    return rows


# ---------------------------------------------------------------------------- outputs
def write_csv(path, rows, fields):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def md_table(rows, cols, headers, aligns):
    head = "| " + " | ".join(headers) + " |"
    sep = "|" + "|".join(aligns) + "|"
    body = ["| " + " | ".join(str(r[c]) for c in cols) + " |" for r in rows]
    return "\n".join([head, sep] + body)


def build_report(pyro_rows, fwi_rows, spin, spin_method, png_pot, png_gate, png_decoup):
    md = os.path.join(REPORT_OUT, f"forecast_3day_provinces_{RUNDATE}.md")
    pdf = os.path.join(REPORT_OUT, f"forecast_3day_provinces_{RUNDATE}.pdf")
    classes = ", ".join(f"{k} {PYROCONVECTION_TYPE_LABEL[k]}" for k in sorted(PYROCONVECTION_TYPE_LABEL))
    pyro_md = md_table(
        pyro_rows,
        ["province", "day", "pot_peak", "pot_pyroCu", "pot_deep", "gate_peak", "gate_pyroCu",
         "decoup_max", "lcl_abl", "abl_m", "rh_top_min"],
        ["Province", "Day", "Pot peak", "Pot %pyroCu", "Pot %deep", "Gate peak", "Gate %pyroCu",
         "Decoup max", "LCL/ABL", "ABL m", "RH-top min"],
        [":--", ":--", "--:", "--:", "--:", "--:", "--:", "--:", "--:", "--:", "--:"])
    fwi_md = md_table(
        fwi_rows,
        ["province", "day", "ffmc", "dmc", "dc", "isi", "bui", "fwi", "danger"],
        ["Province", "Day", "FFMC", "DMC", "DC", "ISI", "BUI", "FWI", "Danger"],
        [":--", ":--", "--:", "--:", "--:", "--:", "--:", "--:", ":--"])
    hours_txt = ", ".join(f"{h:02d}Z" for h in HOURS)
    with open(md, "w") as f:
        f.write(f"""---
title: "Tuscany Pyroconvection + FWI -- 3-Day Provincial Forecast -- run {RUNDATE} {RUN:02d}Z"
subtitle: "Valid {DAYS[0]} .. {DAYS[2]}. Hybrid: ICON-EU model levels + ICON-2I 2.2 km gate. Pyroconvection after Castellnou et al. (2022); Canadian FWI System (Van Wagner 1987)."
geometry: a4paper, landscape, margin=1.0cm
fontsize: 8pt
---

## Overview

Three valid days from the single **{RUNDATE} {RUN:02d}Z** run (day 0, +1, +2), using the
day+1/+2 forecast steps in the same ICON-EU + ICON-2I datasets. Each day is a complete 24 h
cycle, 3-hourly, **starting at 00Z of the run day** ({hours_txt}); panels are the Tuscany
domain with ISTAT province borders. The atmospheric pyroconvection classification is the
hybrid product (ICON-EU native model levels for the profile -- bulk-Richardson ABL, Bolton
LCL, measured mixed-layer dtheta/dz, cap gamma-theta, ABL-top RH -- with the 10 MW/m fuel
gate on ICON-2I's 2.2 km surface fields over the Tuscany .lcp fuels). Three map series
follow; per-province tables summarise the pyroconvection potential, the dry-pyrocloud
decoupling and the Canadian fire-weather danger for each of the ten Tuscany provinces.

## Fuel-gated pyroconvection (expected)

![gated]({os.path.basename(png_gate)}){{width=100%}}

## Potential pyroconvection (atmospheric upper bound)

![potential]({os.path.basename(png_pot)}){{width=100%}}

## Dry-pyrocloud decoupling -- DIAGNOSTIC (no class label)

![decoupling]({os.path.basename(png_decoup)}){{width=100%}}

The decoupling ratio fireABL / ABL is how high a reference intense fire
({int(_REFERENCE_FIRE_FLUX_W_M2)} W/m^2 convective flux) would grow its own boundary layer by
*sensible heat alone*, divided by the ambient ABL. It is the **dry** counterpart to the two
class maps (Castellnou et al. 2022; Castellnou Ribau et al. 2024): ratios well above 1 mark
deep, hot, dry columns where a fire can decouple from the surface *even where the moist ladder
scores low*. Like the potential map it assumes a fire everywhere (a fixed reference flux), so
it is an upper bound and not a calibrated class -- no dry/moist LCL split is applied. Columns
the classifier rejects (ABL below {int(ABL_MIN_M)} m, or no usable profile) are left blank.

## Per-province pyroconvection potential metrics

Peak-of-day (09-18Z) statistics per province. *Pot peak* / *Gate peak* are the highest
class reached (0 surface plume -> 4 deep pyroCu/pyroCb); *%pyroCu* is the maximum daytime
areal coverage of class >= 2 (pyrocumulus), *%deep* of class 4. *Decoup max* is the peak
daytime fireABL/ABL ratio over the province -- read it alongside the classes, since a high
ratio with a low class is the dry-decoupling case the moist ladder does not score. LCL/ABL,
ABL depth and the minimum ABL-top RH are the noon column diagnostics over the province's
classified land. Classes: {classes}.

{pyro_md}

## Per-province Canadian FWI System (12Z)

Province-mean of the full FWI System at solar noon on each valid day: the moisture **codes**
FFMC (fine fuels), DMC (duff), DC (deep drought); the **indices** ISI (spread) and BUI
(buildup); and the final **FWI**. The slow DMC/DC drought memory is seeded from
{spin_method}, giving a run-day antecedent state of **FFMC {spin.ffmc:.1f} / DMC
{spin.dmc:.1f} / DC {spin.dc:.1f}**; that state is then carried forward through the three
forecast days on the run's own 12Z weather. Forecast rain is taken as 0 (dry-spell
assumption); a wetting event would lower FFMC/DMC and the indices. Danger classes follow the
EFFIS FWI scale.

{fwi_md}

Generated by scripts/forecast_3day_province_report.py. Pyroconvection: pyflam.pyroconvection_type
(hybrid). FWI: pyflam.fwi + pyflam.fire_weather_history (ERA5/SIR spin-up). CSV tables accompany
this report in the same folder.
""")
    try:
        # run from REPORT_OUT so the basename image references resolve to the copied PNGs
        subprocess.run(["pandoc", os.path.basename(md), "-o", os.path.basename(pdf),
                        "--pdf-engine=tectonic"], check=True, capture_output=True,
                       timeout=420, cwd=REPORT_OUT)
        return md, pdf
    except Exception as e:
        sys.stderr.write(f"PDF build skipped ({e}); MD + CSVs still written\n")
        return md, None


def main():
    if RUN != 0:
        sys.exit(f"run {RUN:02d}Z: this product starts at 00Z of the run day, so day 0 would be "
                 f"missing its first {RUN} h. Use the 00Z run.")
    os.makedirs(REPORT_OUT, exist_ok=True)
    run_base_pipeline()

    gdf = load_provinces()
    sample = os.path.join(_raster_dir(DAYS[0]), f"pyroconv_potential_{NOON:02d}Z.tif")
    lbl = province_label_grid(sample, gdf)

    print("[report] pyroconvection potential metrics ...")
    pyro_rows = pyroconv_metrics(gdf, lbl)
    print("[report] gridded FWI (ERA5 spin-up + 3-day carry) ...")
    grids, land, spin, spin_method = fwi_grids()
    fwi_rows = fwi_metrics(gdf, lbl, grids, land)

    write_csv(os.path.join(REPORT_OUT, f"pyroconv_metrics_{RUNDATE}.csv"), pyro_rows,
              ["province", "day", "pot_peak", "pot_pyroCu", "pot_deep", "gate_peak",
               "gate_pyroCu", "decoup_max", "lcl_abl", "abl_m", "rh_top_min"])
    write_csv(os.path.join(REPORT_OUT, f"fwi_provinces_{RUNDATE}.csv"), fwi_rows,
              ["province", "day", "ffmc", "dmc", "dc", "isi", "bui", "fwi", "danger"])

    # bring the three figure series into the separate folder so the report is self-contained
    png_pot = os.path.join(REPORT_OUT, f"forecast_3day_potential_{RUNDATE}.png")
    png_gate = os.path.join(REPORT_OUT, f"forecast_3day_gated_{RUNDATE}.png")
    png_decoup = os.path.join(REPORT_OUT, f"forecast_3day_decoupling_{RUNDATE}.png")
    shutil.copyfile(os.path.join(BASE_OUT, f"forecast_3day_potential_{RUNDATE}.png"), png_pot)
    shutil.copyfile(os.path.join(BASE_OUT, f"forecast_3day_gated_{RUNDATE}.png"), png_gate)
    shutil.copyfile(os.path.join(BASE_OUT, f"forecast_3day_decoupling_{RUNDATE}.png"), png_decoup)

    md, pdf = build_report(pyro_rows, fwi_rows, spin, spin_method, png_pot, png_gate, png_decoup)
    print(f"\nOK 3-day provincial report -> {REPORT_OUT}")
    for p in (md, pdf,
              os.path.join(REPORT_OUT, f"pyroconv_metrics_{RUNDATE}.csv"),
              os.path.join(REPORT_OUT, f"fwi_provinces_{RUNDATE}.csv")):
        if p:
            print(f"  {p}")


if __name__ == "__main__":
    main()
