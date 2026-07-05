# SPDX-License-Identifier: AGPL-3.0-or-later
"""Operational propagation run for the wildfire ignited ~13:00 local on 2026-07-03
at 43.93561 N, 11.09635 E (Calvana ridge, Prato, Tuscany).

Physics: ICON-2I 2.2 km surface wind/moisture + real vertical sounding -> Cruz-2005
crown fire on the GEDI-derived canopy LCP over the Tuscany SB40 fuel grid, coupled
to the fire's own buoyant plume (OpenFOAM RANS re-solved every 30 min) with the
pyroconvection loft scaling. Ember spotting from the crown fireline intensity.

Start 2026-07-03 14:00 CEST (12:00 UTC), 6 h, dt = 30 min. Outputs: arrival-time +
fire-type GeoTIFFs, 30-min isochrone GeoJSON, a map PNG, and a metrics report.

Fire-weather chain: a multi-week FWI spin-up (SIR Toscana observed rain + ERA5
noon T/RH/wind) builds the drought state; a fire-day gridded FWI seeds the dead
fuel moisture (FFMC) and cures the live herb (BUI). See pyflam.fwi /
pyflam.fire_weather_history / pyflam.sir_toscana.

Known caveats (marked ``TODO`` in the code, to confirm/solve):
  * SIR Toscana endpoint is a best-effort guess (only a 429 seen so far, never a
    verified 200) -- confirm it, or use pyflam.sir_toscana.read_sir_csv meanwhile.
  * ERA5 has ~5-day latency: days near the present repeat the last available step.
  * The default LCP is a user-specific absolute path -- override with PYFLAM_LCP.
"""
from __future__ import annotations
import datetime as dt
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.dirname(__file__))

import pyflam
from pyflam import units, landscape as L
from pyflam.io_lcp import read_lcp
from _icon2i_provider import build_icon2i, profile_at

# ----------------------------------------------------------------------------- config
LAT, LON = 43.93561, 11.09635
START = dt.datetime(2026, 7, 3, 12, tzinfo=dt.timezone.utc)     # 14:00 CEST
TOTAL = float(os.environ.get("PYFLAM_TOTAL", "360"))
DT = float(os.environ.get("PYFLAM_DT", "30"))
FOLIAR = 100.0                       # % live foliar moisture (summer broadleaf/conifer)
# Spin the time-lag dead-fuel model up over the hours before the burn window so the
# 10/100-h classes remember the humid overnight/morning air. 12 h reaches back to the
# 00Z ICON-2I run (the earliest valid time available from it).
FM_SPINUP_H = float(os.environ.get("PYFLAM_FM_SPINUP_H", "12"))
# Fire Weather Index spin-up: weeks of history before the fire (SIR Toscana observed
# rain + ERA5 T/RH/wind) build the drought state carried into the fire-day gridded
# FWI. Set PYFLAM_NO_SIR to use ERA5 precipitation only; ERA5 needs CDS credentials.
FWI_SPINUP_DAYS = int(os.environ.get("PYFLAM_FWI_DAYS", "28"))
USE_SIR = os.environ.get("PYFLAM_NO_SIR") is None
AOI_HALF_KM = float(os.environ.get("PYFLAM_HALF_KM", "18"))   # half-width box (km)
DOWNSAMPLE = int(os.environ.get("PYFLAM_DS", "2"))            # 100 m * DS -> grid res
# asymmetric margins around the ignition (km); default to the symmetric half-width.
MARGIN = {s: float(os.environ.get(f"PYFLAM_MARGIN_{s}", AOI_HALF_KM))
          for s in ("W", "E", "N", "S")}
# TODO(solve): the default is a user-specific absolute path. Ship a small sample
# .lcp (or require PYFLAM_LCP) so the pipeline runs unmodified on other machines.
LCP = os.environ.get(
    "PYFLAM_LCP",
    "/Users/cristianofoderi/DATI/FUEL_TOS/pyflam_canopy_tuscany/canopy_tuscany.lcp")
OUT = os.path.join(os.path.dirname(__file__), "..", "docs", "fire_calvana_2026-07-03")
CFD_KW = dict(iterations=int(os.environ.get("PYFLAM_ITERS", "350")),
              nz=int(os.environ.get("PYFLAM_NZ", "12")), n_processors=6)
os.makedirs(OUT, exist_ok=True)

FT_M = 0.3048
def _write_raster(ls, data, path, *, dtype="float32", nodata=-9999.0):
    """Write a landscape-aligned GeoTIFF: pyflam_gui if present, else rasterio."""
    try:
        from pyflam_gui.core import outputs as gout
        gout.write_geotiff(ls, np.asarray(data, dtype), path, dtype=dtype, nodata=nodata)
        return
    except Exception:
        pass
    import rasterio
    from rasterio.transform import from_origin
    a = np.asarray(data, dtype=dtype)
    with rasterio.open(path, "w", driver="GTiff", height=a.shape[0], width=a.shape[1],
                       count=1, dtype=dtype, crs=ls.crs or "EPSG:3035",
                       transform=from_origin(ls.west, ls.north, ls.cellsize_x, ls.cellsize_y),
                       nodata=nodata, compress="deflate") as dst:
        dst.write(a, 1)


def r_mmin(a): return np.asarray(a) * FT_M                        # ft/min -> m/min
def fli_kwm(a): return units.btu_per_ft_s_to_kw_per_m(np.asarray(a))
def fl_m(a): return np.asarray(a) * FT_M


def log(*a):
    print(*a, flush=True)


# ----------------------------------------------------------------------------- landscape
def clip_landscape():
    from pyproj import Transformer
    F = DOWNSAMPLE
    d = read_lcp(LCP)
    tx = Transformer.from_crs("EPSG:4326", "EPSG:3035", always_xy=True)
    x, y = tx.transform(LON, LAT)
    c0 = int((x - MARGIN["W"] * 1000 - d.west) / d.cellsize_x)
    c1 = int((x + MARGIN["E"] * 1000 - d.west) / d.cellsize_x)
    r0 = int((d.north - (y + MARGIN["N"] * 1000)) / d.cellsize_y)
    r1 = int((d.north - (y - MARGIN["S"] * 1000)) / d.cellsize_y)
    c0, r0 = max(0, c0), max(0, r0); c1, r1 = min(d.ncols, c1), min(d.nrows, r1)
    # trim to a multiple of F so the block-downsample is exact
    c1 -= (c1 - c0) % F; r1 -= (r1 - r0) % F
    b = {k: np.asarray(v)[r0:r1, c0:c1] for k, v in d.bands.items()}

    def mean_ds(a):
        a = a.astype(float)
        return a.reshape(a.shape[0] // F, F, a.shape[1] // F, F).mean(axis=(1, 3))

    def near_ds(a):                                   # class-preserving (fuel/aspect)
        return a[::F, ::F]

    csx, csy = d.cellsize_x * F, d.cellsize_y * F
    west = d.west + c0 * d.cellsize_x
    north = d.north - r0 * d.cellsize_y
    ls = L.Landscape(
        fuel_model=near_ds(b["fuel_model"]), slope=mean_ds(b["slope"]),
        elevation=mean_ds(b["elevation"]) if b.get("elevation") is not None else None,
        aspect=near_ds(b["aspect"]) if b.get("aspect") is not None else None,
        canopy_cover=mean_ds(b["canopy_cover"]) if b.get("canopy_cover") is not None else None,
        canopy_height=mean_ds(b["canopy_height"]) if b.get("canopy_height") is not None else None,
        canopy_base_height=mean_ds(b["canopy_base_height"]) if b.get("canopy_base_height") is not None else None,
        canopy_bulk_density=mean_ds(b["canopy_bulk_density"]) if b.get("canopy_bulk_density") is not None else None,
        cellsize_x=csx, cellsize_y=csy, west=west, north=north,
        slope_units=d.slope_units, crs="EPSG:3035")
    irow = int((north - y) / csy); icol = int((x - west) / csx)
    irow = min(max(irow, 0), ls.shape[0] - 1); icol = min(max(icol, 0), ls.shape[1] - 1)
    return ls, (irow, icol), (x, y)


# ----------------------------------------------------------------------------- main
def main():
    log("== clip landscape ==")
    ls, ign, xy = clip_landscape()
    log(f"AOI shape {ls.shape}  ignition cell {ign}  fuel {int(np.asarray(ls.fuel_model)[ign])}")

    log("== build ICON-2I provider + sounding ==")
    prov, ds, levels, files = build_icon2i(START.date(), run=0)
    st0 = prov.state_at(LAT, LON, START)
    prof = profile_at(prov, levels, st0, LAT, LON, START)
    log(f"start state: wind {st0.wind_speed:.1f} m/s @ {st0.wind_direction:.0f}  "
        f"T {st0.temperature:.1f}C  RH {st0.relative_humidity:.0f}%")

    # --- Fire Weather Index: a multi-week spin-up (SIR Toscana observed rain + ERA5
    # noon T/RH/wind) builds the drought state, then a fire-day gridded FWI over the
    # AOI. It feeds the spread two ways: the FFMC fine fuel moisture seeds the dead-
    # fuel model, and the drought (BUI) cures the live herb. Degrades gracefully
    # (spring-startup FWI) when ERA5/SIR are unavailable, so the run always proceeds.
    from pyflam.fuel_conditioning import (live_fuel_moisture, live_fuel_aridity,
                                          drought_curing)
    ffmc_fm = None; fwi_grid = None; bui_mean = 0.0
    try:
        from pyflam import fire_weather_history as fwh, fwi as _fwi
        fwi_prev = fwh.fire_weather_spinup(START.date(), LAT, LON,
                                           days=FWI_SPINUP_DAYS, use_sir=USE_SIR, log=log)
        fwi_grid = _fwi.gridded_fwi_on(prov.field_on(ls, START),
                                       month=START.month, previous=fwi_prev)
        ffmc_fm = float(np.nanmean(np.asarray(fwi_grid["fine_fuel_moisture"]))) / 100.0
        bui_mean = float(np.nanmean(np.asarray(fwi_grid["bui"])))
        log(f"FWI fire-day (AOI mean): FFMC fine-fuel {ffmc_fm*100:.1f}%  "
            f"BUI {bui_mean:.0f}  FWI {float(np.nanmean(np.asarray(fwi_grid['fwi']))):.1f}")
        try:                                  # gridded FWI GeoTIFFs over the AOI
            for k in ("fwi", "bui", "isi", "ffmc"):
                _write_raster(ls, np.asarray(fwi_grid[k], float),
                              os.path.join(OUT, f"fwi_{k}.tif"))
            log("wrote gridded FWI GeoTIFFs (fwi/bui/isi/ffmc)")
        except Exception as e:
            log("FWI geotiff export skipped:", e)
    except Exception as e:
        log("FWI spin-up/grid skipped:", e)

    # Live herb/woody moisture: a seasonal (phenological) estimate for the ignition
    # date -- it cannot come from a single ICON snapshot -- cured by whichever is drier
    # of the burn-window VPD and the multi-week FWI drought (BUI).
    arid = max(live_fuel_aridity(st0.temperature, st0.relative_humidity),
               drought_curing(bui_mean))
    live = live_fuel_moisture(START.timetuple().tm_yday, aridity=arid)
    log(f"live fuel moisture (est.): herb {live['m_live_herb']*100:.0f}%  "
        f"woody {live['m_live_woody']*100:.0f}%  (aridity {arid:.2f})")

    march_kw = dict(wind_relax=0.5, max_wind_factor=3.0, return_history=True)
    if os.environ.get("PYFLAM_NOCFD"):
        log("== march: NO-CFD validation mode (uniform ambient wind provider) ==")
        from pyflam.pyroconvection import _uniform_wind_field
        march_kw["wind_provider"] = lambda ls_, inten, active, spd, dirn: _uniform_wind_field(ls_, spd, dirn)
    else:
        log("== march: ICON-2I + Cruz crown + CFD plume + pyroconvection (return_history) ==")
        march_kw.update(CFD_KW)
    out = pyflam.fire_atmosphere_march(
        ls, [ign], total_time=TOTAL, dt=DT, atmosphere=prov, location=(LAT, LON),
        start_time=START, crown=True, foliar_moisture=FOLIAR, crown_spread="cruz2005",
        fuel_moisture_lag=True, fuel_moisture_spinup_hours=FM_SPINUP_H,
        fuel_moisture_lag_init=ffmc_fm,       # seed dead fuel from the FWI-spun FFMC
        m_live_herb=live["m_live_herb"], m_live_woody=live["m_live_woody"],
        pyroconvection=True, profile=prof, **march_kw)
    arrival = np.asarray(out["arrival_time"], dtype=float)
    fields = out["fields"]; winds = out["winds"]; times = out["times"]
    pf_hist = out.get("plume_factor", [])
    m1h_hist = out.get("m_1h", [])            # lagged 1-h dead fuel moisture per step
    log(f"march done. burned cells {int(np.isfinite(arrival).sum())}  "
        f"max arrival {np.nanmax(np.where(np.isfinite(arrival), arrival, np.nan)):.0f} min")

    # -------- checkpoint the expensive march output before post-processing --------
    def _stack(getter):
        return np.array([np.asarray(getter(f)) for f in fields])
    try:
        np.savez_compressed(
            os.path.join(OUT, "march_checkpoint.npz"),
            arrival=arrival, times=np.array(times), pf_hist=np.array(pf_hist),
            ros_max=_stack(lambda f: f.ros_max),
            fli=_stack(lambda f: f.fireline_intensity),
            heading=_stack(lambda f: f.heading),
            wspeed=np.array([np.asarray(w.speed_ft_per_min()) for w in winds]),
            wdir=np.array([np.asarray(w.direction) for w in winds]),
            fire_type=np.asarray(out.get("fire_type")) if out.get("fire_type") is not None else np.zeros(0),
            crown_fraction=np.asarray(out.get("crown_fraction_burned")) if out.get("crown_fraction_burned") is not None else np.zeros(0),
            ignition=np.array(ign), cellsize=np.array([ls.cellsize_x, ls.cellsize_y]))
        log("== checkpoint saved: march_checkpoint.npz")
    except Exception as e:
        log("checkpoint save failed:", type(e).__name__, e)

    # -------- per-30-min metrics --------
    from pyflam import FirebrandPhysics
    rng = np.random.default_rng(42)
    cell_ha = ls.cellsize_x * ls.cellsize_y / 1e4
    fb = FirebrandPhysics()
    rows = []
    prev_spots = 0
    for i, t in enumerate(times):
        fld = fields[min(i, len(fields) - 1)]
        wf = winds[min(i, len(winds) - 1)]
        ann = np.isfinite(arrival) & (arrival > t - DT) & (arrival <= t)     # burned this step
        burned_cum = np.isfinite(arrival) & (arrival <= t)
        ros = r_mmin(np.asarray(fld.ros_max))                 # heading ROS, ft/min->m/min
        fli = fli_kwm(np.asarray(fld.fireline_intensity))     # Btu/ft/s -> kW/m
        fl = 0.0775 * np.power(np.clip(fli, 0.0, None), 0.46)  # Byram flame length (m)
        state = prov.state_at(LAT, LON, START + dt.timedelta(minutes=t))
        # lagged 1-h moisture the march actually used (falls back to instantaneous)
        m1h = (float(m1h_hist[min(i, len(m1h_hist) - 1)]) if m1h_hist
               else pyflam.atmosphere.dead_fuel_moisture(state)["m_1h"])
        rec = {
            "t_min": int(t),
            "clock_local": (START + dt.timedelta(minutes=t) + dt.timedelta(hours=2)).strftime("%H:%M"),
            "burned_ha": round(float(burned_cum.sum()) * cell_ha, 1),
            "step_ha": round(float(ann.sum()) * cell_ha, 1),
            "ros_head_p95_mmin": round(float(np.nanpercentile(ros[ann], 95)) if ann.any() else 0.0, 1),
            "ros_mean_mmin": round(float(np.nanmean(ros[ann])) if ann.any() else 0.0, 1),
            "fli_max_kwm": round(float(np.nanmax(fli[ann])) if ann.any() else 0.0, 0),
            "fli_mean_kwm": round(float(np.nanmean(fli[ann])) if ann.any() else 0.0, 0),
            "flame_max_m": round(float(np.nanmax(fl[ann])) if ann.any() else 0.0, 1),
            "wind_ms": round(float(state.wind_speed), 1),
            "wind_dir": round(float(state.wind_direction), 0),
            "rh_pct": round(float(state.relative_humidity), 0),
            "m1h_pct": round(float(m1h) * 100, 1),
            "plume_factor": round(float(pf_hist[i]) if i < len(pf_hist) else 1.0, 2),
        }
        # spotting (cumulative to t) from crown fireline intensity
        try:
            w20 = float(np.nanmean(np.asarray(wf.speed_ft_per_min())))   # scalar 20-ft wind
            spots = fb.generate_spots(fld, arrival, wind_20ft=w20,
                                      wind_direction=float(state.wind_direction),
                                      max_time=t, rng=rng, fuel_moisture=float(m1h))
            n_new = max(0, len(spots) - prev_spots)
            maxd = 0.0
            if spots:
                dr = np.array([(sr - ign[0]) * ls.cellsize_y for sr, sc, _ in spots])
                dc = np.array([(sc - ign[1]) * ls.cellsize_x for sr, sc, _ in spots])
                maxd = float(np.hypot(dr, dc).max())
            rec["spots_cum"] = len(spots)
            rec["spot_maxdist_m"] = round(maxd, 0)
            prev_spots = len(spots)
        except Exception as e:
            rec["spots_cum"] = None; rec["spot_maxdist_m"] = None
            log("  spotting failed:", type(e).__name__, e)
        rows.append(rec)
        log(f"  t={int(t):3d} {rec['clock_local']}  area {rec['burned_ha']:6.1f} ha  "
            f"ROS_head {rec['ros_head_p95_mmin']:5.1f} m/min  FLI_max {rec['fli_max_kwm']:6.0f} kW/m  "
            f"flame {rec['flame_max_m']:4.1f} m  spots {rec['spots_cum']}")

    # -------- fire type & crown fraction --------
    ft = out.get("fire_type"); cfb = out.get("crown_fraction_burned")
    ft_summary = {}
    burned = np.isfinite(arrival)
    if ft is not None:
        fta = np.asarray(ft)
        if fta.dtype.kind in "iu" or fta.dtype.kind == "f":
            uniq, cnt = np.unique(fta[burned], return_counts=True)
            ft_summary = {str(int(u)): int(c) for u, c in zip(uniq, cnt)}
        else:
            uniq, cnt = np.unique(fta[burned].astype(str), return_counts=True)
            ft_summary = {str(u): int(c) for u, c in zip(uniq, cnt)}
    crown_ha = None
    if cfb is not None:
        crown_ha = round(float((np.asarray(cfb)[burned] > 0.1).sum()) * cell_ha, 1)

    # -------- convective diagnostics --------
    from pyflam import continuous_haines, inverted_v, lcl_height_m, pyroconvection_potential
    ch = float(continuous_haines(prof))
    iv, depr, midrh = inverted_v(prof)
    lcl = float(lcl_height_m(st0.temperature, st0.relative_humidity))
    conv = {
        "continuous_haines": round(ch, 1), "inverted_v": bool(iv),
        "sfc_dewpoint_depression_C": round(float(depr), 1), "mid_rh_pct": round(float(midrh), 0),
        "lcl_m": round(lcl, 0),
        "plume_factor_max": round(float(np.nanmax(pf_hist)) if len(pf_hist) else 1.0, 2),
    }
    pot = out.get("pyroconvection")
    if pot is not None:
        conv["pyroconvection_type"] = getattr(pot, "type", getattr(pot, "label", str(pot)))
    if "pyrocb_firepower_threshold" in out:
        conv["pyrocb_firepower_threshold_GW"] = round(float(out["pyrocb_firepower_threshold"]) / 1e9, 2)

    # -------- operative sector analysis at final time --------
    op_summary = None
    try:
        prod_state = prov.state_at(LAT, LON, START + dt.timedelta(minutes=TOTAL))
        si = pyflam.spread_inputs_from_state(prod_state)
        op = pyflam.analyze_perimeter(
            ls, arrival, TOTAL, wind_midflame=si["wind_midflame"],
            wind_direction=si["wind_direction"], m_1h=si["m_1h"], m_10h=si["m_10h"],
            m_100h=si["m_100h"], load_factor=1.0)
        op_summary = op.summary()
    except Exception as e:
        log("operative analysis failed:", type(e).__name__, e)

    # -------- write rasters + isochrones + map --------
    write_outputs(ls, arrival, ft, cfb, times)

    result = {
        "meta": {"lat": LAT, "lon": LON, "start_utc": START.isoformat(),
                 "start_local": "2026-07-03 14:00 CEST", "total_min": TOTAL, "dt_min": DT,
                 "aoi_shape": list(ls.shape), "ignition_cell": list(ign),
                 "fuel_model_at_ignition": int(np.asarray(ls.fuel_model)[ign]),
                 "physics": "ICON-2I + Cruz2005 crown + OpenFOAM RANS plume + pyroconvection + spotting"},
        "intervals": rows, "fire_type_cells": ft_summary, "crown_area_ha": crown_ha,
        "convection": conv, "operative": op_summary,
        "total_burned_ha": round(float(burned.sum()) * cell_ha, 1),
    }
    with open(os.path.join(OUT, "metrics.json"), "w") as fh:
        json.dump(result, fh, indent=2)
    log("== wrote", os.path.join(OUT, "metrics.json"))
    write_report(result)
    log("== DONE ==")


def write_outputs(ls, arrival, ft, cfb, times):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # GeoTIFFs
    try:
        from pyflam_gui.core import outputs as gout
    except Exception:
        gout = None
    arr = np.where(np.isfinite(arrival), arrival, np.nan)
    if gout is not None:
        try:
            gout.write_geotiff(ls, arr, os.path.join(OUT, "arrival_time.tif"))
            if cfb is not None:
                gout.write_geotiff(ls, np.asarray(cfb, float), os.path.join(OUT, "crown_fraction.tif"),
                                   dtype="float32", nodata=0)
        except Exception as e:
            log("geotiff failed:", e)
    # isochrone GeoJSON via operative.perimeter_rings at each 30-min
    try:
        from pyproj import Transformer
        from pyflam.operative import perimeter_rings
        tx = Transformer.from_crs("EPSG:3035", "EPSG:4326", always_xy=True).transform
        feats = []
        for t in times:
            mask = np.isfinite(arrival) & (arrival <= t)
            rings = perimeter_rings(mask, ls, transform=tx) or []
            for ring in rings:
                feats.append({"type": "Feature",
                              "geometry": {"type": "LineString", "coordinates": ring},
                              "properties": {"t_min": int(t)}})
        with open(os.path.join(OUT, "isochrones_30min.geojson"), "w") as fh:
            json.dump({"type": "FeatureCollection", "features": feats}, fh)
    except Exception as e:
        log("isochrone geojson failed:", e)
    # map PNG: isochrones over hillshade
    try:
        fig, ax = plt.subplots(figsize=(8, 8))
        elev = np.asarray(ls.elevation, float)
        from matplotlib.colors import LightSource
        ls_sh = LightSource(azdeg=315, altdeg=45)
        ax.imshow(ls_sh.hillshade(elev, vert_exag=2, dx=ls.cellsize_x, dy=ls.cellsize_y),
                  cmap="gray", alpha=0.7)
        am = np.where(np.isfinite(arrival), arrival, np.nan)
        cs = ax.contour(am, levels=list(times), cmap="inferno", linewidths=1.2)
        ax.clabel(cs, inline=True, fontsize=7, fmt=lambda v: f"{int(v)}m")
        ir, ic = np.unravel_index(np.nanargmin(np.where(np.isfinite(arrival), arrival, np.inf)), arrival.shape)
        ax.plot(ic, ir, "*", color="red", ms=16, mec="white")
        ax.set_title("Calvana fire 2026-07-03 — 30-min isochrones (ICON-2I + CFD plume + crown)")
        ax.set_xticks([]); ax.set_yticks([])
        fig.tight_layout(); fig.savefig(os.path.join(OUT, "isochrones_map.png"), dpi=140)
        log("wrote isochrones_map.png")
    except Exception as e:
        log("map png failed:", e)


def write_report(r):
    m = r["meta"]
    lines = [f"# Calvana wildfire — pyflam propagation run",
             f"", f"**Ignition** {m['lat']}, {m['lon']} · **start** {m['start_local']} "
             f"({m['start_utc']}) · **duration** {int(m['total_min'])} min · dt {int(m['dt_min'])} min",
             f"**Physics** {m['physics']}",
             f"**AOI** {m['aoi_shape'][0]}×{m['aoi_shape'][1]} @100 m · fuel model at ignition "
             f"SB40 #{m['fuel_model_at_ignition']}", f"",
             f"**Total burned area:** {r['total_burned_ha']} ha · **crown-involved:** {r['crown_area_ha']} ha",
             f"", f"## Time-of-arrival / behaviour every 30 min", f"",
             f"| t (min) | local | area ha | ROS head p95 (m/min) | FLI max (kW/m) | flame (m) | wind m/s | RH % | 1-h % | plume× | spots |",
             f"|--:|:--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|"]
    for x in r["intervals"]:
        lines.append(f"| {x['t_min']} | {x['clock_local']} | {x['burned_ha']} | {x['ros_head_p95_mmin']} | "
                     f"{x['fli_max_kwm']:.0f} | {x['flame_max_m']} | {x['wind_ms']} | {x['rh_pct']:.0f} | "
                     f"{x['m1h_pct']} | {x['plume_factor']} | {x['spots_cum']} |")
    c = r["convection"]
    lines += ["", "## Convective / pyro behaviour", "",
              f"- Continuous-Haines **{c['continuous_haines']}**, inverted-V {c['inverted_v']}, "
              f"surface dewpoint depression {c['sfc_dewpoint_depression_C']} °C, mid-level RH "
              f"{c['mid_rh_pct']:.0f}%, LCL {c['lcl_m']:.0f} m",
              f"- Plume loft factor (max) **{c['plume_factor_max']}**",
              f"- Fire-type cell counts: {r['fire_type_cells']}"]
    if r.get("operative"):
        lines += ["", "## Operative sector analysis (final)", "", "```", r["operative"], "```"]
    with open(os.path.join(OUT, "report.md"), "w") as fh:
        fh.write("\n".join(lines))
    log("wrote report.md")


if __name__ == "__main__":
    main()
