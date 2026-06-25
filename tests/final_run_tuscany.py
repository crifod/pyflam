#!/usr/bin/env python
"""End-to-end pyflam run on the real Tuscany landscape -- a complete-function test.

Scenario: 300 random ignitions across the region, 28 June 2026, fire running for
36 hours. Drives the whole pipeline and writes each phase's output into its own
subfolder:

    00_landscape/        band rasters + summary of the loaded landscape
    01_weather/          atmospheric state for the scenario (live GFS, else manual)
    02_fuel_moisture/    per-cell conditioned dead fuel moisture (terrain + canopy)
    03_surface_behavior/ Rothermel surface rasters (ROS, flame length, intensity)
    04_wind/             terrain wind (mass-consistent solver) midflame field
    05_spread_field/     directional elliptical spread template (ros_max/heading/ecc)
    06_ignitions/        the 300 ignition points (CSV, row/col + world x/y + lat/lon)
    07_burn_probability/ BP + connected metrics over 36 h (the FlamMap randig analog)

Usage (from the pyflam dir):
    PYTHONPATH=src python tests/final_run_tuscany.py
    PYTHONPATH=src python tests/final_run_tuscany.py --quick      # fast smoke run
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime

import numpy as np

import pyflam
from pyflam import fuel_models
from pyflam.units import mph_to_ft_per_min, ft_per_min_to_m_per_min

DATA = "/Users/cristianofoderi/DATI/FUEL_TOS/Toscana_Fuel_data"
GEOTIFFS = {
    "fuel_model": f"{DATA}/Fuel_Model_40_tos_epsg3035.tif",
    "slope": f"{DATA}/slope_deg_3035.tif",
    "aspect": f"{DATA}/aspect_deg_azimut_3035.tif",
    "elevation": f"{DATA}/dem_3035_sm.tif",
    "canopy_cover": f"{DATA}/cancov.tif",
}
LIVE_MOIST = dict(m_live_herb=0.60, m_live_woody=0.90)


def _log(fh, msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    fh.write(line + "\n")
    fh.flush()


def _step(out, name):
    d = os.path.join(out, name)
    os.makedirs(d, exist_ok=True)
    return d


def _burnable(ls):
    burn = np.zeros(ls.shape, dtype=bool)
    for num in np.unique(ls.fuel_model):
        try:
            if fuel_models.get(int(num)).is_burnable:
                burn |= ls.fuel_model == int(num)
        except KeyError:
            pass
    return burn


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="/Users/cristianofoderi/DATI/FUEL_TOS/"
                    "pyflam_final_run_2026-06-28")
    ap.add_argument("--nfires", type=int, default=300)
    ap.add_argument("--hours", type=float, default=36.0)
    ap.add_argument("--datetime", default="2026-06-28 12:00", help="scenario time (UTC)")
    ap.add_argument("--gfs-run", default="2026-06-25 00:00", help="GFS cycle to forecast from")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--quick", action="store_true", help="tiny smoke run (no GFS)")
    ap.add_argument("--terrain-wind", action="store_true",
                    help="run the mass-consistent wind solver (heavy on 5M cells)")
    args = ap.parse_args(argv)

    if args.quick:
        args.nfires, args.hours = 8, 2.0
    os.makedirs(args.out, exist_ok=True)
    fh = open(os.path.join(args.out, "run_log.txt"), "w")
    summary = {"scenario": args.datetime, "nfires": args.nfires, "hours": args.hours,
               "seed": args.seed}
    rng = np.random.default_rng(args.seed)
    when = datetime.strptime(args.datetime, "%Y-%m-%d %H:%M")
    max_time = args.hours * 60.0
    t0 = time.time()
    _log(fh, f"pyflam complete-function run | {args.nfires} fires | {args.hours} h | "
             f"{args.datetime}")

    # --- 00 landscape ---------------------------------------------------------
    d = _step(args.out, "00_landscape")
    ls = pyflam.Landscape.from_geotiffs(GEOTIFFS)
    # Real-world GeoTIFFs carry nodata NaNs (edges / water); fill them so the
    # per-cell terrain math (insolation, slope) stays finite. Nonburnable anyway.
    for band in ("slope", "aspect", "elevation", "canopy_cover"):
        a = getattr(ls, band, None)
        if a is not None:
            setattr(ls, band, np.nan_to_num(np.asarray(a, dtype=float), nan=0.0))
    burnable = _burnable(ls)
    for band in ("fuel_model", "slope", "aspect", "elevation", "canopy_cover"):
        ls.to_geotiff(os.path.join(d, f"{band}.tif"),
                      np.asarray(getattr(ls, band), dtype=float))
    info = {"shape": list(ls.shape), "cellsize_m": round(ls.cellsize_x, 2),
            "crs": str(ls.crs), "burnable_cells": int(burnable.sum()),
            "burnable_pct": round(100 * burnable.mean(), 1),
            "n_fuel_models": int(np.unique(ls.fuel_model).size)}
    json.dump(info, open(os.path.join(d, "landscape.json"), "w"), indent=2)
    _log(fh, f"00 landscape {ls.shape} EPSG:3035; burnable {info['burnable_pct']}%")
    ll = pyflam.atmosphere.latlon_grid(ls)
    clat = float(np.nanmean(ll[0])); clon = float(np.nanmean(ll[1]))

    # --- 01 weather -----------------------------------------------------------
    d = _step(args.out, "01_weather")
    prov, source = None, "manual"
    if not args.quick:
        try:
            run = datetime.strptime(args.gfs_run, "%Y-%m-%d %H:%M")
            fxx = int(round((when - run).total_seconds() / 3600.0))
            _log(fh, f"01 fetching GFS {args.gfs_run} +{fxx}h -> {args.datetime} ...")
            prov = pyflam.atmosphere.fetch_gfs(run=run, fxx=fxx)
            source = f"GFS {args.gfs_run} +{fxx}h"
        except Exception as exc:
            _log(fh, f"   GFS fetch failed ({exc}); falling back to manual weather")
    if prov is not None:
        st = prov.state_at(clat, clon, when)
        temp_c, rh = st.temperature, st.relative_humidity
        wind_mph = st.wind_speed * 2.23694
        wind_dir = st.wind_direction
    else:
        temp_c, rh, wind_mph, wind_dir = 30.0, 25.0, 16.0, 225.0   # hot dry SW wind
    wx = {"source": source, "temperature_C": round(float(temp_c), 1),
          "relative_humidity_pct": round(float(rh), 1),
          "wind_mph": round(float(wind_mph), 1), "wind_dir_deg": round(float(wind_dir), 1),
          "site_lat": round(clat, 3), "site_lon": round(clon, 3)}
    json.dump(wx, open(os.path.join(d, "weather.json"), "w"), indent=2)
    summary["weather"] = wx
    _log(fh, f"01 weather [{source}]: T={wx['temperature_C']}C RH={wx['relative_humidity_pct']}%"
             f" wind {wx['wind_mph']}mph @ {wx['wind_dir_deg']}deg")

    # --- 02 fuel moisture (conditioned per cell) ------------------------------
    d = _step(args.out, "02_fuel_moisture")
    moist = pyflam.condition_from_weather(
        ls, time=when, atmosphere=prov, temperature=temp_c, relative_humidity=rh,
        latitude=clat, longitude=clon, model="emc")
    for k in ("m_1h", "m_10h", "m_100h"):
        ls.to_geotiff(os.path.join(d, f"{k}.tif"), 100.0 * moist[k])   # percent
    m1 = 100.0 * moist["m_1h"][burnable]
    fm_sum = {"m_1h_pct": [round(float(np.nanmin(m1)), 1), round(float(np.nanmean(m1)), 1),
                           round(float(np.nanmax(m1)), 1)]}
    json.dump(fm_sum, open(os.path.join(d, "moisture.json"), "w"), indent=2)
    _log(fh, f"02 conditioned 1-h moisture (burnable): {fm_sum['m_1h_pct']} % (min/mean/max)")

    midflame = pyflam.midflame_field(ls, mph_to_ft_per_min(wind_mph))

    # --- 03 surface behavior --------------------------------------------------
    d = _step(args.out, "03_surface_behavior")
    sb = pyflam.basic_fire_behavior(ls, wind_midflame=midflame, **moist, **LIVE_MOIST)
    for k in ("rate_of_spread", "flame_length", "fireline_intensity", "reaction_intensity"):
        ls.to_geotiff(os.path.join(d, f"{k}.tif"), np.asarray(sb[k], dtype=float))
    ros = np.asarray(sb["rate_of_spread"])[burnable]
    _log(fh, f"03 surface ROS (burnable) mean {ft_per_min_to_m_per_min(ros.mean()):.1f} "
             f"m/min, max {ft_per_min_to_m_per_min(ros.max()):.1f} m/min")

    # --- 04 terrain wind (mass-consistent) ------------------------------------
    d = _step(args.out, "04_wind")
    wfl_dir = wind_dir
    if args.terrain_wind:
        try:
            tw = time.time()
            wf = pyflam.wind_field_from_landscape(
                ls, speed=mph_to_ft_per_min(wind_mph), direction=wind_dir)
            ls.to_geotiff(os.path.join(d, "wind_20ft_ftmin.tif"),
                          np.asarray(wf.speed_ft_per_min(), dtype=float))
            midflame = wf.speed_ft_per_min() * 0.4
            wfl_dir = wf.direction
            _log(fh, f"04 terrain wind: mass-consistent solver ok ({time.time()-tw:.1f}s)")
        except Exception as exc:
            _log(fh, f"04 terrain wind solver failed ({exc}); uniform wind used")
            json.dump({"note": f"solver failed: {exc}"},
                      open(os.path.join(d, "wind.json"), "w"), indent=2)
    else:
        ls.to_geotiff(os.path.join(d, "midflame_ftmin.tif"),
                      np.broadcast_to(np.asarray(midflame, dtype=float),
                                      ls.shape).copy())
        json.dump({"note": "ambient wind -> canopy/fuel wind-reduced midflame "
                   "(pass --terrain-wind for the mass-consistent terrain solver)",
                   "wind_mph": round(float(wind_mph), 1)},
                  open(os.path.join(d, "wind.json"), "w"), indent=2)
        _log(fh, "04 wind: wind-reduced midflame (terrain solver off; --terrain-wind to enable)")

    # --- 05 spread field ------------------------------------------------------
    d = _step(args.out, "05_spread_field")
    field = pyflam.spread_field(ls, wind_midflame=midflame, wind_direction=wfl_dir,
                                **moist, **LIVE_MOIST)
    ls.to_geotiff(os.path.join(d, "ros_max_ftmin.tif"), np.asarray(field.ros_max))
    ls.to_geotiff(os.path.join(d, "heading_deg.tif"), np.asarray(field.heading))
    ls.to_geotiff(os.path.join(d, "eccentricity.tif"), np.asarray(field.eccentricity))
    _log(fh, "05 directional spread field built")

    # --- 06 ignitions ---------------------------------------------------------
    d = _step(args.out, "06_ignitions")
    rows, cols = np.where(burnable)
    pick = rng.choice(rows.size, size=min(args.nfires, rows.size), replace=False)
    ign_rc = list(zip(rows[pick].tolist(), cols[pick].tolist()))
    lat2d, lon2d = ll
    with open(os.path.join(d, "ignitions.csv"), "w") as f:
        f.write("idx,row,col,x,y,lat,lon\n")
        for i, (r, c) in enumerate(ign_rc):
            x = ls.west + (c + 0.5) * ls.cellsize_x
            y = ls.north - (r + 0.5) * ls.cellsize_y
            f.write(f"{i},{r},{c},{x:.1f},{y:.1f},{lat2d[r,c]:.5f},{lon2d[r,c]:.5f}\n")
    _log(fh, f"06 {len(ign_rc)} random ignitions written")

    # --- 07 burn probability + connected metrics (36 h) -----------------------
    d = _step(args.out, "07_burn_probability")
    _log(fh, f"07 building travel-time graph + running {len(ign_rc)} fires of "
             f"{args.hours} h ...")
    tg = time.time()
    graph = pyflam.mtt.build_traveltime_graph(field)
    _log(fh, f"   graph built in {time.time()-tg:.1f}s")
    tb = time.time()
    res = pyflam.burn_probability(field, ign_rc, max_time=max_time, graph=graph,
                                  return_metrics=True, rng=rng)
    _log(fh, f"   {res.n_fires} fires done in {time.time()-tb:.1f}s")
    ls.to_geotiff(os.path.join(d, "burn_probability.tif"), res.burn_prob)
    ls.to_geotiff(os.path.join(d, "conditional_flame_length_ft.tif"),
                  np.nan_to_num(res.conditional_flame_length))
    ls.to_geotiff(os.path.join(d, "conditional_intensity_btu_ft_s.tif"),
                  np.nan_to_num(res.conditional_intensity))
    np.savetxt(os.path.join(d, "fire_sizes_ha.csv"), res.fire_sizes / 1e4,
               header="burned_area_ha", comments="")
    burned = res.burn_prob > 0
    sizes_ha = res.fire_sizes / 1e4
    bp_sum = {
        "n_fires": res.n_fires,
        "mean_burn_prob": round(float(res.burn_prob[burnable].mean()), 5),
        "cells_ever_burned_pct": round(100 * float((burned & burnable).mean()
                                                    / max(burnable.mean(), 1e-9)), 2),
        "fire_size_ha": {"mean": round(float(sizes_ha.mean()), 1),
                         "median": round(float(np.median(sizes_ha)), 1),
                         "max": round(float(sizes_ha.max()), 1)},
        "conditional_flame_length_ft_mean":
            round(float(np.nanmean(res.conditional_flame_length[burned])), 2),
    }
    json.dump(bp_sum, open(os.path.join(d, "burn_probability.json"), "w"), indent=2)
    summary["burn_probability"] = bp_sum
    _log(fh, f"07 mean BP {bp_sum['mean_burn_prob']}; fire size mean "
             f"{bp_sum['fire_size_ha']['mean']} ha, max {bp_sum['fire_size_ha']['max']} ha")

    summary["runtime_s"] = round(time.time() - t0, 1)
    json.dump(summary, open(os.path.join(args.out, "SUMMARY.json"), "w"), indent=2)
    _log(fh, f"DONE in {summary['runtime_s']}s -> {args.out}")
    fh.close()


if __name__ == "__main__":
    main()
