#!/usr/bin/env python
"""Validate pyflam's spread/perimeter (MTT) engine against a FlamMap MTT run.

The Tuscany dataset includes a FlamMap "MTT Random Ignitions" burn-probability
run (`burn_prob_test_*`). Its `run_log.txt` pins every input: 2000 random-ignition
fires, 60-min sim time, wind 30 km/h @ 10 m from 45 deg, moisture 6/7/8/60/90.
Two of its outputs exercise the spread engine:

* ``MAX_SPRE_DIR.tif`` (radians) — the per-cell maximum-spread direction. This is
  deterministic (wind + slope vector combination), so it's a clean, exact check
  of pyflam's ``spread_field(...).heading``.
* ``BURN_PROB.tif`` (fraction) — burn probability. We reproduce it statistically
  with ``pyflam.burn_probability`` (random ignitions, 60-min fires, one prebuilt
  graph). It cannot match cell-for-cell: FlamMap's random seed differs and it
  models ember spotting (spot prob 0.1) that pyflam does not. We compare the
  pattern after block-aggregation (which averages out Monte-Carlo noise) and the
  mean burn probability / typical fire size.

Usage (from the pyflam dir):
    PYTHONPATH=src python tests/validate_flammap_mtt.py --direction
    PYTHONPATH=src python tests/validate_flammap_mtt.py --burnprob --nfires 2000
    # vary the weather fire-to-fire (FlamMap-like) over a 7-scenario ensemble:
    PYTHONPATH=src python tests/validate_flammap_mtt.py --burnprob --ensemble 7
    # condition dead fuel moisture per cell (terrain insolation + canopy):
    PYTHONPATH=src python tests/validate_flammap_mtt.py --burnprob --condition \
        --air-temp 30 --air-rh 20 --datetime "2025-07-29 14:00"
    # ... or derive the weather from a GFS run (needs --longitude + network):
    PYTHONPATH=src python tests/validate_flammap_mtt.py --burnprob --condition \
        --gfs "2025-07-29 12:00" --gfs-fxx 2 --longitude 11 --latitude 43
    # compare the MTT and fast_marching engines on the real terrain:
    PYTHONPATH=src python tests/validate_flammap_mtt.py --compare-methods --window 200
"""

from __future__ import annotations

import argparse
import sys
import time
from datetime import datetime, timedelta

import numpy as np

import pyflam
from pyflam import validate
from pyflam.units import btu_per_ft_s_to_kw_per_m, mph_to_ft_per_min

# Reuse the landscape/moisture loaders from the ROS validation script.
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from validate_flammap_ros import (  # noqa: E402
    DATA, DEFAULT_LANDSCAPE, DEFAULT_FMS, load_landscape, read_fms_default,
)

DEFAULT_RUN = f"{DATA}/FlamMap_output_v1/burn_prob_test_6___forse_buono"

# From run_log.txt: wind 30 km/h @ 10 m. The ROS validation found this matches a
# ~16 mph 20-ft wind (10 m -> 20 ft log reduction); direction 45 deg FROM.
RUN_WIND_MPH = 16.0
RUN_WIND_DIR = 45.0
RUN_SIM_TIME = 60.0     # MTT_SIM_TIME, minutes
RUN_NFIRES = 2000


def build_field(ls, moist, wind_mph, wind_dir, load_factor=1.0):
    midflame = pyflam.midflame_field(ls, mph_to_ft_per_min(wind_mph))
    return pyflam.spread_field(
        ls, wind_midflame=midflame, wind_direction=wind_dir,
        load_factor=load_factor, **moist)


def build_ensemble(ls, moist, wind_mph, wind_dir, n_scenarios, wind_sd, dir_sd,
                   rng, load_factor=1.0):
    """A weather ensemble: ``n_scenarios`` equal-weight spread fields whose wind
    speed and direction are drawn around the nominal values.

    FlamMap's random-ignition run varies the weather fire-to-fire; modelling that
    here (rather than one deterministic weather) is what makes ``burn_probability``
    return a non-degenerate flame-length distribution (FLP) and spreads burned area
    across more of the landscape. Returns ``(scenarios, speeds, dirs)`` where
    ``scenarios`` is the list of scenario dicts :func:`burn_probability` consumes;
    each dict carries its own ``wind_20ft`` / ``wind_direction`` so that with
    ``--spotting`` the ember transport uses that scenario's wind too. One
    travel-time graph is built (lazily) per scenario, so keep ``n_scenarios`` small
    on big landscapes.
    """
    speeds = np.clip(rng.normal(wind_mph, wind_sd, n_scenarios), 1.0, None)
    dirs = rng.normal(wind_dir, dir_sd, n_scenarios) % 360.0
    scenarios = [
        {"weight": 1.0,
         "field": build_field(ls, moist, s, d, load_factor=load_factor),
         "wind_20ft": mph_to_ft_per_min(s), "wind_direction": d}
        for s, d in zip(speeds, dirs)
    ]
    return scenarios, speeds, dirs


def weather_provider(args):
    """Build a meteo provider from --gfs if requested, else None (manual T/RH)."""
    if not args.gfs:
        return None, "manual T/RH"
    try:
        provider = pyflam.atmosphere.fetch_gfs(run=args.gfs, fxx=args.gfs_fxx)
    except Exception as exc:                      # network / missing deps / bad run
        print(f"  (could not fetch GFS [{exc}]; falling back to manual T/RH)")
        return None, "manual T/RH (GFS fetch failed)"
    return provider, f"GFS {args.gfs} +{args.gfs_fxx}h"


def conditioned_moisture(ls, moist, args):
    """Replace the scalar dead 1/10/100-h moistures with per-cell conditioned ones.

    Derives dead fuel moisture from weather for a date/time/location -- pulled from
    a meteo provider (``--gfs``) when available, else from the manual --air-temp /
    --air-rh -- then spreads it across the landscape with terrain insolation +
    canopy shading (:func:`pyflam.condition_from_weather`). Sun-exposed (south, open)
    cells come out drier and shaded (north, canopy) cells moister, instead of one
    landscape-wide value. The .fms live-fuel moistures are kept.
    """
    if args.datetime:
        when = args.datetime
    else:
        when = datetime(2025, 1, 1) + timedelta(
            days=args.day_of_year - 1, hours=args.hour)
    provider, source = weather_provider(args)
    cond = pyflam.condition_from_weather(
        ls, time=when, atmosphere=provider,
        temperature=args.air_temp, relative_humidity=args.air_rh,
        latitude=args.latitude, longitude=args.longitude, timezone=args.timezone,
        model=args.moisture_model, dt_sun=args.dt_sun, canopy=not args.no_canopy)

    burnable = pyflam.spread_field(ls, wind_midflame=0.0, **moist).ros_max > 0.0
    m1 = 100.0 * cond["m_1h"][burnable]
    has_canopy = ls.canopy_cover is not None
    print(f"CONDITIONING dead fuel moisture ({args.moisture_model}; weather: {source}; "
          f"time {when}, lat {args.latitude}"
          f"{'' if has_canopy else '; no canopy band'}):")
    print(f"  per-cell 1-h moisture  min {m1.min():.1f}%  mean {m1.mean():.1f}%  "
          f"max {m1.max():.1f}%  (was scalar {100*moist['m_1h']:.1f}%)\n")
    out = dict(moist)
    out.update(cond)                    # m_1h/m_10h/m_100h now per-cell arrays
    return out


def read_band(path):
    import rasterio
    with rasterio.open(path) as ds:
        a = ds.read(1).astype(float)
        nd = ds.nodata
    valid = np.isfinite(a)
    if nd is not None:
        valid &= a != nd
    return a, valid


def validate_direction(ls, moist, run_dir, wind_mph, wind_dir):
    path = f"{run_dir}/MAX_SPRE_DIR.tif"
    ref_rad, valid = read_band(path)
    ref_deg = np.degrees(np.where(valid, ref_rad, np.nan)) % 360.0

    field = build_field(ls, moist, wind_mph, wind_dir)
    # Direction is only meaningful where the fire actually spreads.
    mask = valid & (field.ros_max > 0.0)
    cmp = validate.compare_directions(field.heading, ref_deg, mask=mask)
    print("MAX SPREAD DIRECTION  (pyflam spread_field.heading vs FlamMap):")
    print(cmp.summary())
    return cmp


def validate_burnprob(ls, moist, run_dir, wind_mph, wind_dir, nfires, sim_time,
                      seed, block, load_factor=1.0, spot_model=None,
                      ensemble=0, wind_sd=3.0, dir_sd=15.0):
    import rasterio
    ref, valid = read_band(f"{run_dir}/BURN_PROB.tif")

    # The nominal field fixes the burnable mask + ignition pool (and is the single
    # weather unless an ensemble is requested).
    field = build_field(ls, moist, wind_mph, wind_dir, load_factor=load_factor)
    burnable = field.ros_max > 0.0

    rng = np.random.default_rng(seed)
    rows, cols = np.where(burnable)
    pick = rng.integers(0, rows.size, size=nfires)
    ignitions = list(zip(rows[pick].tolist(), cols[pick].tolist()))

    spot_desc = type(spot_model).__name__ if spot_model is not None else "no spotting"
    if ensemble:
        scenarios, speeds, dirs = build_ensemble(
            ls, moist, wind_mph, wind_dir, ensemble, wind_sd, dir_sd, rng,
            load_factor=load_factor)
        bp_field, bp_graph = scenarios, None
        weather = (f"weather ensemble of {ensemble} scenarios "
                   f"(wind {speeds.min():.0f}-{speeds.max():.0f} mph, "
                   f"dir {dirs.min():.0f}-{dirs.max():.0f} deg)")
    else:
        bp_field, bp_graph = field, pyflam.mtt.build_traveltime_graph(field)
        weather = "single weather"
    print(f"running {nfires} fires of {sim_time:.0f} min "
          f"({spot_desc}, {weather}, load_factor={load_factor}) ...")
    t = time.time()
    res = pyflam.burn_probability(
        bp_field, ignitions, max_time=sim_time, graph=bp_graph, spotting=spot_model,
        wind_20ft=mph_to_ft_per_min(wind_mph), wind_direction=wind_dir,
        fuel_moisture=moist["m_1h"], rng=rng, return_metrics=True)
    prob, nf = res.burn_prob, res.n_fires
    print(f"  done in {time.time()-t:.1f}s ({nf} fires)\n")

    m = valid & burnable
    pf, rf = prob[m], ref[m]
    note = ("spotting on, calibrated params" if spot_model is not None
            else "no spotting")
    print(f"BURN PROBABILITY  (statistical — different seeds/fire-count; {note}):")
    print(f"  mean burn prob   pyflam {pf.mean():.5f}   FlamMap {rf.mean():.5f}")
    print(f"  cells ever burned pyflam {100*(pf>0).mean():.2f}%   "
          f"FlamMap {100*(rf>0).mean():.2f}%")
    if res.fire_sizes.size:
        sizes_ha = res.fire_sizes / 1e4   # landscape units are metres -> hectares
        print(f"  fire size (ha)   pyflam mean {sizes_ha.mean():.1f}  "
              f"median {np.median(sizes_ha):.1f}  max {sizes_ha.max():.1f}")
    # Block-aggregate to average out Monte-Carlo noise, then correlate.
    pb = _block_mean(prob, block)
    rb = _block_mean(np.where(valid, ref, 0.0), block)
    vb = _block_mean(m.astype(float), block) > 0.5
    if vb.sum() > 1:
        r = float(np.corrcoef(pb[vb], rb[vb])[0, 1])
        print(f"  spatial correlation at {block*100:.0f} m blocks: r = {r:.3f}")

    # Connected metric: conditional fireline intensity vs FlamMap FIRE_LINE_INT.
    fli_path = f"{run_dir}/FIRE_LINE_INT.tif"
    try:
        fli_ref, fli_valid = read_band(fli_path)
    except (rasterio.errors.RasterioIOError, OSError):
        return
    # pyflam intensity is Btu/ft/s; FlamMap metric output is kW/m.
    cint_kw = btu_per_ft_s_to_kw_per_m(res.conditional_intensity)
    cm = m & fli_valid & np.isfinite(res.conditional_intensity) & (fli_ref > 0)
    if cm.sum() > 1:
        print("\nCONDITIONAL FIRELINE INTENSITY  (given burned; kW/m):")
        print(f"  mean             pyflam {cint_kw[cm].mean():.0f}   "
              f"FlamMap {fli_ref[cm].mean():.0f}")
        print(f"  conditional flame length (ft) pyflam mean "
              f"{np.nanmean(res.conditional_flame_length[m]):.2f}")
        cb = _block_mean(np.where(cm, cint_kw, 0.0), block)
        rbf = _block_mean(np.where(cm, fli_ref, 0.0), block)
        vbf = _block_mean(cm.astype(float), block) > 0.5
        if vbf.sum() > 1:
            rr = float(np.corrcoef(cb[vbf], rbf[vbf])[0, 1])
            print(f"  spatial correlation at {block*100:.0f} m blocks: r = {rr:.3f}")

    # Connected metric: the flame-length probability distribution (FlamMap FLP).
    # Averaged over burned cells, this is only non-degenerate when the weather
    # varies fire-to-fire (the ensemble) -- under one weather every burned cell
    # sits in a single class.
    burned = res.burn_prob > 0
    if burned.any():
        centers = res.flame_length_class_centers()
        # Mean FLP profile over burned cells (weight classes by how often a cell burns).
        w = res.burn_prob[burned]
        profile = (res.flp[:, burned] * w).sum(axis=1) / w.sum()
        occupied = int((profile > 1e-4).sum())
        print(f"\nFLAME-LENGTH PROBABILITY (FLP, mean over burned cells; "
              f"{occupied}/{centers.size} classes occupied):")
        top = np.argsort(profile)[::-1][:5]
        for k in sorted(top):
            lo = res.flame_length_classes[k]
            hi = res.flame_length_classes[k + 1]
            label = f"{lo:.0f}+ ft" if not np.isfinite(hi) else f"{lo:.0f}-{hi:.0f} ft"
            print(f"  {label:>9} : {100*profile[k]:5.1f}%")


def compare_methods(ls, moist, wind_mph, wind_dir, sim_time, window,
                    load_factor=1.0):
    """Head-to-head of the two propagation engines on the real landscape.

    The Tuscany dataset has no FlamMap time-of-arrival raster to validate against,
    so this compares pyflam's own ``method="mtt"`` (Dijkstra) and
    ``method="fast_marching"`` (anisotropic Eikonal) on the *same* spread field --
    a single fire grown for ``sim_time`` minutes on a ``window``-cell crop around a
    central ignition (the Eikonal solver sweeps the whole window, so it is bounded
    for tractability). Reports timing-agreement and burned-area overlap.
    """
    from scipy.ndimage import uniform_filter
    field = build_field(ls, moist, wind_mph, wind_dir, load_factor=load_factor)
    nr, nc = field.shape
    burnable = field.ros_max > 0.0
    # Centre the window on the most fuel-dense region so the fire grows through
    # contiguous fuel (the landscape centre is often fragmented / barrier-heavy).
    dens = uniform_filter(burnable.astype(np.float32), size=window, mode="constant")
    cr, cc = np.unravel_index(int(np.argmax(dens)), dens.shape)
    h = window // 2
    r0, r1 = max(0, cr - h), min(nr, cr + h)
    c0, c1 = max(0, cc - h), min(nc, cc + h)
    sub = pyflam.SpreadField(
        ros_max=np.ascontiguousarray(field.ros_max[r0:r1, c0:c1]),
        eccentricity=np.ascontiguousarray(field.eccentricity[r0:r1, c0:c1]),
        heading=np.ascontiguousarray(field.heading[r0:r1, c0:c1]),
        cellsize_x=field.cellsize_x, cellsize_y=field.cellsize_y)
    if sub.ros_max.max() <= 0.0:
        print("COMPARE ENGINES: no burnable fuel in the window; skipping\n")
        return
    sir, sic = np.unravel_index(int(np.argmax(sub.ros_max)), sub.shape)
    ign = [(int(sir), int(sic))]
    print(f"COMPARE ENGINES  (single fire, {sub.shape} window, {sim_time:.0f} min; "
          f"no FlamMap TOA raster -> the two engines vs each other):")

    t = time.time()
    t_mtt = pyflam.minimum_travel_time(sub, ign, max_time=sim_time, method="mtt")
    dt_mtt = time.time() - t
    t = time.time()
    t_fm = pyflam.minimum_travel_time(
        sub, ign, max_time=sim_time, method="fast_marching")
    dt_fm = time.time() - t

    bm, bf = t_mtt <= sim_time, t_fm <= sim_time
    print(f"  burned cells   MTT {int(bm.sum())} ({dt_mtt*1e3:.0f} ms)   "
          f"fast_marching {int(bf.sum())} ({dt_fm*1e3:.0f} ms)")
    both = bm & bf
    n = int(both.sum())
    if n > 1:
        d = np.abs(t_fm[both] - t_mtt[both])            # arrival-time |Delta|, minutes
        r = float(np.corrcoef(t_fm[both], t_mtt[both])[0, 1])
        print(f"  arrival time (min) over {n} cells both burn: "
              f"mean |dt| {d.mean():.2f}  median {np.median(d):.2f}  max {d.max():.2f}"
              f"   correlation r = {r:.3f}")
    perim = validate.compare_perimeters(
        bf, bm, cellsize_x=field.cellsize_x, cellsize_y=field.cellsize_y)
    print(f"  burned-area overlap: Jaccard {perim.jaccard:.3f}  Dice {perim.dice:.3f}  "
          f"Hausdorff {perim.hausdorff:.0f} m")
    print("  note: MTT prunes via its heap+limit; the Eikonal sweep covers the whole "
          "window, so it is slower for a small bounded fire.")


def _block_mean(a, k):
    nr, nc = a.shape
    nr2, nc2 = (nr // k) * k, (nc // k) * k
    return a[:nr2, :nc2].reshape(nr2 // k, k, nc2 // k, k).mean(axis=(1, 3))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--landscape", default=DEFAULT_LANDSCAPE)
    ap.add_argument("--run-dir", default=DEFAULT_RUN)
    ap.add_argument("--fms", default=DEFAULT_FMS)
    ap.add_argument("--wind-mph", type=float, default=RUN_WIND_MPH)
    ap.add_argument("--wind-dir", type=float, default=RUN_WIND_DIR)
    ap.add_argument("--direction", action="store_true", help="check MAX_SPR_DIR")
    ap.add_argument("--burnprob", action="store_true", help="check BURN_PROB")
    ap.add_argument("--compare-methods", action="store_true",
                    help="MTT vs fast_marching engine head-to-head on a window")
    ap.add_argument("--window", type=int, default=200,
                    help="crop size (cells) for --compare-methods")
    ap.add_argument("--nfires", type=int, default=RUN_NFIRES)
    ap.add_argument("--sim-time", type=float, default=RUN_SIM_TIME)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--block", type=int, default=20, help="block size (cells)")
    ap.add_argument("--load-factor", type=float, default=1.0,
                    help="fuel load multiplier (e.g. 1.3 for +30%%)")
    ap.add_argument("--spotting", action="store_true", help="enable ember spotting")
    ap.add_argument("--physics", action="store_true",
                    help="use the stochastic physics-based firebrand model")
    ap.add_argument("--front-length", type=float, default=50.0,
                    help="coherent flaming-front length (m) for the plume source")
    ap.add_argument("--embers-per-source", type=float, default=3.0)
    ap.add_argument("--spot-prob", type=float, default=0.1,
                    help="spot probability (run_log: 0.1)")
    ap.add_argument("--loft-coeff", type=float, default=12.0)
    ap.add_argument("--terminal-velocity", type=float, default=250.0)
    ap.add_argument("--launch-fraction", type=float, default=0.02)
    ap.add_argument("--spot-delay", type=float, default=2.0,
                    help="minutes (run_log: 120 s)")
    ap.add_argument("--ensemble", type=int, default=0, metavar="N",
                    help="vary weather fire-to-fire over N scenarios (0 = single "
                         "weather); makes the FLP distribution non-degenerate")
    ap.add_argument("--wind-sd", type=float, default=3.0,
                    help="ensemble wind-speed spread (mph, 1 sigma)")
    ap.add_argument("--dir-sd", type=float, default=15.0,
                    help="ensemble wind-direction spread (deg, 1 sigma)")
    # --- per-cell dead fuel moisture conditioning (terrain + canopy) ---
    ap.add_argument("--condition", action="store_true",
                    help="condition dead fuel moisture per cell from T/RH + "
                         "slope/aspect/elevation/canopy + sun position")
    ap.add_argument("--air-temp", type=float, default=28.0,
                    help="air temperature for conditioning (C)")
    ap.add_argument("--air-rh", type=float, default=25.0,
                    help="air relative humidity for conditioning (%%)")
    ap.add_argument("--latitude", type=float, default=43.0,
                    help="site latitude for the sun position (deg, Tuscany ~43)")
    ap.add_argument("--day-of-year", type=int, default=210)
    ap.add_argument("--hour", type=float, default=14.0,
                    help="solar hour, or local clock hour if --longitude+--tz set")
    ap.add_argument("--datetime", default=None,
                    help="'YYYY-MM-DD HH:MM' run time (overrides --day-of-year/--hour)")
    ap.add_argument("--gfs", default=None, metavar="RUN",
                    help="pull T/RH from a GFS run ('YYYY-MM-DD HH:MM') instead of "
                         "--air-temp/--air-rh; needs --longitude and network")
    ap.add_argument("--gfs-fxx", type=int, default=0,
                    help="GFS forecast hour for --gfs")
    ap.add_argument("--longitude", type=float, default=None,
                    help="site longitude (deg E); with --tz, --hour is clock time")
    ap.add_argument("--tz", dest="timezone", type=float, default=None,
                    help="UTC offset hours (e.g. 2 for CEST)")
    ap.add_argument("--moisture-model", choices=("emc", "vpd"), default="emc")
    ap.add_argument("--dt-sun", type=float, default=17.0,
                    help="max near-fuel solar heating over air (C)")
    ap.add_argument("--no-canopy", action="store_true",
                    help="ignore canopy shading in conditioning")
    args = ap.parse_args(argv)

    try:
        import rasterio  # noqa: F401
    except ImportError:
        sys.exit("This script needs rasterio: pip install 'pyflam[geo]'")

    if not (args.direction or args.burnprob or args.compare_methods):
        args.direction = args.burnprob = True   # default: both reference checks

    moist = read_fms_default(args.fms)
    ls = load_landscape(args.landscape)
    print(f"landscape {ls.shape}; wind {args.wind_mph} mph / {args.wind_dir} deg\n")

    if args.compare_methods:
        compare_methods(ls, moist, args.wind_mph, args.wind_dir, args.sim_time,
                        args.window, load_factor=args.load_factor)
        print()
    if args.direction:
        validate_direction(ls, moist, args.run_dir, args.wind_mph, args.wind_dir)
        print()
    if args.burnprob:
        spot_model = None
        if args.spotting and args.physics:
            spot_model = pyflam.FirebrandPhysics(
                front_length=args.front_length,
                embers_per_source=args.embers_per_source,
                launch_fraction=args.launch_fraction, spot_delay=args.spot_delay)
        elif args.spotting:
            spot_model = pyflam.SpottingModel(
                spot_probability=args.spot_prob, loft_coeff=args.loft_coeff,
                terminal_velocity=args.terminal_velocity,
                launch_fraction=args.launch_fraction, spot_delay=args.spot_delay)
        bp_moist = conditioned_moisture(ls, moist, args) if args.condition else moist
        validate_burnprob(ls, bp_moist, args.run_dir, args.wind_mph, args.wind_dir,
                          args.nfires, args.sim_time, args.seed, args.block,
                          load_factor=args.load_factor, spot_model=spot_model,
                          ensemble=args.ensemble, wind_sd=args.wind_sd,
                          dir_sd=args.dir_sd)


if __name__ == "__main__":
    main()
