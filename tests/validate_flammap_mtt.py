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
"""

from __future__ import annotations

import argparse
import sys
import time

import numpy as np

import pyflam
from pyflam import validate
from pyflam.units import mph_to_ft_per_min

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
                      seed, block, load_factor=1.0, spot_model=None):
    import rasterio
    ref, valid = read_band(f"{run_dir}/BURN_PROB.tif")

    field = build_field(ls, moist, wind_mph, wind_dir, load_factor=load_factor)
    burnable = field.ros_max > 0.0

    rng = np.random.default_rng(seed)
    rows, cols = np.where(burnable)
    pick = rng.integers(0, rows.size, size=nfires)
    ignitions = list(zip(rows[pick].tolist(), cols[pick].tolist()))

    spot_desc = type(spot_model).__name__ if spot_model is not None else "no spotting"
    print(f"building graph + running {nfires} fires of {sim_time:.0f} min "
          f"({spot_desc}, load_factor={load_factor}) ...")
    t = time.time()
    graph = pyflam.mtt.build_traveltime_graph(field)
    prob, nf = pyflam.burn_probability(
        field, ignitions, max_time=sim_time, graph=graph, spotting=spot_model,
        wind_20ft=mph_to_ft_per_min(wind_mph), wind_direction=wind_dir,
        fuel_moisture=moist["m_1h"], rng=rng)
    print(f"  done in {time.time()-t:.1f}s ({nf} fires)\n")

    m = valid & burnable
    pf, rf = prob[m], ref[m]
    note = ("spotting on, calibrated params" if spot_model is not None
            else "no spotting")
    print(f"BURN PROBABILITY  (statistical — different seeds/fire-count; {note}):")
    print(f"  mean burn prob   pyflam {pf.mean():.5f}   FlamMap {rf.mean():.5f}")
    print(f"  cells ever burned pyflam {100*(pf>0).mean():.2f}%   "
          f"FlamMap {100*(rf>0).mean():.2f}%")
    # Block-aggregate to average out Monte-Carlo noise, then correlate.
    pb = _block_mean(prob, block)
    rb = _block_mean(np.where(valid, ref, 0.0), block)
    vb = _block_mean(m.astype(float), block) > 0.5
    if vb.sum() > 1:
        r = float(np.corrcoef(pb[vb], rb[vb])[0, 1])
        print(f"  spatial correlation at {block*100:.0f} m blocks: r = {r:.3f}")


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
    args = ap.parse_args(argv)

    try:
        import rasterio  # noqa: F401
    except ImportError:
        sys.exit("This script needs rasterio: pip install 'pyflam[geo]'")

    if not (args.direction or args.burnprob):
        args.direction = args.burnprob = True   # default: both

    moist = read_fms_default(args.fms)
    ls = load_landscape(args.landscape)
    print(f"landscape {ls.shape}; wind {args.wind_mph} mph / {args.wind_dir} deg\n")

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
        validate_burnprob(ls, moist, args.run_dir, args.wind_mph, args.wind_dir,
                          args.nfires, args.sim_time, args.seed, args.block,
                          load_factor=args.load_factor, spot_model=spot_model)


if __name__ == "__main__":
    main()
