#!/usr/bin/env python
"""Validate pyflam's MTT growth against a FlamMap single-fire perimeter / arrival time.

This is the clean growth-engine acceptance test: run FlamMap MTT for a *single
fire with spotting OFF*, export its arrival-time (time-of-arrival) raster, then
diff pyflam's MTT against it -- perimeter overlap (Jaccard / Dice / Hausdorff) and
arrival-time agreement where both fires reach.

The bundled Tuscany dataset only has spotting-*on* burn-probability runs, so this
script is parameterized for the export you produce:

    PYTHONPATH=src python tests/validate_flammap_perimeter.py \
        --toa flammap_arrival.tif --toa-units minutes \
        --ignition-xy 4350000 2350000 --wind-mph 16 --wind-dir 45 \
        --perimeter-time 120

It loads the Tuscany landscape (or ``--landscape``), runs pyflam MTT (no spotting)
at the same inputs, and prints the perimeter + arrival-time comparison.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

import pyflam
from pyflam import validate
from pyflam.units import mph_to_ft_per_min

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from validate_flammap_ros import (  # noqa: E402
    DEFAULT_LANDSCAPE, DEFAULT_FMS, load_landscape, read_fms_default,
)


def load_arrival(path, units="minutes"):
    """Read a FlamMap time-of-arrival raster -> arrival minutes (+ valid mask)."""
    import rasterio
    with rasterio.open(path) as ds:
        toa = ds.read(1).astype(float)
        nodata = ds.nodata
    valid = np.isfinite(toa)
    if nodata is not None:
        valid &= toa != nodata
    valid &= toa >= 0
    factor = 60.0 if units.startswith("hour") else 1.0
    arrival = np.where(valid, toa * factor, np.inf)
    return arrival, valid


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--landscape", default=DEFAULT_LANDSCAPE)
    ap.add_argument("--fms", default=DEFAULT_FMS)
    ap.add_argument("--toa", required=True, help="FlamMap time-of-arrival raster")
    ap.add_argument("--toa-units", default="minutes", choices=["minutes", "hours"])
    ap.add_argument("--ignition-xy", type=float, nargs=2, metavar=("X", "Y"),
                    help="ignition world coordinates")
    ap.add_argument("--ignition-rc", type=int, nargs=2, metavar=("ROW", "COL"))
    ap.add_argument("--wind-mph", type=float, default=16.0)
    ap.add_argument("--wind-dir", type=float, default=45.0)
    ap.add_argument("--perimeter-time", type=float, default=None,
                    help="minutes at which to compare perimeters (default: max TOA)")
    args = ap.parse_args(argv)

    try:
        import rasterio  # noqa: F401
    except ImportError:
        sys.exit("needs rasterio: pip install 'pyflam[geo]'")

    moist = read_fms_default(args.fms)
    ls = load_landscape(args.landscape)
    ref_arrival, ref_valid = load_arrival(args.toa, args.toa_units)
    if ref_arrival.shape != ls.shape:
        sys.exit(f"grid mismatch: TOA {ref_arrival.shape} vs landscape {ls.shape}")

    if args.ignition_rc:
        ign = tuple(args.ignition_rc)
    elif args.ignition_xy:
        ign = pyflam.ignition_from_xy(ls, *args.ignition_xy)
    else:
        # default: the earliest-arrival cell of the FlamMap run
        ign = np.unravel_index(np.argmin(np.where(ref_valid, ref_arrival, np.inf)),
                               ref_arrival.shape)
    print(f"ignition cell {ign}")

    t_end = (args.perimeter_time
             if args.perimeter_time is not None
             else float(ref_arrival[np.isfinite(ref_arrival)].max()))

    field = pyflam.spread_field(
        ls, wind_midflame=pyflam.midflame_field(ls, mph_to_ft_per_min(args.wind_mph)),
        wind_direction=args.wind_dir, **moist)
    pf_arrival = pyflam.minimum_travel_time(field, [ign], max_time=t_end)

    pf_burned = np.isfinite(pf_arrival) & (pf_arrival <= t_end)
    ref_burned = np.isfinite(ref_arrival) & (ref_arrival <= t_end)
    print(f"\nperimeter at {t_end:g} min:")
    print(validate.compare_perimeters(
        pf_burned, ref_burned, cellsize_x=ls.cellsize_x,
        cellsize_y=ls.cellsize_y).summary())
    print("\narrival-time agreement (cells both reach):")
    print(validate.compare_arrival_times(
        pf_arrival, ref_arrival, max_time=t_end).summary(units="min"))


if __name__ == "__main__":
    main()
