#!/usr/bin/env python
"""Validate pyflam surface ROS against a real FlamMap ``ROS.tif``.

Wires the user's Tuscany dataset into :mod:`pyflam.validate`. FlamMap's "Rate of
Spread" output is the maximum (head-fire) spread rate, combining wind and slope,
so the matching pyflam quantity is ``spread_field(...).ros_max`` (the vector
wind+slope combination), not the scalar ``basic_fire_behavior``.

The landscape is read from the 8-band ``Toscana.tif`` FlamMap actually ran on
(Elevation, Slope[deg], Aspect[deg], Fuel Model, Canopy Cover, then zeroed crown
bands — so this validates *surface* ROS only). FlamMap ROS is in chains/hour;
we convert to ft/min. Fuel moisture is read from the ``.fms`` file.

The one input a FlamMap project file doesn't expose in an easily parsed form is
the wind. So by default this runs at a stated wind AND offers ``--scan`` to sweep
wind speeds and report which best reproduces FlamMap — a diagnostic to pin down
the run's setting (clearly distinct from the validation itself).

Usage (from the pyflam dir):
    PYTHONPATH=src python tests/validate_flammap_ros.py --wind-mph 12 --wind-dir 45
    PYTHONPATH=src python tests/validate_flammap_ros.py --scan 2,4,6,8,12,16,20,25,30
    PYTHONPATH=src python tests/validate_flammap_ros.py --wind-mph 12 --diff-out diff.tif
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

import pyflam
from pyflam import validate
from pyflam.units import chains_per_hour_to_ft_per_min, mph_to_ft_per_min

DATA = "/Users/cristianofoderi/DATI/FUEL_TOS/Toscana_Fuel_data"
DEFAULT_LANDSCAPE = f"{DATA}/FlamMap_output_v1/Toscana.tif"
DEFAULT_ROS = f"{DATA}/FlamMap_output_v1/burn_prob_test_6___forse_buono/ROS.tif"
DEFAULT_FMS = f"{DATA}/FlamMap_output_v1/moisture_prova_1.fms"

_INT_MIN = -2147483648  # FlamMap/GDAL nodata fill in the integer bands


def read_fms_default(path: str) -> dict:
    """Read fuel moisture from a FlamMap ``.fms`` file (percent -> fraction).

    Uses the fuel-0 (default) row if present, else the first row. Columns are
    ``fuel 1h 10h 100h live_herb live_woody [...]``.
    """
    rows = {}
    with open(path) as fh:
        for line in fh:
            p = line.split()
            if len(p) >= 6:
                rows[int(float(p[0]))] = [float(x) for x in p[1:6]]
    vals = rows.get(0) or next(iter(rows.values()))
    m1, m10, m100, herb, woody = vals
    return dict(m_1h=m1 / 100, m_10h=m10 / 100, m_100h=m100 / 100,
                m_live_herb=herb / 100, m_live_woody=woody / 100)


def load_landscape(path: str):
    """Build a pyflam Landscape from the 8-band FlamMap Toscana.tif."""
    import rasterio
    with rasterio.open(path) as ds:
        slope = ds.read(2).astype(float)
        aspect = ds.read(3).astype(float)
        fuel = ds.read(4).astype(np.int32)
        cover = ds.read(5).astype(float)
        height = ds.read(6).astype(float)   # 0 throughout this dataset
        tr = ds.transform
        crs = ds.crs
    # Integer-band nodata -> benign values; such cells are nonburnable anyway.
    slope[slope == _INT_MIN] = 0.0
    aspect[aspect == _INT_MIN] = 0.0
    cover[cover == _INT_MIN] = 0.0
    ls = pyflam.Landscape(
        fuel_model=fuel, slope=slope, aspect=aspect, canopy_cover=cover,
        canopy_height=height,
        cellsize_x=tr.a, cellsize_y=-tr.e, west=tr.c, north=tr.f,
        slope_units="degrees", crs=crs,
    )
    return ls


def load_flammap_ros(path: str):
    """Read FlamMap ROS.tif -> (ros_ft_per_min array, valid mask)."""
    import rasterio
    with rasterio.open(path) as ds:
        ros = ds.read(1).astype(float)
        nodata = ds.nodata
    valid = np.isfinite(ros)
    if nodata is not None:
        valid &= ros != nodata
    ros_ftmin = chains_per_hour_to_ft_per_min(np.where(valid, ros, np.nan))
    return ros_ftmin, valid


def pyflam_ros(ls, moist, wind_mph, wind_dir, wrf, waf_mode="auto"):
    """pyflam max-spread ROS (ft/min) for a 20-ft wind reduced to midflame.

    ``waf_mode="flat"`` uses the single ``wrf`` everywhere; ``"auto"`` uses the
    per-cell Albini & Baughman / Andrews wind adjustment factor (unsheltered from
    fuel depth here, since the canopy height band is zero).
    """
    wind_20ft = mph_to_ft_per_min(wind_mph)
    if waf_mode == "flat":
        midflame = wind_20ft * wrf
    else:
        midflame = pyflam.midflame_field(ls, wind_20ft)
    field = pyflam.spread_field(
        ls, wind_midflame=midflame, wind_direction=wind_dir, **moist)
    return field.ros_max


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--landscape", default=DEFAULT_LANDSCAPE)
    ap.add_argument("--ros", default=DEFAULT_ROS, help="FlamMap ROS.tif")
    ap.add_argument("--fms", default=DEFAULT_FMS, help="moisture .fms file")
    ap.add_argument("--wind-mph", type=float, default=12.0,
                    help="20-ft wind speed (mi/h)")
    ap.add_argument("--wind-dir", type=float, default=45.0,
                    help="wind direction (deg FROM, met)")
    ap.add_argument("--wrf", type=float, default=0.4,
                    help="flat wind reduction factor (only with --waf flat)")
    ap.add_argument("--waf", choices=["auto", "flat"], default="auto",
                    help="per-cell WAF (auto) or a single flat factor")
    ap.add_argument("--scan", default=None,
                    help="comma-separated wind speeds (mph) to sweep for best fit")
    ap.add_argument("--diff-out", default=None,
                    help="write a (pyflam - FlamMap) diff GeoTIFF here")
    args = ap.parse_args(argv)

    try:
        import rasterio  # noqa: F401
    except ImportError:
        sys.exit("This script needs rasterio: pip install 'pyflam[geo]'")

    moist = read_fms_default(args.fms)
    print(f"moisture (fractions): {moist}")
    ls = load_landscape(args.landscape)
    ros_ref, valid = load_flammap_ros(args.ros)
    if ros_ref.shape != ls.shape:
        sys.exit(f"grid mismatch: ROS {ros_ref.shape} vs landscape {ls.shape}")
    print(f"landscape {ls.shape}, {valid.sum():,} valid FlamMap cells\n")

    if args.scan:
        speeds = [float(x) for x in args.scan.split(",")]
        best, best_cmp, table = validate.scan_parameter(
            lambda u: pyflam_ros(ls, moist, u, args.wind_dir, args.wrf, args.waf),
            speeds, ros_ref, mask=valid,
        )
        print("wind-speed scan (diagnostic — finds the best-fit wind):")
        print(f"  {'mph':>6} {'log_rmse':>9} {'med_ratio':>10} {'within25%':>10}")
        for u, c in table:
            mark = "  <-- best" if u == best else ""
            print(f"  {u:6.1f} {c.log_rmse:9.4f} {c.median_ratio:10.3f} "
                  f"{100*c.within_25pct:9.1f}%{mark}")
        print(f"\nbest-fit wind ~ {best} mph (WAF={args.waf}); "
              f"comparison at that wind:\n")
        print(best_cmp.summary())
        return

    cmp = validate.compare_fields(
        pyflam_ros(ls, moist, args.wind_mph, args.wind_dir, args.wrf, args.waf),
        ros_ref, mask=valid,
    )
    waf_desc = f"WAF {args.wrf} flat" if args.waf == "flat" else "WAF auto (per-fuel)"
    print(f"pyflam at {args.wind_mph} mph / {args.wind_dir} deg / {waf_desc}:")
    print(cmp.summary())

    if args.diff_out:
        import rasterio
        pf = pyflam_ros(ls, moist, args.wind_mph, args.wind_dir, args.wrf, args.waf)
        diff = np.where(valid, pf - ros_ref, np.nan).astype("float32")
        with rasterio.open(args.ros) as src:
            profile = src.profile
        profile.update(dtype="float32", count=1, nodata=np.nan)
        with rasterio.open(args.diff_out, "w", **profile) as dst:
            dst.write(diff, 1)
        print(f"\nwrote diff raster -> {args.diff_out}")


if __name__ == "__main__":
    main()
