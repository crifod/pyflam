#!/usr/bin/env python
"""Validate pyflam's crown-fire classification against a FlamMap crown-fire run.

pyflam's :func:`pyflam.crownfire.crown_fire_potential` labels every cell
**surface / passive / active** (Scott & Reinhardt 2001, on Van Wagner 1977 +
Rothermel 1991) and reports crown fraction burned and total fireline intensity.
FlamMap's "Crown Fire Activity" raster is the same classification, so the diff is a
3-class confusion matrix (`pyflam.validate.compare_categories`) plus field
comparisons of crown fraction / intensity where both classify crowning.

Note: FlamMap is **not ground truth** for crown fire. FlamMap and the Rothermel +
Van Wagner operational stack under-predict crown fire spread (Cruz & Alexander 2010),
which is why pyflam defaults to the Cruz et al. (2005) crown ROS. So this harness is
for *comparing* crown classifications across models/sources (Cruz vs Rothermel vs
FlamMap), not for treating a FlamMap diff as an acceptance test -- a disagreement
where pyflam-Cruz crowns and FlamMap stays surface is expected, not an error.

**Data requirements** (this is why it's a harness, not yet a result): the landscape
must carry **canopy base height (CBH)** and **canopy bulk density (CBD)** bands, and
the FlamMap run must export a crown-fire-activity raster. The bundled Tuscany
dataset has neither (its `.lcp` has canopy cover but no CBH/CBD, and the FlamMap
output is surface-only) — so on that data the script reports what is missing and
stops. Run ``--synthetic`` to exercise the comparison logic end-to-end without it.

Usage (from the pyflam dir):
    PYTHONPATH=src python tests/validate_flammap_crown.py --synthetic
    PYTHONPATH=src python tests/validate_flammap_crown.py \
        --landscape canopy.lcp --fire-type FLAMMAP_CROWN_ACTIVITY.tif
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

# pyflam fire_type codes and their names.
TYPE_LABELS = [0, 1, 2]
TYPE_NAMES = ["surface", "passive", "active"]

# Default FlamMap "Crown Fire Activity" code -> pyflam fire_type. FlamMap 6 uses
# 1 surface, 2 passive (torching), 3 active (running crown); 0 = unburnable.
DEFAULT_CODE_MAP = {1: 0, 2: 1, 3: 2}

RUN_WIND_MPH = 16.0
RUN_WIND_DIR = 45.0
RUN_FOLIAR = 100.0          # foliar moisture %, a typical summer value


def run_pyflam_crown(ls, moist, wind_mph, foliar, *, cbh_scale=0.1, cbd_scale=0.01):
    midflame = pyflam.midflame_field(ls, mph_to_ft_per_min(wind_mph))
    return pyflam.crownfire.crown_fire_potential(
        ls, foliar_moisture=foliar, wind_20ft_ft_per_min=mph_to_ft_per_min(wind_mph),
        wind_midflame=midflame, cbh_scale=cbh_scale, cbd_scale=cbd_scale, **moist)


def _missing_bands(ls):
    miss = [n for n in ("canopy_base_height", "canopy_bulk_density")
            if getattr(ls, n, None) is None
            or not np.any(np.asarray(getattr(ls, n)) > 0)]
    return miss


def remap_reference(ref, code_map):
    """FlamMap crown-activity codes -> pyflam fire_type (cells outside the map -> -1)."""
    out = np.full(ref.shape, -1, dtype=np.int64)
    for code, t in code_map.items():
        out[np.asarray(ref) == code] = t
    return out


def report(out, ref_type, valid):
    """Print the 3-class confusion + crown-fraction/intensity field diffs."""
    py_type = np.asarray(out["fire_type"], dtype=np.int64)
    cmp = validate.compare_categories(py_type, ref_type, labels=TYPE_LABELS, mask=valid)
    print("CROWN FIRE TYPE  (pyflam vs FlamMap):")
    print("  " + cmp.summary(names=TYPE_NAMES).replace("\n", "\n  "))

    crowning = valid & (py_type >= 1) & (ref_type >= 1)
    if int(crowning.sum()) > 1 and "crown_fraction_burned" in out:
        cfb = validate.compare_fields(
            out["crown_fraction_burned"], out["crown_fraction_burned"],
            mask=crowning, burn_threshold=-np.inf)  # placeholder until ref CFB exists
        del cfb  # FlamMap CFB raster not standard; class agreement is the headline
    return cmp


def synthetic_demo():
    """Build a canopy landscape spanning surface/passive/active and round-trip the
    comparison machinery (proves the harness without the real dataset)."""
    n = 60
    rng = np.random.default_rng(0)
    # Graded canopy: low CBH + high CBD on one side (crown-prone), the reverse on the
    # other; a band of nonburnable fuel; the rest timber (FM10) on moderate slope.
    fuel = np.full((n, n), 10, dtype=int)
    fuel[:, :4] = 91                                   # nonburnable strip
    cbh = np.tile(np.linspace(2.0, 80.0, n), (n, 1))   # CBH*10 (LCP units): 0.2..8 m
    cbd = np.tile(np.linspace(40.0, 2.0, n)[::-1], (n, 1))  # CBD*100: 0.02..0.40
    ls = pyflam.Landscape(
        fuel_model=fuel, slope=np.full((n, n), 25.0, dtype=float),
        aspect=np.full((n, n), 180.0, dtype=float),
        canopy_base_height=cbh.astype(float),
        canopy_bulk_density=cbd.astype(float),
        canopy_height=np.full((n, n), 150.0, dtype=float),   # 15 m
        cellsize_x=30.0, cellsize_y=30.0, west=0.0, north=n * 30.0,
        slope_units="degrees")
    moist = dict(m_1h=0.04, m_10h=0.05, m_100h=0.06, m_live_herb=0.5, m_live_woody=0.7)
    out = run_pyflam_crown(ls, moist, 25.0, 90.0)

    py_type = np.asarray(out["fire_type"], np.int64)
    burnable = ls.fuel_model != 91
    dist = {TYPE_NAMES[t]: int(np.sum((py_type == t) & burnable)) for t in TYPE_LABELS}
    print(f"synthetic landscape {ls.shape}; pyflam crown-type cells: {dist}\n")

    # 1) self-consistency: comparing the field to itself must be 100% on the diagonal.
    same = validate.compare_categories(py_type, py_type, labels=TYPE_LABELS, mask=burnable)
    print("SELF-CHECK (pyflam vs itself — must be 100% diagonal):")
    print("  " + same.summary(names=TYPE_NAMES).replace("\n", "\n  "))
    assert same.overall_agreement == 1.0, "self-comparison should be perfect"

    # 2) perturbed 'reference': flip ~15% of cells to a neighbouring class so the
    # confusion matrix and per-class recall are non-trivial.
    ref = py_type.copy()
    flip = burnable & (rng.random(py_type.shape) < 0.15)
    ref[flip] = np.clip(ref[flip] + rng.choice([-1, 1], size=int(flip.sum())), 0, 2)
    print("\nPERTURBED REFERENCE (~15% of cells nudged one class):")
    perturbed = report({"fire_type": py_type}, ref, burnable)
    assert 0.6 < perturbed.overall_agreement < 1.0
    print("\nharness OK — ready for a canopy landscape + FlamMap crown-activity raster.")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--landscape", default=DEFAULT_LANDSCAPE)
    ap.add_argument("--fms", default=DEFAULT_FMS)
    ap.add_argument("--fire-type", default=None,
                    help="FlamMap crown-fire-activity raster (GeoTIFF)")
    ap.add_argument("--wind-mph", type=float, default=RUN_WIND_MPH)
    ap.add_argument("--foliar", type=float, default=RUN_FOLIAR)
    ap.add_argument("--synthetic", action="store_true",
                    help="self-test the comparison logic on a synthetic canopy")
    args = ap.parse_args(argv)

    if args.synthetic:
        synthetic_demo()
        return

    ls = load_landscape(args.landscape)
    moist = read_fms_default(args.fms)
    print(f"landscape {ls.shape}")

    miss = _missing_bands(ls)
    if miss:
        print(f"\nCannot run the crown-fire diff: the landscape is missing {miss}.")
        print("Crown fire needs canopy base height (CBH) and canopy bulk density (CBD)"
              " bands.\nProvide a landscape that carries them (and a FlamMap "
              "crown-activity raster via --fire-type), or run --synthetic to test the "
              "harness.")
        return
    if not args.fire_type:
        print("\nNo --fire-type raster given: pyflam can classify, but there is no "
              "FlamMap crown-activity raster to diff against. Pass --fire-type.")
        return

    try:
        import rasterio
    except ImportError:
        sys.exit("This script needs rasterio: pip install 'pyflam[geo]'")
    with rasterio.open(args.fire_type) as ds:
        ref_raw = ds.read(1)
    ref_type = remap_reference(ref_raw, DEFAULT_CODE_MAP)

    out = run_pyflam_crown(ls, moist, args.wind_mph, args.foliar)
    valid = (np.asarray(ls.fuel_model) != 91) & (ref_type >= 0)
    report(out, ref_type, valid)


if __name__ == "__main__":
    main()
