#!/usr/bin/env python
"""Generate a synthetic canopy landscape (``.lcp``) and exercise the crown path.

The bundled Tuscany dataset has no canopy base-height / bulk-density bands, so the
crown-fire diff (``validate_flammap_crown.py``) and the crown / plume coupling can't
be driven by it. This builds a *synthetic* landscape that DOES carry the full canopy
fuel stack -- a hill of conifer timber with spatially varying canopy base height and
bulk density (low/dense in the south, high/sparse in the north), a nonburnable river
and a grass corner -- so the whole crown pipeline can be run end to end on a real
``.lcp`` file: read -> crown classification (Cruz 2005 vs Rothermel 1991) -> a
crown-aware spread field -> the plume-coupled crown march.

Usage (from the pyflam dir):
    PYTHONPATH=src python tests/make_synthetic_canopy_lcp.py --out canopy.lcp
"""

from __future__ import annotations

import argparse

import numpy as np

import pyflam
from pyflam import validate
from pyflam.units import mph_to_ft_per_min

RUN = dict(m_1h=0.05, m_10h=0.06, m_100h=0.07, m_live_herb=0.6, m_live_woody=0.9)
TYPE_NAMES = ["surface", "passive", "active"]


def make_synthetic_canopy_landscape(n: int = 120, cellsize: float = 30.0, seed: int = 0):
    """A conifer hill with the full canopy stack, designed to span all crown types."""
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:n, 0:n].astype(float)
    cx, cy = n * 0.45, n * 0.55

    elev = 800.0 + 600.0 * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * (n * 0.3) ** 2))
    gy, gx = np.gradient(elev, cellsize)
    slope = np.clip(np.degrees(np.arctan(np.hypot(gx, gy))), 0.0, 45.0)
    aspect = np.degrees(np.arctan2(gx, -gy)) % 360.0          # downslope azimuth

    fuel = np.full((n, n), 10, dtype=int)                     # TL/timber default
    fuel[:, n // 2 - 1:n // 2 + 1] = 91                       # nonburnable river
    fuel[: n // 4, : n // 4] = 1                              # grass corner
    forest = fuel == 10

    cover = np.where(forest, rng.integers(55, 95, (n, n)), 0)
    # CBH (LCP units = m*10): low (crown-prone) south -> high north, + noise.
    cbh = np.clip(15 + 0.9 * yy + rng.normal(0, 8, (n, n)), 5, 140)
    # CBD (LCP units = (kg/m^3)*100): dense on the hill, sparse at the edges.
    cbd = np.clip(8 + 30 * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * (n * 0.35) ** 2))
                  + rng.normal(0, 4, (n, n)), 2, 40)
    ch = np.clip(130 + 0.5 * yy + rng.normal(0, 15, (n, n)), 90, 250)   # height m*10

    z = np.zeros((n, n))
    return pyflam.Landscape(
        fuel_model=fuel, slope=slope.astype(int), aspect=aspect.astype(int),
        elevation=elev.astype(int),
        canopy_cover=np.where(forest, cover, 0).astype(int),
        canopy_base_height=np.where(forest, cbh, z).astype(int),
        canopy_bulk_density=np.where(forest, cbd, z).astype(int),
        canopy_height=np.where(forest, ch, z).astype(int),
        cellsize_x=cellsize, cellsize_y=cellsize, west=0.0, north=n * cellsize,
        slope_units="degrees")


def _crown(ls, model):
    return pyflam.crownfire.crown_fire_potential(
        ls, foliar_moisture=100.0, wind_20ft_ft_per_min=mph_to_ft_per_min(30),
        wind_midflame=pyflam.midflame_field(ls, mph_to_ft_per_min(30)),
        crown_spread=model, **RUN)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="canopy.lcp", help="output .lcp path")
    ap.add_argument("--n", type=int, default=120)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)

    ls = make_synthetic_canopy_landscape(args.n, seed=args.seed)
    ls.to_lcp(args.out)
    back = pyflam.Landscape.from_lcp(args.out)
    print(f"wrote {args.out}  ({back.shape}; bands round-trip: "
          f"CBH {back.canopy_base_height is not None}, "
          f"CBD {back.canopy_bulk_density is not None})\n")

    burnable = back.fuel_model != 91
    cruz = _crown(back, "cruz2005")
    roth = _crown(back, "rothermel1991")
    for name, out in (("Cruz 2005", cruz), ("Rothermel 1991", roth)):
        ft = out["fire_type"][burnable]
        dist = {TYPE_NAMES[t]: int((ft == t).sum()) for t in (0, 1, 2)}
        print(f"{name:<16} crown types (burnable cells): {dist}")

    cmp = validate.compare_categories(
        cruz["fire_type"], roth["fire_type"], labels=[0, 1, 2], mask=burnable)
    print("\nCruz vs Rothermel classification on the synthetic landscape:")
    print("  " + cmp.summary(names=TYPE_NAMES).replace("\n", "\n  "))

    # Exercise the crown-aware field + the plume-coupled crown march end to end.
    caf = pyflam.crown_spread_field(
        back, wind_midflame=pyflam.midflame_field(back, mph_to_ft_per_min(30)),
        wind_direction=45.0, wind_20ft_ft_per_min=mph_to_ft_per_min(30),
        foliar_moisture=100.0, crown_spread="cruz2005", **RUN)
    r, c = np.unravel_index(int(np.argmax(caf.field.ros_max)), caf.field.shape)
    march = pyflam.fire_atmosphere_march(
        back, [(int(r), int(c))], total_time=40, dt=10, speed=8.0, direction=45.0,
        wind_provider=lambda ls_, i, a, s, d: pyflam.pyroconvection._uniform_wind_field(ls_, s, d),
        crown=True, foliar_moisture=100.0, max_wind_factor=4.0, **RUN)
    print(f"\ncrown march burned {int(np.isfinite(march['arrival_time']).sum())} cells; "
          f"final fire-type cells "
          f"{ {TYPE_NAMES[t]: int((march['fire_type']==t).sum()) for t in (0,1,2)} }")


if __name__ == "__main__":
    main()
