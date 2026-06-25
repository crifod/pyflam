#!/usr/bin/env python
"""Lattice-bias benchmark: MTT (Dijkstra) vs the fast_marching (Eikonal) backend.

On a *uniform* spread field the arrival time has the closed form
``T(x) = |x| / R(bearing(x))`` (straight geodesics), so we can measure each
solver's angular (lattice) discretization error exactly. The semi-Lagrangian
anisotropic-Eikonal backend should carry markedly less bias than Dijkstra-on-a-
lattice at the same resolution, especially for wind-driven (elliptical) spread.

Usage (from the pyflam dir):
    PYTHONPATH=src python tests/benchmark_propagation.py
    PYTHONPATH=src python tests/benchmark_propagation.py --n 61 --wind 12
"""

from __future__ import annotations

import argparse
import time

import numpy as np

import pyflam
from pyflam.units import mph_to_ft_per_min

SCENARIO = dict(m_1h=0.06, m_10h=0.07, m_100h=0.08,
                m_live_herb=0.60, m_live_woody=0.90)


def _grid(n, cell=30.0):
    return pyflam.Landscape(
        fuel_model=np.full((n, n), 1, dtype=int),
        slope=np.zeros((n, n)), aspect=np.zeros((n, n)),
        cellsize_x=cell, cellsize_y=cell, west=0.0, north=n * cell,
        slope_units="degrees")


def _analytic(field, src):
    nr, nc = field.shape
    rr, cc = np.mgrid[0:nr, 0:nc]
    dx = (cc - src[1]) * field.cellsize_x
    dy = -(rr - src[0]) * field.cellsize_y
    dist = np.hypot(dx, dy)
    bearing = np.degrees(np.arctan2(dx, dy)) % 360.0
    with np.errstate(divide="ignore", invalid="ignore"):
        t = dist / field.directional_ros(bearing)
    t[src] = 0.0
    return t, dist


def _errors(T, field, src):
    Tex, dist = _analytic(field, src)
    rmax = dist.max()
    m = (dist > 0.4 * rmax) & (dist < 0.7 * rmax) & np.isfinite(T) & (Tex > 0)
    rel = np.abs(T[m] - Tex[m]) / Tex[m]
    return 100 * rel.mean(), 100 * rel.max()


def _run(label, field, src):
    print(f"\n{label}")
    print(f"  {'method':<22}{'mean err %':>11}{'max err %':>11}{'time (s)':>10}")
    for name, kw in (("MTT (Dijkstra, ring 2)", dict(method="mtt", ring=2)),
                     ("MTT (Dijkstra, ring 3)", dict(method="mtt", ring=3)),
                     ("fast_marching (Eikonal)", dict(method="fast_marching"))):
        t0 = time.time()
        T = pyflam.minimum_travel_time(field, [src], **kw)
        dt = time.time() - t0
        mean_e, max_e = _errors(T, field, src)
        print(f"  {name:<22}{mean_e:>11.2f}{max_e:>11.2f}{dt:>10.3f}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=41, help="grid size (cells)")
    ap.add_argument("--wind", type=float, default=10.0, help="midflame wind (mph)")
    ap.add_argument("--wind-dir", type=float, default=270.0)
    args = ap.parse_args(argv)

    src = (args.n // 2, args.n // 2)
    calm = pyflam.spread_field(_grid(args.n), wind_midflame=0.0, **SCENARIO)
    windy = pyflam.spread_field(
        _grid(args.n), wind_midflame=mph_to_ft_per_min(args.wind),
        wind_direction=args.wind_dir, **SCENARIO)

    print(f"propagation bias benchmark  (grid {args.n}x{args.n}, analytic ground truth)")
    _run("CALM  (true shape: circle)", calm, src)
    _run(f"WIND {args.wind:.0f} mph  (true shape: ellipse)", windy, src)

    # Bounded-fire timing: where the heap backend's max_time pruning pays off.
    big_n, mt = 400, 30.0
    big = pyflam.spread_field(
        _grid(big_n), wind_midflame=mph_to_ft_per_min(args.wind),
        wind_direction=args.wind_dir, **SCENARIO)
    bsrc = (big_n // 2, big_n // 2)
    print(f"\nBOUNDED fire timing  (grid {big_n}x{big_n}, {mt:.0f} min — heap prunes)")
    print(f"  {'method':<30}{'time (s)':>10}{'burned cells':>14}")
    for name, fn in (
        ("MTT (Dijkstra)",
         lambda: pyflam.minimum_travel_time(big, [bsrc], max_time=mt, method="mtt")),
        ("fast_marching heap (FMM+bbox)",
         lambda: pyflam.anisotropic_eikonal(big, [bsrc], max_time=mt, backend="heap")),
        ("fast_marching full sweep",
         lambda: pyflam.anisotropic_eikonal(big, [bsrc], max_time=mt, backend="numba")),
    ):
        t0 = time.time(); T = fn(); dt = time.time() - t0
        print(f"  {name:<30}{dt:>10.3f}{int(np.isfinite(T).sum()):>14}")


if __name__ == "__main__":
    main()
