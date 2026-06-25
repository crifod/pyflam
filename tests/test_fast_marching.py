"""Tests for the anisotropic-Eikonal (fast_marching) propagation backend.

Property tests plus the headline claim: the semi-Lagrangian Eikonal solver has
lower angular (lattice) bias than Dijkstra MTT against the analytic arrival time
on uniform fields (calm circle, wind-driven ellipse).
"""

from __future__ import annotations

import numpy as np
import pytest

import pyflam
from pyflam.units import mph_to_ft_per_min

SCENARIO = dict(m_1h=0.06, m_10h=0.07, m_100h=0.08,
                m_live_herb=0.60, m_live_woody=0.90)


def _grid(n, cellsize=30.0):
    return pyflam.Landscape(
        fuel_model=np.full((n, n), 1, dtype=int),
        slope=np.zeros((n, n)), aspect=np.zeros((n, n)),
        cellsize_x=cellsize, cellsize_y=cellsize, west=0.0, north=n * cellsize,
        slope_units="degrees")


def _analytic(field, src):
    """Exact uniform-field arrival time: distance / R(bearing) (straight geodesics)."""
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


def _mean_rel_error(T, field, src):
    Tex, dist = _analytic(field, src)
    rmax = dist.max()
    m = (dist > 0.4 * rmax) & (dist < 0.7 * rmax) & np.isfinite(T) & (Tex > 0)
    return float((np.abs(T[m] - Tex[m]) / Tex[m]).mean())


# --- basic correctness --------------------------------------------------------

def test_ignition_is_time_zero():
    f = pyflam.spread_field(_grid(21), wind_midflame=0.0, **SCENARIO)
    T = pyflam.anisotropic_eikonal(f, [(10, 10)])
    assert T[10, 10] == 0.0


def test_arrival_increases_with_distance():
    f = pyflam.spread_field(_grid(31), wind_midflame=0.0, **SCENARIO)
    T = pyflam.anisotropic_eikonal(f, [(15, 15)])
    assert T[15, 16] < T[15, 20] < T[15, 25]


def test_no_burnable_ignition_returns_inf():
    ls = _grid(15)
    ls.fuel_model[7, 7] = 91
    f = pyflam.spread_field(ls, wind_midflame=0.0, **SCENARIO)
    T = pyflam.anisotropic_eikonal(f, [(7, 7)])
    assert not np.isfinite(T).any()


def test_barrier_blocks_direct_spread():
    ls = _grid(21)
    ls.fuel_model[:, 11] = 91                       # full nonburnable wall
    f = pyflam.spread_field(ls, wind_midflame=0.0, **SCENARIO)
    T = pyflam.anisotropic_eikonal(f, [(10, 5)])
    assert np.isfinite(T[10, 5])
    assert not np.isfinite(T[10, 15])               # other side is unreachable


def test_max_time_bounds_growth():
    f = pyflam.spread_field(_grid(41), wind_midflame=0.0, **SCENARIO)
    T = pyflam.anisotropic_eikonal(f, [(20, 20)], max_time=20.0)
    assert np.nanmax(np.where(np.isfinite(T), T, -1)) <= 20.0
    assert not np.isfinite(T).all()                 # not the whole grid


# --- dispatch -----------------------------------------------------------------

def test_method_dispatch_matches_direct_call():
    f = pyflam.spread_field(
        _grid(25), wind_midflame=mph_to_ft_per_min(6), wind_direction=270.0, **SCENARIO)
    a = pyflam.minimum_travel_time(f, [(12, 12)], method="fast_marching")
    b = pyflam.anisotropic_eikonal(f, [(12, 12)])
    np.testing.assert_array_equal(a, b)


def test_invalid_method_raises():
    f = pyflam.spread_field(_grid(11), wind_midflame=0.0, **SCENARIO)
    with pytest.raises(ValueError):
        pyflam.minimum_travel_time(f, [(5, 5)], method="nope")


@pytest.mark.skipif(not pyflam.mtt._HAVE_NUMBA, reason="needs numba")
def test_all_backends_agree():
    """heap (FMM+bbox), numba sweep and numpy sweep must give the same field."""
    f = pyflam.spread_field(
        _grid(31), wind_midflame=mph_to_ft_per_min(8), wind_direction=300.0, **SCENARIO)
    a = pyflam.anisotropic_eikonal(f, [(15, 15)], backend="numpy")
    for backend in ("numba", "heap"):
        b = pyflam.anisotropic_eikonal(f, [(15, 15)], backend=backend)
        both = np.isfinite(a) & np.isfinite(b)
        assert np.max(np.abs(a[both] - b[both])) < 1e-9, backend
        assert np.array_equal(np.isfinite(a), np.isfinite(b)), backend


@pytest.mark.skipif(not pyflam.mtt._HAVE_NUMBA, reason="needs numba")
def test_heap_backend_prunes_with_max_time():
    """The heap backend must match the sweep on a bounded fire (and only touch it)."""
    f = pyflam.spread_field(
        _grid(121), wind_midflame=mph_to_ft_per_min(8), wind_direction=270.0, **SCENARIO)
    heap = pyflam.anisotropic_eikonal(f, [(60, 60)], max_time=8.0, backend="heap")
    sweep = pyflam.anisotropic_eikonal(f, [(60, 60)], max_time=8.0, backend="numba")
    both = np.isfinite(heap) & np.isfinite(sweep)
    assert np.max(np.abs(heap[both] - sweep[both])) < 1e-9
    assert np.array_equal(np.isfinite(heap), np.isfinite(sweep))
    assert not np.isfinite(heap).all()              # bounded: not the whole grid


def test_spread_perimeter_accepts_method():
    out = pyflam.spread_perimeter(
        _grid(21), [(10, 10)], wind_midflame=0.0, method="fast_marching", **SCENARIO)
    assert out["arrival_time"][10, 10] == 0.0


# --- the headline: lower lattice bias than MTT --------------------------------

def test_fast_marching_beats_mtt_on_calm_circle():
    f = pyflam.spread_field(_grid(41), wind_midflame=0.0, **SCENARIO)
    e_mtt = _mean_rel_error(pyflam.minimum_travel_time(f, [(20, 20)], ring=2), f, (20, 20))
    e_fm = _mean_rel_error(
        pyflam.minimum_travel_time(f, [(20, 20)], method="fast_marching"), f, (20, 20))
    assert e_fm < e_mtt


def test_fast_marching_beats_mtt_on_wind_ellipse():
    f = pyflam.spread_field(
        _grid(41), wind_midflame=mph_to_ft_per_min(10), wind_direction=270.0, **SCENARIO)
    e_mtt = _mean_rel_error(pyflam.minimum_travel_time(f, [(20, 20)], ring=2), f, (20, 20))
    e_fm = _mean_rel_error(
        pyflam.minimum_travel_time(f, [(20, 20)], method="fast_marching"), f, (20, 20))
    assert e_fm < e_mtt


def test_calm_fire_is_more_circular_than_mtt():
    """Std of arrival time around a ring measures faceting; FM should be smoother."""
    f = pyflam.spread_field(_grid(41), wind_midflame=0.0, **SCENARIO)
    nr, nc = f.shape
    rr, cc = np.mgrid[0:nr, 0:nc]
    dist = np.hypot((cc - 20) * 30.0, (rr - 20) * 30.0)
    ring = (dist > 9 * 30.0) & (dist < 11 * 30.0)

    def spread_std(T):
        return float(np.std(T[ring] / dist[ring]))     # normalise out radius

    s_mtt = spread_std(pyflam.minimum_travel_time(f, [(20, 20)], ring=2))
    s_fm = spread_std(pyflam.minimum_travel_time(f, [(20, 20)], method="fast_marching"))
    assert s_fm < s_mtt
