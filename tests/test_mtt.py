"""Tests for the directional spread + Minimum Travel Time engine.

Physics/property tests (the fire must spread faster downwind, arrival time must
grow with distance, barriers must block, a no-wind fire must be radially
symmetric) plus a small golden-master lock.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import pyflam
from pyflam import mtt
from pyflam.units import mph_to_ft_per_min

SCENARIO = dict(m_1h=0.06, m_10h=0.07, m_100h=0.08,
                m_live_herb=0.60, m_live_woody=0.90)


def _grid(n, fuel=1, slope=0.0, aspect=0.0, cellsize=30.0):
    return pyflam.Landscape(
        fuel_model=np.full((n, n), fuel, dtype=int),
        slope=np.full((n, n), slope, dtype=float),
        aspect=np.full((n, n), aspect, dtype=float),
        cellsize_x=cellsize, cellsize_y=cellsize, west=0.0, north=n * cellsize,
        slope_units="degrees",
    )


# --- Spread field (directional elliptical spread, step 3) ---------------------

def test_ros_max_matches_surface_max_spread():
    """With wind only (no slope), heading ROS equals the Rothermel max spread."""
    ls = _grid(5)
    field = pyflam.spread_field(
        ls, wind_midflame=mph_to_ft_per_min(5), wind_direction=270.0, **SCENARIO)
    expected = pyflam.spread(
        pyflam.get_fuel_model(1), wind_midflame=mph_to_ft_per_min(5),
        slope=0.0, **SCENARIO).rate_of_spread
    assert field.ros_max[2, 2] == pytest.approx(expected, rel=1e-9)


def test_heading_points_downwind():
    """Wind from the west (270 deg FROM) drives spread toward the east (90)."""
    ls = _grid(5)
    field = pyflam.spread_field(
        ls, wind_midflame=mph_to_ft_per_min(8), wind_direction=270.0, **SCENARIO)
    assert field.heading[2, 2] == pytest.approx(90.0, abs=1e-6)


def test_calm_fire_is_a_circle():
    ls = _grid(5)
    field = pyflam.spread_field(ls, wind_midflame=0.0, **SCENARIO)
    assert np.allclose(field.eccentricity, 0.0)


def test_directional_ros_max_at_heading():
    ls = _grid(5)
    field = pyflam.spread_field(
        ls, wind_midflame=mph_to_ft_per_min(10), wind_direction=270.0, **SCENARIO)
    head = field.heading[2, 2]
    r_head = float(field.directional_ros(head)[2, 2])
    r_back = float(field.directional_ros((head + 180) % 360)[2, 2])
    r_flank = float(field.directional_ros((head + 90) % 360)[2, 2])
    assert r_head == pytest.approx(field.ros_max[2, 2], rel=1e-9)
    assert r_head > r_flank > r_back > 0.0


# --- Minimum Travel Time ------------------------------------------------------

def test_ignition_is_time_zero():
    ls = _grid(21)
    at = pyflam.minimum_travel_time(
        pyflam.spread_field(ls, wind_midflame=mph_to_ft_per_min(5), **SCENARIO),
        [(10, 10)])
    assert at[10, 10] == 0.0


def test_arrival_time_increases_with_distance():
    ls = _grid(31)
    field = pyflam.spread_field(ls, wind_midflame=0.0, **SCENARIO)
    at = pyflam.minimum_travel_time(field, [(15, 15)])
    assert at[15, 16] < at[15, 20] < at[15, 25]


def test_downwind_faster_than_upwind():
    ls = _grid(41)
    field = pyflam.spread_field(
        ls, wind_midflame=mph_to_ft_per_min(6), wind_direction=270.0, **SCENARIO)
    at = pyflam.minimum_travel_time(field, [(20, 20)])
    east = at[20, 35]   # downwind
    west = at[20, 5]    # upwind
    assert east < west


def test_calm_fire_radially_symmetric():
    ls = _grid(31)
    field = pyflam.spread_field(ls, wind_midflame=0.0, **SCENARIO)
    at = pyflam.minimum_travel_time(field, [(15, 15)])
    # Equal-distance cells along the axes should arrive at the same time.
    assert at[15, 25] == pytest.approx(at[15, 5], rel=1e-9)
    assert at[15, 25] == pytest.approx(at[5, 15], rel=1e-9)
    assert at[25, 15] == pytest.approx(at[5, 15], rel=1e-9)


def test_calm_arrival_matches_distance_over_ros():
    """Along a primitive lattice axis, arrival ~ distance / ROS (no wind)."""
    ls = _grid(31, cellsize=30.0)
    field = pyflam.spread_field(ls, wind_midflame=0.0, **SCENARIO)
    ros = field.ros_max[15, 15]            # ft/min, isotropic
    at = pyflam.minimum_travel_time(field, [(15, 15)])
    dist_ft = 10 * 30.0                    # 10 cells east
    assert at[15, 25] == pytest.approx(dist_ft / ros, rel=1e-6)


def test_nonburnable_barrier_blocks_spread():
    ls = _grid(21, fuel=1)
    ls.fuel_model[:, 10] = 91              # NB1 wall down the middle column
    ls.fuel_model[0, 10] = 1               # ...with no gap except the top edge
    field = pyflam.spread_field(ls, wind_midflame=0.0, **SCENARIO)
    at = pyflam.minimum_travel_time(field, [(20, 5)])
    # Ignite bottom-left; the only way across is around the top -> the far side
    # is reachable but takes much longer than the same distance with no wall.
    open_field = pyflam.spread_field(_grid(21), wind_midflame=0.0, **SCENARIO)
    at_open = pyflam.minimum_travel_time(open_field, [(20, 5)])
    assert at[20, 15] > at_open[20, 15]
    assert math.isinf(at[5, 10])           # the barrier cell itself never burns


def test_perimeter_and_area_monotonic_in_time():
    ls = _grid(41)
    field = pyflam.spread_field(
        ls, wind_midflame=mph_to_ft_per_min(5), wind_direction=270.0, **SCENARIO)
    at = pyflam.minimum_travel_time(field, [(20, 20)])
    a10 = mtt.burned_area(at, 10.0, 30.0, 30.0)
    a30 = mtt.burned_area(at, 30.0, 30.0, 30.0)
    assert a30 >= a10 > 0.0
    assert pyflam.perimeter_mask(at, 30.0).sum() >= pyflam.perimeter_mask(at, 10.0).sum()


def test_ignition_from_xy_roundtrips():
    ls = _grid(21, cellsize=30.0)
    # cell (5, 8) center in world coords
    x = ls.west + (8 + 0.5) * ls.cellsize_x
    y = ls.north - (5 + 0.5) * ls.cellsize_y
    assert pyflam.ignition_from_xy(ls, x, y) == (5, 8)


def test_max_time_limits_growth():
    ls = _grid(41)
    field = pyflam.spread_field(ls, wind_midflame=0.0, **SCENARIO)
    full = pyflam.minimum_travel_time(field, [(20, 20)])
    capped = pyflam.minimum_travel_time(field, [(20, 20)], max_time=5.0)
    assert np.isfinite(capped).sum() < np.isfinite(full).sum()
    assert np.nanmax(capped[np.isfinite(capped)]) <= 5.0 + 1e-9


# --- Scalable backend equivalence & size --------------------------------------

def test_scipy_matches_python_reference():
    """The vectorized/SciPy solver must match the pure-Python heap Dijkstra."""
    ls = _grid(35, slope=10.0, aspect=180.0)
    field = pyflam.spread_field(
        ls, wind_midflame=mph_to_ft_per_min(7), wind_direction=225.0, **SCENARIO)
    fast = pyflam.minimum_travel_time(field, [(17, 17)])
    ref = mtt._mtt_python(field, [(17, 17)])
    both = np.isfinite(fast) & np.isfinite(ref)
    assert np.allclose(fast[both], ref[both], rtol=1e-9, atol=1e-9)
    assert np.array_equal(np.isfinite(fast), np.isfinite(ref))


def test_scipy_matches_python_with_barrier_and_multi_ignition():
    ls = _grid(25)
    ls.fuel_model[5:20, 12] = 91          # nonburnable wall with gaps at ends
    field = pyflam.spread_field(
        ls, wind_midflame=mph_to_ft_per_min(4), wind_direction=270.0, **SCENARIO)
    igns = [(3, 3), (22, 4)]
    fast = pyflam.minimum_travel_time(field, igns)
    ref = mtt._mtt_python(field, igns)
    both = np.isfinite(fast) & np.isfinite(ref)
    assert np.allclose(fast[both], ref[both], rtol=1e-9, atol=1e-9)
    assert np.array_equal(np.isfinite(fast), np.isfinite(ref))


def test_large_grid_runs():
    """A quarter-million-cell grid solves quickly via the sparse-graph backend."""
    n = 500
    ls = _grid(n, cellsize=10.0)          # meter-ish scale, 250k cells
    field = pyflam.spread_field(
        ls, wind_midflame=mph_to_ft_per_min(6), wind_direction=270.0, **SCENARIO)
    at = pyflam.minimum_travel_time(field, [(n // 2, n // 2)])
    assert at[n // 2, n // 2] == 0.0
    assert np.isfinite(at).all()          # homogeneous burnable grid -> all reached
    # Downwind (east) still arrives before upwind (west).
    assert at[n // 2, n // 2 + 200] < at[n // 2, n // 2 - 200]


def test_build_traveltime_graph_shape():
    ls = _grid(10)
    field = pyflam.spread_field(ls, wind_midflame=0.0, **SCENARIO)
    g = mtt.build_traveltime_graph(field)
    assert g.shape == (100, 100)
    assert g.nnz > 0
    assert (g.data > 0).all()             # travel times are strictly positive


def test_chunked_graph_identical_to_single_shot():
    """The chunked build must produce the same graph as the single-shot build."""
    ls = _grid(37, slope=12.0, aspect=135.0)
    ls.fuel_model[10:25, 18] = 91         # a barrier, to exercise dropped edges
    field = pyflam.spread_field(
        ls, wind_midflame=mph_to_ft_per_min(6), wind_direction=300.0, **SCENARIO)
    g_single = mtt.build_traveltime_graph(field, chunk_rows=0).tocsr()
    for ch in (1, 5, 8, 100):
        g_chunk = mtt.build_traveltime_graph(field, chunk_rows=ch).tocsr()
        assert g_chunk.shape == g_single.shape
        assert g_chunk.nnz == g_single.nnz
        # Same edges and weights regardless of how rows are banded.
        d = (g_chunk - g_single)
        assert abs(d).max() < 1e-12 if d.nnz else True


def test_chunked_mtt_matches_single_shot():
    ls = _grid(40)
    ls.fuel_model[:, 20] = 91
    ls.fuel_model[0, 20] = 1              # leave a gap so the far side is reachable
    field = pyflam.spread_field(
        ls, wind_midflame=mph_to_ft_per_min(5), wind_direction=270.0, **SCENARIO)
    igns = [(20, 5)]
    a_single = pyflam.minimum_travel_time(field, igns, chunk_rows=0)
    a_chunk = pyflam.minimum_travel_time(field, igns, chunk_rows=7)
    both = np.isfinite(a_single) & np.isfinite(a_chunk)
    assert np.array_equal(np.isfinite(a_single), np.isfinite(a_chunk))
    assert np.allclose(a_single[both], a_chunk[both], rtol=1e-12, atol=1e-12)


# --- Burn probability ---------------------------------------------------------

def test_burn_probability_basic():
    ls = _grid(31)
    field = pyflam.spread_field(
        ls, wind_midflame=mph_to_ft_per_min(5), wind_direction=270.0, **SCENARIO)
    # Three fires from the same cell, bounded to 20 min: every fire footprint is
    # identical, so burn prob is exactly 0 or 1 and the ignition is always 1.
    prob, nf = pyflam.burn_probability(field, [(15, 15)] * 3, max_time=20.0)
    assert nf == 3
    assert prob[15, 15] == 1.0
    assert set(np.unique(prob)).issubset({0.0, 1.0})
    # Bounded growth -> not the whole grid burns.
    assert 0.0 < prob.mean() < 1.0


def test_burn_probability_accumulates_across_ignitions():
    ls = _grid(41)
    field = pyflam.spread_field(ls, wind_midflame=0.0, **SCENARIO)
    igns = [(10, 10), (30, 30)]            # two far-apart fires
    prob, nf = pyflam.burn_probability(field, igns, max_time=15.0)
    assert nf == 2
    # Each ignition burns its own neighborhood with prob 1/2 unless they overlap.
    assert prob[10, 10] == pytest.approx(0.5)
    assert prob[30, 30] == pytest.approx(0.5)


def test_burn_probability_skips_nonburnable_ignition():
    ls = _grid(21)
    ls.fuel_model[10, 10] = 91            # nonburnable ignition
    field = pyflam.spread_field(ls, wind_midflame=0.0, **SCENARIO)
    prob, nf = pyflam.burn_probability(field, [(10, 10), (5, 5)], max_time=10.0)
    assert nf == 1                        # only the burnable ignition counts


# --- Golden-master ------------------------------------------------------------

def test_golden_master_arrival():
    """Lock the no-wind FM1 arrival time 10 cells east of the ignition."""
    ls = _grid(31, cellsize=30.0)
    field = pyflam.spread_field(ls, wind_midflame=0.0, **SCENARIO)
    at = pyflam.minimum_travel_time(field, [(15, 15)])
    assert at[15, 25] == pytest.approx(65.15, abs=0.01)
