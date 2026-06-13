"""Tests for ember spotting and its coupling into the MTT growth engine."""

from __future__ import annotations

import numpy as np
import pytest

import pyflam
from pyflam import spotting
from pyflam.units import mph_to_ft_per_min

SCENARIO = dict(m_1h=0.06, m_10h=0.07, m_100h=0.08,
                m_live_herb=0.60, m_live_woody=0.90)


def _grid(n=60, fuel=102):
    return pyflam.Landscape(
        fuel_model=np.full((n, n), fuel, dtype=int), slope=np.zeros((n, n)),
        aspect=np.zeros((n, n)), cellsize_x=30.0, cellsize_y=30.0,
        west=0.0, north=n * 30.0, slope_units="degrees")


# --- Firebrand physics --------------------------------------------------------

def test_flame_length_from_intensity():
    assert float(spotting.byram_flame_length(0.0)) == 0.0
    # Byram: L = 0.45 I^0.46, monotonic in intensity.
    assert (float(spotting.byram_flame_length(500.0))
            > float(spotting.byram_flame_length(50.0)) > 0.0)


def test_max_spot_distance_scales_with_intensity_and_wind():
    m = pyflam.SpottingModel()
    near = m.max_spot_distance(50.0, mph_to_ft_per_min(10))
    far_intensity = m.max_spot_distance(500.0, mph_to_ft_per_min(10))
    far_wind = m.max_spot_distance(50.0, mph_to_ft_per_min(30))
    assert far_intensity > near          # more energy -> farther
    assert far_wind > near               # more wind -> farther
    assert m.max_spot_distance(50.0, 0.0) == pytest.approx(0.0)  # no wind, no drift


def test_lower_terminal_velocity_flies_farther():
    light = pyflam.SpottingModel(terminal_velocity=100.0)
    heavy = pyflam.SpottingModel(terminal_velocity=400.0)
    u = mph_to_ft_per_min(15)
    assert light.max_spot_distance(300.0, u) > heavy.max_spot_distance(300.0, u)


# --- Spot ignition generation -------------------------------------------------

def test_spots_land_downwind():
    ls = _grid()
    field = pyflam.spread_field(
        ls, wind_midflame=mph_to_ft_per_min(10), wind_direction=270.0, **SCENARIO)
    # A wide burned patch in the west half, wind from the west -> spots to the east.
    arrival = np.full(field.shape, np.inf)
    arrival[20:40, 5:20] = 1.0
    model = pyflam.SpottingModel(spot_probability=1.0, launch_fraction=1.0,
                                 loft_coeff=40.0, terminal_velocity=120.0,
                                 min_intensity=0.0)
    rng = np.random.default_rng(0)
    spots = spotting.generate_spot_ignitions(
        field, arrival, wind_20ft=mph_to_ft_per_min(10), wind_direction=270.0,
        max_time=120.0, model=model, rng=rng)
    assert len(spots) > 0
    src_cols = np.arange(5, 20)
    # Landing columns should be east of (greater than) the source patch.
    land_cols = [c for _, c, _ in spots]
    assert np.median(land_cols) > src_cols.max()


def test_spot_ignition_times_after_delay():
    ls = _grid()
    field = pyflam.spread_field(
        ls, wind_midflame=mph_to_ft_per_min(12), wind_direction=270.0, **SCENARIO)
    arrival = np.full(field.shape, np.inf)
    arrival[30, 10] = 5.0
    model = pyflam.SpottingModel(spot_probability=1.0, launch_fraction=1.0,
                                 spot_delay=3.0, min_intensity=0.0,
                                 loft_coeff=40.0, terminal_velocity=120.0)
    rng = np.random.default_rng(1)
    spots = spotting.generate_spot_ignitions(
        field, arrival, wind_20ft=mph_to_ft_per_min(12), wind_direction=270.0,
        max_time=120.0, model=model, rng=rng)
    for _, _, t in spots:
        assert t == pytest.approx(5.0 + 3.0)   # source arrival + delay


# --- Coupling into MTT growth -------------------------------------------------

def test_spotting_increases_burned_area():
    ls = _grid()
    field = pyflam.spread_field(
        ls, wind_midflame=mph_to_ft_per_min(10), wind_direction=270.0, **SCENARIO)
    model = pyflam.SpottingModel(spot_probability=1.0, launch_fraction=0.3,
                                 loft_coeff=30.0, terminal_velocity=120.0)
    no_spot = pyflam.minimum_travel_time(field, [(30, 10)], max_time=90.0)
    with_spot = pyflam.spread_with_spotting(
        field, [(30, 10)], max_time=90.0, wind_20ft=mph_to_ft_per_min(10),
        wind_direction=270.0, model=model, rng=np.random.default_rng(0))
    assert np.isfinite(with_spot).sum() > np.isfinite(no_spot).sum()


def test_spotting_crosses_a_barrier():
    n = 60
    fuel = np.full((n, n), 102, dtype=int)
    fuel[:, 28:33] = 91                  # a 5-cell-thick nonburnable wall
    ls = pyflam.Landscape(
        fuel_model=fuel, slope=np.zeros((n, n)), aspect=np.zeros((n, n)),
        cellsize_x=30.0, cellsize_y=30.0, west=0.0, north=n * 30.0,
        slope_units="degrees")
    field = pyflam.spread_field(
        ls, wind_midflame=mph_to_ft_per_min(12), wind_direction=270.0, **SCENARIO)
    model = pyflam.SpottingModel(spot_probability=1.0, launch_fraction=0.5,
                                 loft_coeff=60.0, terminal_velocity=80.0)
    no_spot = pyflam.minimum_travel_time(field, [(30, 8)], max_time=150.0)
    with_spot = pyflam.spread_with_spotting(
        field, [(30, 8)], max_time=150.0, wind_20ft=mph_to_ft_per_min(12),
        wind_direction=270.0, model=model, rng=np.random.default_rng(0))
    # The thick wall blocks contiguous spread; spotting jumps it.
    east = np.s_[:, 40:]
    assert np.isfinite(no_spot[east]).sum() == 0
    assert np.isfinite(with_spot[east]).sum() > 0


def test_spotting_is_deterministic_with_seed():
    ls = _grid()
    field = pyflam.spread_field(
        ls, wind_midflame=mph_to_ft_per_min(10), wind_direction=270.0, **SCENARIO)
    model = pyflam.SpottingModel(spot_probability=0.5, launch_fraction=0.2)
    a = pyflam.spread_with_spotting(
        field, [(30, 10)], max_time=90.0, wind_20ft=mph_to_ft_per_min(10),
        wind_direction=270.0, model=model, rng=np.random.default_rng(42))
    b = pyflam.spread_with_spotting(
        field, [(30, 10)], max_time=90.0, wind_20ft=mph_to_ft_per_min(10),
        wind_direction=270.0, model=model, rng=np.random.default_rng(42))
    fin = np.isfinite(a)
    assert np.array_equal(fin, np.isfinite(b))
    assert np.allclose(a[fin], b[fin])


def test_burn_probability_with_spotting_higher():
    ls = _grid(n=50)
    field = pyflam.spread_field(
        ls, wind_midflame=mph_to_ft_per_min(10), wind_direction=270.0, **SCENARIO)
    igns = [(25, 10), (25, 12), (20, 8)]
    model = pyflam.SpottingModel(spot_probability=1.0, launch_fraction=0.3,
                                 loft_coeff=30.0, terminal_velocity=120.0)
    plain, _ = pyflam.burn_probability(field, igns, max_time=60.0)
    spot, _ = pyflam.burn_probability(
        field, igns, max_time=60.0, spotting=model,
        wind_20ft=mph_to_ft_per_min(10), wind_direction=270.0,
        rng=np.random.default_rng(0))
    assert spot.mean() > plain.mean()


# --- Physics-based stochastic model (FirebrandPhysics) ------------------------

def test_terminal_velocity_increases_with_diameter():
    fp = pyflam.FirebrandPhysics()
    wt = fp.terminal_velocity(np.array([0.002, 0.005, 0.01, 0.02]))
    assert np.all(np.diff(wt) > 0.0)
    assert np.all(wt > 0.0)


def test_burnout_time_scales_with_diameter_squared():
    fp = pyflam.FirebrandPhysics()
    assert fp.burnout_time(0.01) == pytest.approx(4.0 * fp.burnout_time(0.005))


def test_loft_height_grows_with_energy_and_shrinks_with_terminal_velocity():
    fp = pyflam.FirebrandPhysics()
    wt = fp.terminal_velocity(0.005)
    assert fp.loft_height(1000.0, wt) > fp.loft_height(100.0, wt)   # more energy
    big_wt = fp.terminal_velocity(0.02)
    assert fp.loft_height(500.0, big_wt) < fp.loft_height(500.0, wt)  # heavier sinks


def test_buoyancy_flux_increases_with_intensity():
    fp = pyflam.FirebrandPhysics()
    assert fp.buoyancy_flux(800.0) > fp.buoyancy_flux(200.0) > 0.0


def test_ignition_probability_falls_with_moisture():
    fp = pyflam.FirebrandPhysics()
    p = fp.ignition_probability(np.array([0.04, 0.12, 0.30]))
    assert np.all(np.diff(p) < 0.0)
    assert np.all((p >= 0.0) & (p <= 1.0))


def _intense_field(n=120, fuel=104, wind_mph=20, load_factor=1.3):
    ls = _grid(n=n, fuel=fuel)
    return pyflam.spread_field(
        ls, wind_midflame=pyflam.midflame_field(ls, mph_to_ft_per_min(wind_mph)),
        wind_direction=270.0, load_factor=load_factor, **SCENARIO)


def test_physics_spots_land_downwind_and_are_stochastic():
    field = _intense_field()
    arrival = np.full(field.shape, np.inf)
    arrival[55:65, 10:30] = 1.0
    fp = pyflam.FirebrandPhysics(launch_fraction=1.0)
    a = fp.generate_spots(field, arrival, wind_20ft=mph_to_ft_per_min(20),
                          wind_direction=270.0, max_time=180.0,
                          rng=np.random.default_rng(1), fuel_moisture=0.06)
    b = fp.generate_spots(field, arrival, wind_20ft=mph_to_ft_per_min(20),
                          wind_direction=270.0, max_time=180.0,
                          rng=np.random.default_rng(2), fuel_moisture=0.06)
    assert len(a) > 0
    assert np.median([c for _, c, _ in a]) > 30      # downwind (east) of source
    assert a != b                                     # stochastic across seeds


def test_physics_higher_intensity_spots_farther():
    # A larger grid so the intense fire's (longer) spots stay in bounds; the
    # source patch sits near the west edge so there is downwind room to the east.
    n = 200
    arrival = np.full((n, n), np.inf)
    arrival[95:105, 15:30] = 1.0
    fp = pyflam.FirebrandPhysics(launch_fraction=1.0)
    weak = fp.generate_spots(
        _intense_field(n=n, wind_mph=8, load_factor=1.0), arrival,
        wind_20ft=mph_to_ft_per_min(8), wind_direction=270.0, max_time=300.0,
        rng=np.random.default_rng(0), fuel_moisture=0.06)
    strong = fp.generate_spots(
        _intense_field(n=n, wind_mph=25, load_factor=1.3), arrival,
        wind_20ft=mph_to_ft_per_min(25), wind_direction=270.0, max_time=300.0,
        rng=np.random.default_rng(0), fuel_moisture=0.06)
    # More intense fire -> embers carried farther downwind (east, higher col).
    assert max(c for _, c, _ in strong) > max(c for _, c, _ in weak)


def test_physics_moist_fuel_ignites_less():
    field = _intense_field()
    arrival = np.full(field.shape, np.inf)
    arrival[55:65, 10:30] = 1.0
    fp = pyflam.FirebrandPhysics(launch_fraction=1.0)
    dry = fp.generate_spots(field, arrival, wind_20ft=mph_to_ft_per_min(20),
                            wind_direction=270.0, max_time=180.0,
                            rng=np.random.default_rng(3), fuel_moisture=0.04)
    wet = fp.generate_spots(field, arrival, wind_20ft=mph_to_ft_per_min(20),
                            wind_direction=270.0, max_time=180.0,
                            rng=np.random.default_rng(3), fuel_moisture=0.30)
    assert len(dry) > len(wet)


def test_physics_model_drives_growth():
    field = _intense_field(n=80)
    fp = pyflam.FirebrandPhysics(launch_fraction=0.3)
    no_spot = pyflam.minimum_travel_time(field, [(40, 12)], max_time=60.0)
    with_spot = pyflam.spread_with_spotting(
        field, [(40, 12)], max_time=60.0, wind_20ft=mph_to_ft_per_min(20),
        wind_direction=270.0, model=fp, rng=np.random.default_rng(0),
        fuel_moisture=0.06)
    assert np.isfinite(with_spot).sum() >= np.isfinite(no_spot).sum()


# --- Calibration against measured spot-distance data --------------------------

def test_spot_distance_distribution_grows_with_intensity_and_wind():
    fp = pyflam.FirebrandPhysics()
    rng = np.random.default_rng(0)
    base = np.percentile(fp.spot_distance_distribution(1000.0, 8.0, rng=rng), 95)
    hi_i = np.percentile(fp.spot_distance_distribution(5000.0, 8.0, rng=rng), 95)
    hi_u = np.percentile(fp.spot_distance_distribution(1000.0, 15.0, rng=rng), 95)
    assert hi_i > base                    # more energy -> farther
    assert hi_u > base                    # more wind -> farther


def test_default_matches_literature_anchors():
    """The calibrated default reproduces literature spot distances within ~2x."""
    rep = spotting.spot_distance_report(rng=np.random.default_rng(0))
    for r in rep:
        assert 0.5 <= r["ratio"] <= 2.0   # order-of-magnitude anchors
    # The three non-extreme anchors should be within ~30%.
    for r in rep[:3]:
        assert 0.7 <= r["ratio"] <= 1.3


def test_calibration_recovers_scaled_data():
    """Calibrating to anchors scaled by 2x should ~double the fitted length."""
    rng = np.random.default_rng(1)
    base_fl = spotting.calibrate_front_length(rng=rng)
    doubled = [(I, U, 2.0 * d) for I, U, d in spotting.LITERATURE_SPOT_ANCHORS]
    scaled_fl = spotting.calibrate_front_length(doubled, rng=rng)
    assert scaled_fl > 1.5 * base_fl      # longer measured distances -> longer L


def test_calibration_is_near_fixed_point_at_default():
    """The default front_length is already ~the fit, so calibration barely moves."""
    fitted = spotting.calibrate_front_length(rng=np.random.default_rng(0))
    default = pyflam.FirebrandPhysics().front_length
    assert 0.8 < fitted / default < 1.25


# --- Particle sub-models vs measured firebrand data ---------------------------

def test_terminal_velocity_matches_measurements():
    """The drag model reproduces measured firebrand settling velocities."""
    fp = pyflam.FirebrandPhysics()
    for d_mm, measured in spotting.FIREBRAND_TERMINAL_VELOCITY_DATA:
        modelled = float(fp.terminal_velocity(d_mm / 1000.0))
        assert modelled == pytest.approx(measured, abs=0.2)


def test_burning_constant_in_tarifa_range():
    fp = pyflam.FirebrandPhysics()
    lo, hi = spotting.TARIFA_BURN_CONSTANT_RANGE
    assert lo <= fp.burning_constant <= hi


def test_burnout_follows_d_squared_law():
    fp = pyflam.FirebrandPhysics()
    # t_burn = d^2 / K, so doubling d quadruples the burnout time.
    assert fp.burnout_time(0.01) == pytest.approx(4.0 * fp.burnout_time(0.005))
    assert fp.burnout_time(0.005) == pytest.approx(0.005 ** 2 / fp.burning_constant)


def test_from_burning_constant_roundtrip():
    fp = pyflam.FirebrandPhysics.from_burning_constant(4.0e-7)
    assert fp.burning_constant == pytest.approx(4.0e-7)


def test_calibrate_burnout_recovers_k():
    # Synthetic measurements from a known K via t = d^2 / K.
    k_true = 6.0e-7
    meas = [(d, d * d / k_true) for d in (0.003, 0.006, 0.012, 0.02)]
    assert spotting.calibrate_burnout(meas) == pytest.approx(k_true, rel=1e-9)
