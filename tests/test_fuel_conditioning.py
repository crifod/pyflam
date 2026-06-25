"""Tests for per-cell dead fuel moisture conditioning (terrain + canopy)."""

from __future__ import annotations

from datetime import datetime

import numpy as np
import pytest

import pyflam
from pyflam import fuel_conditioning as fc
from pyflam.atmosphere import AtmosphericState, ConstantAtmosphere


def _terrain(n=20, *, slope=20.0, aspect=180.0, cover=0.0, cellsize=30.0):
    """A uniform-aspect slope with optional canopy cover."""
    return pyflam.Landscape(
        fuel_model=np.full((n, n), 1, dtype=int),
        slope=np.full((n, n), slope, dtype=float),
        aspect=np.full((n, n), aspect, dtype=float),
        canopy_cover=np.full((n, n), cover, dtype=float),
        elevation=np.zeros((n, n), dtype=float),
        cellsize_x=cellsize, cellsize_y=cellsize, west=0.0, north=n * cellsize,
        slope_units="degrees",
    )


# --- solar geometry -----------------------------------------------------------

def test_solar_noon_points_south_northern_hemisphere():
    zen, az = fc.solar_position(latitude=45.0, day_of_year=172, hour=12.0)
    assert az == pytest.approx(180.0, abs=1.0)        # sun due south at solar noon
    assert zen == pytest.approx(45.0 - 23.45, abs=2.0)  # zenith ~ lat - declination


def test_solar_night_is_below_horizon():
    zen, _ = fc.solar_position(latitude=45.0, day_of_year=172, hour=0.0)
    assert zen >= 90.0


def test_equation_of_time_in_expected_range():
    eot = [fc.equation_of_time(d) for d in range(1, 366)]
    assert -16.0 < min(eot) < -12.0          # early-November-ish minimum ~ -14 min
    assert 14.0 < max(eot) < 18.0            # early-November maximum ~ +16 min


def test_clock_time_correction_shifts_solar_noon():
    """A site west of its DST time-zone meridian sees solar noon after clock noon.

    Tuscany ~11 deg E on CEST (UTC+2, zone meridian 30 deg E) is well west of its
    meridian, so true solar noon falls ~1.3 h after 12:00 clock time (~13:17).
    """
    lat, doy = 43.0, 172
    zen_clock_noon, _ = fc.solar_position(
        lat, doy, 12.0, longitude=11.0, timezone=2.0)
    zeniths = {h: fc.solar_position(lat, doy, h, longitude=11.0, timezone=2.0)[0]
               for h in np.arange(12.0, 14.01, 0.25)}
    h_min = min(zeniths, key=zeniths.get)
    assert 13.0 <= h_min <= 13.5              # minimum zenith ~1.3 h after clock noon
    # And passing longitude/timezone actually changes the result vs solar-time mode.
    assert zen_clock_noon != fc.solar_position(lat, doy, 12.0)[0]


# --- terrain insolation -------------------------------------------------------

def test_flat_insolation_is_cos_zenith():
    ls = _terrain(slope=0.0)
    zen, az = 30.0, 180.0
    fac = fc.terrain_insolation_factor(ls, zen, az)
    assert np.allclose(fac, np.cos(np.radians(zen)))


def test_south_slope_brighter_than_north_at_solar_noon():
    south = _terrain(aspect=180.0)
    north = _terrain(aspect=0.0)
    zen, az = fc.solar_position(latitude=45.0, day_of_year=172, hour=12.0)
    fs = fc.terrain_insolation_factor(south, zen, az).mean()
    fn = fc.terrain_insolation_factor(north, zen, az).mean()
    assert fs > fn


def test_night_insolation_is_zero():
    ls = _terrain()
    fac = fc.terrain_insolation_factor(ls, zenith_deg=95.0, azimuth_deg=180.0)
    assert np.all(fac == 0.0)


# --- canopy -------------------------------------------------------------------

def test_canopy_transmission_endpoints():
    assert np.allclose(fc.canopy_transmission(_terrain(cover=0.0)), 1.0)
    assert np.allclose(fc.canopy_transmission(_terrain(cover=100.0)), 0.0)
    assert np.allclose(fc.canopy_transmission(_terrain(cover=50.0)), 0.5)


# --- VPD model ----------------------------------------------------------------

def test_vpd_zero_at_saturation():
    assert fc.vapour_pressure_deficit(20.0, 100.0) == pytest.approx(0.0)
    assert fc.vapour_pressure_deficit(20.0, 0.0) > 0.0    # dry air -> positive VPD


def test_vpd_moisture_decreases_with_drier_air():
    moist_air = fc.dead_fuel_moisture_vpd(30.0, 60.0)
    dry_air = fc.dead_fuel_moisture_vpd(30.0, 10.0)
    assert dry_air < moist_air                            # drier air -> drier fuel


def test_vpd_moisture_matches_nolan_coefficients():
    a, b, c = fc.NOLAN_VPD_COEFFS
    vpd = fc.vapour_pressure_deficit(30.0, 25.0)
    assert fc.dead_fuel_moisture_vpd(30.0, 25.0) == pytest.approx(
        a + b * np.exp(-c * vpd))


# --- full conditioning --------------------------------------------------------

WX = dict(temperature=30.0, relative_humidity=20.0)


def test_zero_insolation_recovers_plain_emc():
    """With no sun the conditioned 1-h moisture is just the ambient NFDRS EMC."""
    ls = _terrain()
    out = fc.condition_dead_fuel_moisture(ls, insolation=np.zeros(ls.shape), **WX)
    emc = pyflam.atmosphere.equilibrium_moisture_content(30.0, 20.0) / 100.0
    assert np.allclose(out["m_1h"], emc)


def test_sun_dries_fuel_below_ambient():
    ls = _terrain()
    sun = fc.condition_dead_fuel_moisture(
        ls, insolation=np.ones(ls.shape), **WX)["m_1h"]
    shade = fc.condition_dead_fuel_moisture(
        ls, insolation=np.zeros(ls.shape), **WX)["m_1h"]
    assert np.all(sun < shade)                        # sun-exposed fuel is drier


def test_south_facing_drier_than_north_at_noon():
    common = dict(latitude=45.0, day_of_year=172, hour=12.0, **WX)
    south = fc.condition_dead_fuel_moisture(_terrain(aspect=180.0), **common)["m_1h"]
    north = fc.condition_dead_fuel_moisture(_terrain(aspect=0.0), **common)["m_1h"]
    assert south.mean() < north.mean()


def test_canopy_keeps_fuel_moister():
    common = dict(latitude=45.0, day_of_year=172, hour=12.0, **WX)
    open_ = fc.condition_dead_fuel_moisture(_terrain(cover=0.0), **common)["m_1h"]
    closed = fc.condition_dead_fuel_moisture(_terrain(cover=90.0), **common)["m_1h"]
    assert closed.mean() > open_.mean()


def test_vpd_model_conditions_per_cell_and_differs_from_emc():
    common = dict(latitude=45.0, day_of_year=172, hour=12.0, **WX)
    emc = fc.condition_dead_fuel_moisture(_terrain(aspect=180.0), **common)
    vpd_south = fc.condition_dead_fuel_moisture(
        _terrain(aspect=180.0), model="vpd", **common)["m_1h"]
    vpd_north = fc.condition_dead_fuel_moisture(
        _terrain(aspect=0.0), model="vpd", **common)["m_1h"]
    # The VPD submodel still responds to terrain (south drier than north) ...
    assert vpd_south.mean() < vpd_north.mean()
    # ... and gives a different field than the EMC submodel.
    assert not np.allclose(vpd_south, emc["m_1h"])


def test_bad_model_raises():
    with pytest.raises(ValueError):
        fc.condition_dead_fuel_moisture(
            _terrain(), insolation=np.zeros(_terrain().shape), model="nope", **WX)


def test_class_offsets_order():
    ls = _terrain()
    out = fc.condition_dead_fuel_moisture(
        ls, latitude=45.0, day_of_year=172, hour=12.0, **WX)
    assert np.all(out["m_10h"] > out["m_1h"])
    assert np.all(out["m_100h"] > out["m_10h"])


def test_elevation_lapse_cools_and_moistens_high_cells():
    n = 10
    elev = np.tile(np.linspace(0.0, 1000.0, n), (n, 1))   # rises west->east
    ls = pyflam.Landscape(
        fuel_model=np.full((n, n), 1, dtype=int),
        slope=np.zeros((n, n)), aspect=np.zeros((n, n)),
        elevation=elev, canopy_cover=np.zeros((n, n)),
        cellsize_x=30.0, cellsize_y=30.0, west=0.0, north=n * 30.0,
        slope_units="degrees")
    out = fc.condition_dead_fuel_moisture(
        ls, insolation=np.zeros(ls.shape), reference_elevation=0.0, **WX)
    # Cooler high cells (constant vapour pressure) -> higher RH -> moister fuel.
    assert out["m_1h"][0, -1] > out["m_1h"][0, 0]


# --- weather-driven setup (condition_from_weather) ----------------------------

def test_condition_from_weather_manual_matches_explicit():
    """The manual fallback equals calling condition_dead_fuel_moisture directly."""
    ls = _terrain()
    when = datetime(2025, 7, 29, 14, 0)
    a = fc.condition_from_weather(
        ls, time=when, temperature=30.0, relative_humidity=20.0, latitude=43.0)
    b = fc.condition_dead_fuel_moisture(
        ls, temperature=30.0, relative_humidity=20.0,
        latitude=43.0, day_of_year=when.timetuple().tm_yday, hour=14.0)
    np.testing.assert_allclose(a["m_1h"], b["m_1h"])


def test_condition_from_weather_accepts_string_time():
    ls = _terrain()
    out = fc.condition_from_weather(
        ls, time="2025-07-29 14:00", temperature=30.0, relative_humidity=20.0,
        latitude=43.0)
    assert out["m_1h"].shape == ls.shape


def test_condition_from_weather_pulls_temperature_and_rh_from_provider():
    """A meteo provider supplies T/RH; manual values are not needed."""
    ls = _terrain()
    state = AtmosphericState(wind_speed=0.0, wind_direction=0.0,
                             temperature=30.0, relative_humidity=20.0)
    provider = ConstantAtmosphere(state)
    from_provider = fc.condition_from_weather(
        ls, time=datetime(2025, 7, 29, 14, 0), atmosphere=provider, latitude=43.0,
        timezone=None)
    # Same as feeding the provider's T/RH in by hand (UTC time, tz=0 default).
    by_hand = fc.condition_from_weather(
        ls, time=datetime(2025, 7, 29, 14, 0), temperature=30.0,
        relative_humidity=20.0, latitude=43.0, timezone=0.0)
    np.testing.assert_allclose(from_provider["m_1h"], by_hand["m_1h"])


def test_condition_from_weather_provider_overrides_and_drives_moisture():
    """Hotter/drier provider weather yields drier conditioned fuel."""
    ls = _terrain()
    dry = ConstantAtmosphere(AtmosphericState(
        wind_speed=0.0, wind_direction=0.0, temperature=35.0, relative_humidity=10.0))
    humid = ConstantAtmosphere(AtmosphericState(
        wind_speed=0.0, wind_direction=0.0, temperature=20.0, relative_humidity=60.0))
    when = datetime(2025, 7, 29, 14, 0)
    m_dry = fc.condition_from_weather(ls, time=when, atmosphere=dry, latitude=43.0)
    m_humid = fc.condition_from_weather(ls, time=when, atmosphere=humid, latitude=43.0)
    assert m_dry["m_1h"].mean() < m_humid["m_1h"].mean()


def test_condition_from_weather_requires_weather_or_provider():
    with pytest.raises(ValueError):
        fc.condition_from_weather(_terrain(), time="2025-07-29 14:00", latitude=43.0)


# --- integration with the spread engine --------------------------------------

def test_conditioned_moisture_feeds_spread_field():
    """South vs north aspects must yield different fuel moisture and spread rate."""
    n = 12
    aspect = np.full((n, n), 0.0)
    aspect[:, n // 2:] = 180.0                      # east half faces south
    ls = pyflam.Landscape(
        fuel_model=np.full((n, n), 1, dtype=int),
        slope=np.full((n, n), 20.0), aspect=aspect,
        canopy_cover=np.zeros((n, n)), elevation=np.zeros((n, n)),
        cellsize_x=30.0, cellsize_y=30.0, west=0.0, north=n * 30.0,
        slope_units="degrees")
    moist = fc.condition_dead_fuel_moisture(
        ls, latitude=45.0, day_of_year=172, hour=12.0, **WX)
    # South half is drier than north half.
    assert moist["m_1h"][:, n // 2:].mean() < moist["m_1h"][:, :n // 2].mean()

    field = pyflam.spread_field(
        ls, wind_midflame=0.0, m_live_herb=0.6, m_live_woody=0.9, **moist)
    assert field.ros_max.shape == ls.shape
    # Drier south fuel spreads at least as fast as the moister north fuel.
    assert field.ros_max[:, n // 2:].mean() >= field.ros_max[:, :n // 2].mean()
