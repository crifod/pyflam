"""Tests for the fire-weather-history blend (pure logic).

The ERA5 fetch and the top-level fire_weather_spinup orchestration touch the
network/CDS and are not exercised here; the blend and its effect on the spun-up
FWI codes are.
"""

from __future__ import annotations

import datetime as dt

import pytest

from pyflam import fire_weather_history as fwh, fwi


def _weather(n=21, *, temp=32.0, rh=25.0, wind_ms=4.0, rain_mm=0.0, month=6):
    """A run of identical dry-day ERA5-style records ending 2026-07-02."""
    end = dt.date(2026, 7, 2)
    return [dict(date=end - dt.timedelta(days=n - 1 - i), temperature=temp,
                 relative_humidity=rh, wind_ms=wind_ms, rain_mm=rain_mm)
            for i in range(n)]


def test_blend_overrides_rain_and_sets_month():
    weather = _weather(3)
    rain = {weather[1]["date"]: 15.0}          # SIR reports rain on the middle day
    out = fwh.blend_records(weather, rain)
    assert out[0]["rain_mm"] == 0.0            # ERA5 value kept
    assert out[1]["rain_mm"] == 15.0           # SIR override
    assert all(r["month"] == r["date"].month for r in out)


def test_blend_defaults_missing_rain_to_zero():
    weather = [dict(date=dt.date(2026, 6, 10), temperature=30.0,
                    relative_humidity=20.0, wind_ms=3.0)]   # no rain_mm
    out = fwh.blend_records(weather)
    assert out[0]["rain_mm"] == 0.0


def test_spinup_state_dry_spell_high_drought():
    st = fwh.spinup_state(_weather(28, rain_mm=0.0))
    assert st.dc > 200.0                        # weeks of drought -> high DC
    assert st.dmc > 20.0


def test_observed_rain_lowers_drought_vs_era5_zero():
    weather = _weather(28, rain_mm=0.0)         # ERA5 says bone dry
    # SIR observed real rain on ~a third of the days
    rain = {}
    for i, w in enumerate(weather):
        if i % 3 == 0:
            rain[w["date"]] = 12.0
    dry = fwh.spinup_state(weather)             # ERA5-only
    wet = fwh.spinup_state(weather, rain)       # blended with observed rain
    assert wet.dc < dry.dc                       # observed rain relaxes the drought
    assert wet.dmc < dry.dmc


def test_scale_rain_to_observed_matches_total():
    weather = _weather(30, rain_mm=1.0)          # ERA5: 1 mm/day = 30 mm over window
    out = fwh.scale_rain_to_observed(weather, 9.9, over_days=30)
    assert sum(r["rain_mm"] for r in out) == pytest.approx(9.9)   # scaled to gauge
    # relative day-to-day structure preserved (all equal here) and inputs untouched
    assert all(r["rain_mm"] == pytest.approx(9.9 / 30) for r in out)
    assert sum(r["rain_mm"] for r in weather) == pytest.approx(30.0)


def test_scale_rain_spreads_when_era5_dry():
    weather = _weather(10, rain_mm=0.0)          # ERA5 bone dry but gauge saw rain
    out = fwh.scale_rain_to_observed(weather, 20.0, over_days=10)
    assert sum(r["rain_mm"] for r in out) == pytest.approx(20.0)
    assert all(r["rain_mm"] == pytest.approx(2.0) for r in out)


def test_scale_rain_none_is_noop():
    weather = _weather(5, rain_mm=3.0)
    out = fwh.scale_rain_to_observed(weather, None)
    assert [r["rain_mm"] for r in out] == [3.0] * 5


def test_observed_scaling_lowers_drought_vs_dry_era5():
    weather = _weather(30, rain_mm=0.0)          # ERA5 dry -> severe drought
    dry = fwh.spinup_state(weather)
    scaled = fwh.spinup_state(fwh.scale_rain_to_observed(weather, 60.0, over_days=30))
    # observed rain relaxes the duff moisture code (DMC, 1.5 mm/day threshold); the
    # deep Drought Code needs >2.8 mm/day so light spread rain leaves it unchanged --
    # correct FWI behaviour.
    assert scaled.dmc < dry.dmc
    assert scaled.dc <= dry.dc


def test_spinup_feeds_gridded_fire_day():
    prev = fwh.spinup_state(_weather(28, rain_mm=0.0))
    import numpy as np
    out = fwi.gridded_fwi(temperature=np.full((2, 2), 33.0),
                          relative_humidity=np.full((2, 2), 20.0),
                          wind_kmh=np.full((2, 2), 18.0), month=7, previous=prev)
    assert out["fwi"].shape == (2, 2)
    assert np.all(out["bui"] > 20.0)            # the spun-up drought raises BUI
