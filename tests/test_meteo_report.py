"""Tests for the near-real-time meteo variation report."""

from __future__ import annotations

from datetime import datetime

import numpy as np
import pytest

import pyflam
from pyflam import atmosphere as atm


class _Diurnal(atm.AtmosphereProvider):
    """Synthetic afternoon: warming, drying, wind rising, atmosphere destabilising."""

    def state_at(self, lat, lon, time):
        h = time.hour + time.minute / 60.0
        return atm.AtmosphericState(
            wind_speed=2.0 + 0.7 * (h - 10), wind_direction=290.0 - 4.0 * (h - 12),
            temperature=18.0 + 1.6 * (h - 10),
            relative_humidity=max(12.0, 65.0 - 5.0 * (h - 10)),
            cape=150.0 * (h - 10), sensible_heat_flux=35.0 * (h - 10))


def _report():
    return pyflam.build_meteo_report(
        _Diurnal(), location=(43.0, 11.0),
        start_time=datetime(2026, 8, 1, 11, 0), total_minutes=360,
        step_minutes=60)


def test_report_length_and_records():
    rep = _report()
    assert len(rep.times) == 7                       # 11:00..17:00 inclusive
    recs = rep.to_records()
    assert len(recs) == 7
    assert {"time", "temperature", "m_1h", "stability"} <= set(recs[0])


def test_variation_trends():
    var = _report().variation()
    assert var["relative_humidity"]["trend"] == "falling"
    assert var["temperature"]["trend"] == "rising"
    assert var["wind_speed"]["trend"] == "rising"
    assert var["cape"]["trend"] == "rising"
    assert var["relative_humidity"]["change"] < 0
    assert var["relative_humidity"]["range"] > 0


def test_fuel_moisture_evolves_with_lag():
    rep = _report()
    m1 = [v for v in rep.series["m_1h"]]
    m100 = [v for v in rep.series["m_100h"]]
    # Drying afternoon: 1-h moisture changes more than 100-h (memory).
    assert abs(m1[-1] - m1[0]) > abs(m100[-1] - m100[0])


def test_stability_tracked():
    rep = _report()
    assert len(rep.stability) == len(rep.times)
    assert all(s in ("unstable", "neutral", "stable") for s in rep.stability)


def test_summary_is_text():
    s = _report().summary()
    assert "relative_humidity" in s and "stability" in s
    assert isinstance(s, str)


def test_end_time_equivalent_to_total_minutes():
    a = _report()
    b = pyflam.build_meteo_report(
        _Diurnal(), location=(43.0, 11.0),
        start_time=datetime(2026, 8, 1, 11, 0),
        end_time=datetime(2026, 8, 1, 17, 0), step_minutes=60)
    assert len(a.times) == len(b.times)


def test_requires_window():
    with pytest.raises(ValueError):
        pyflam.build_meteo_report(_Diurnal(), location=(43.0, 11.0),
                                  start_time=datetime(2026, 8, 1, 11, 0))
