"""Tests for the near-real-time run product (weather -> spread -> reports)."""

from __future__ import annotations

from datetime import datetime

import numpy as np
import pytest

import pyflam
from pyflam import atmosphere as atm


class _Diurnal(atm.AtmosphereProvider):
    def state_at(self, lat, lon, time):
        h = time.hour + time.minute / 60.0
        return atm.AtmosphericState(
            wind_speed=3.0 + 0.7 * (h - 10), wind_direction=270.0 - 3.0 * (h - 10),
            temperature=20.0 + 1.5 * (h - 10),
            relative_humidity=max(12.0, 55.0 - 5.0 * (h - 10)),
            cape=150.0 * (h - 10), sensible_heat_flux=40.0 * (h - 10))


def _landscape(n=81):
    fuel = np.full((n, n), 104, dtype=int)
    fuel[:, 50:] = 147
    slope = np.zeros((n, n))
    for r in range(n):
        slope[r, :] = max(0.0, 30 - r)
    return pyflam.Landscape(
        fuel_model=fuel, slope=slope, aspect=np.full((n, n), 180.0),
        cellsize_x=100.0, cellsize_y=100.0, west=4_300_000.0, north=2_370_000.0,
        slope_units="degrees", crs="EPSG:3035")


def _run():
    return pyflam.run_realtime(
        _landscape(), [(40, 40)], atmosphere=_Diurnal(), location=(43.0, 11.0),
        start_time=datetime(2026, 8, 1, 11, 0), total_time=120, dt=30,
        m_live_herb=0.6, m_live_woody=0.9)


def test_run_product_components():
    prod = _run()
    assert isinstance(prod, pyflam.RunProduct)
    assert isinstance(prod.meteo, pyflam.MeteoReport)
    assert isinstance(prod.operative, pyflam.OperativeReport)
    assert np.isfinite(prod.arrival_time).sum() > 1


def test_run_product_summary_and_geojson():
    prod = _run()
    s = prod.summary()
    assert "Near-real-time run" in s and "Meteo variation" in s \
        and "Operative analysis" in s
    gj = prod.operative_geojson(_landscape())
    assert gj["type"] == "FeatureCollection"
    assert any(f["properties"]["kind"] == "perimeter" for f in gj["features"])


def test_run_product_meteo_window_matches():
    prod = _run()
    # 120 min at 60-min steps -> 11:00, 12:00, 13:00
    assert len(prod.meteo.times) == 3


def test_run_product_write_geojson(tmp_path):
    import json
    prod = _run()
    path = tmp_path / "nrt.geojson"
    prod.write_geojson(_landscape(), str(path))
    gj = json.loads(path.read_text())
    assert gj["features"]


def test_run_no_fire_has_no_operative():
    """A nonburnable landscape yields no perimeter; product still returns."""
    n = 21
    ls = pyflam.Landscape(
        fuel_model=np.full((n, n), 91, dtype=int),   # NB1 nonburnable
        slope=np.zeros((n, n)), aspect=np.zeros((n, n)),
        cellsize_x=100.0, cellsize_y=100.0, west=0.0, north=n * 100.0,
        slope_units="degrees")
    prod = pyflam.run_realtime(
        ls, [(10, 10)], atmosphere=_Diurnal(), location=(43.0, 11.0),
        start_time=datetime(2026, 8, 1, 11, 0), total_time=60, dt=30,
        m_live_herb=0.6, m_live_woody=0.9)
    assert prod.operative is None
    assert "no perimeter" in prod.summary()
