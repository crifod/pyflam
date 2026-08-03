"""Tests for the SIR Toscana rain-gauge client's pure parsing / selection logic.

The network paths (daily_rainfall / station_catalog HTTP) are not exercised here --
they hit an unverified live endpoint. What is tested is the robust parsing of the
export shapes and the nearest-gauge selection.
"""

from __future__ import annotations

import datetime as dt

import pytest

from pyflam import sir_toscana as sir


def test_parse_semicolon_italian_decimals():
    text = ("data;pioggia_mm\n"
            "05/06/2026;0,0\n"
            "06/06/2026;12,6\n"
            "07/06/2026;0,2\n"
            "Totale;12,8\n")
    out = sir.parse_rain_table(text)
    assert out == {dt.date(2026, 6, 5): 0.0, dt.date(2026, 6, 6): 12.6,
                   dt.date(2026, 6, 7): 0.2}


def test_parse_html_table():
    text = ("<table><tr><th>Data</th><th>mm</th></tr>"
            "<tr><td>01/07/2026</td><td>0.0</td></tr>"
            "<tr><td>02/07/2026</td><td>3.5</td></tr></table>")
    out = sir.parse_rain_table(text)
    assert out[dt.date(2026, 7, 1)] == 0.0
    assert out[dt.date(2026, 7, 2)] == 3.5


def test_parse_iso_and_whitespace():
    text = "2026-06-30 1.2\n2026-07-01 0.0\n"
    out = sir.parse_rain_table(text)
    assert out[dt.date(2026, 6, 30)] == 1.2
    assert out[dt.date(2026, 7, 1)] == 0.0


def test_parse_raises_on_garbage():
    with pytest.raises(sir.SIRError):
        sir.parse_rain_table("no data here\nstill nothing\n")


def test_read_sir_csv(tmp_path):
    p = tmp_path / "rain.csv"
    p.write_text("data;valore\n05/06/2026;0,0\n06/06/2026;4,8\n", encoding="utf-8")
    out = sir.read_sir_csv(p)
    assert out == {dt.date(2026, 6, 5): 0.0, dt.date(2026, 6, 6): 4.8}


def test_nearest_station_picks_closest():
    stations = [
        sir.SIRStation("far", "A", 40.0, 8.0),
        sir.SIRStation("near", "B", 43.94, 11.10),
        sir.SIRStation("mid", "C", 43.0, 11.0),
    ]
    got = sir.nearest_station(stations, 43.936, 11.096)
    assert got.name == "near"


def test_nearest_station_empty_raises():
    with pytest.raises(sir.SIRError):
        sir.nearest_station([], 43.0, 11.0)


# --- confirmed pluvio_men schema (real rows from the live 200) -----------------

PLUVIO_MEN_FIXTURE = (
    "var VALUES = new Array();\n"
    'VALUES[0] = new Array("TOS11000503","Vaiano acquedotto","Vaiano","PO","B",'
    '"0.0","05/07 10.15","0","0","<b>2.4</b>","<b>9.9</b>","<b>9.9</b>",'
    '"<b>9.9</b>","<b>9.9</b>","<b>6</b>","322","1");\n'
    'VALUES[1] = new Array("TOS01001205","Prato Università","Prato","PO","B",'
    '"0.0","05/07 10.15","0","0","<b>2</b>","<b>29.3</b>","<b>29.3</b>",'
    '"<b>29.3</b>","<b>29.9</b>","<b>6</b>","65","1");\n'
)


def test_parse_pluvio_men_schema():
    obs = sir.parse_pluvio_men(PLUVIO_MEN_FIXTURE)
    assert set(obs) == {"TOS11000503", "TOS01001205"}
    v = obs["TOS11000503"]
    assert v.name == "Vaiano acquedotto" and v.comune == "Vaiano" and v.province == "PO"
    assert v.elevation_m == 322.0 and v.dry_days == 6 and v.today_mm == 0.0
    # cumulative windows in day order, <b> stripped, Italian-safe
    assert v.cumulative == {1: 0.0, 2: 0.0, 5: 2.4, 7: 9.9, 10: 9.9, 15: 9.9, 30: 9.9}
    assert obs["TOS01001205"].cumulative[30] == 29.9


def test_parse_pluvio_men_raises_when_empty():
    with pytest.raises(sir.SIRError):
        sir.parse_pluvio_men("<html>no VALUES array here</html>")


def test_station_catalog_has_real_codes_near_calvana():
    cat = sir.station_catalog()
    codes = {s.code for s in cat}
    assert "TOS11000503" in codes            # Vaiano acquedotto (real SIR code)
    near = sir.nearest_station(cat, 43.936, 11.096)
    assert near.code in codes and near.latitude and near.longitude
