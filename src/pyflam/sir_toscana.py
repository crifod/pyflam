# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (c) 2026 Cristiano Foderi <cristiano.foderi@gmail.com>

"""Best-effort client for SIR Toscana (Servizio Idrologico Regionale) rain gauges.

Pulls **daily cumulative rainfall** (mm/24 h) from the Tuscany regional hydrological
network for a point + date window -- the observed precipitation that drives the
FWI System's DMC/DC (which ERA5 captures poorly over complex Tuscan terrain).

    from pyflam import sir_toscana as sir
    rain = sir.daily_rainfall(43.936, 11.096, date(2026, 6, 5), date(2026, 7, 3))
    # -> {date(2026, 6, 5): 0.0, ..., date(2026, 7, 3): 0.0}

**Endpoint caveat.** SIR does not publish a documented time-series API; the
historical archive (``ricerca-dati``) is a session-based form and the monitoring
pages (``stazioni.php``) are HTML tables. This module targets those observed
structures on a best-effort basis -- the exact endpoint/markup may change, so
:func:`daily_rainfall` raises :class:`SIRError` on any deviation rather than
returning wrong data, and callers should fall back (e.g. to ERA5 precipitation).
The parsing (:func:`parse_rain_table`) and station selection
(:func:`nearest_station`) are pure and unit-tested; for a guaranteed path, export
the series from the SIR site and load it with :func:`read_sir_csv`.

Data © Servizio Idrologico Regionale della Toscana; see
https://www.sir.toscana.it and https://dati.toscana.it/dataset/pluviometri.
"""

from __future__ import annotations

import csv as _csv
import datetime as _dt
import math
import re
from dataclasses import dataclass

# Public monitoring host (mirrors www.sir.toscana.it). Overridable for testing.
SIR_MONITOR_BASE = "https://www.sir.toscana.it"
# A small built-in catalogue of gauges near the Calvana / Prato area, used when the
# online pluviometer catalogue cannot be fetched. (name, code, lat, lon) -- approx.
FALLBACK_STATIONS = (
    ("Prato (Galceti)", "TOS11000108", 43.9075, 11.0906),
    ("Vaiano", "TOS11000110", 43.9600, 11.1200),
    ("Vernio", "TOS11000112", 44.0430, 11.1490),
    ("Cantagallo", "TOS11000114", 44.0130, 11.0700),
    ("Calenzano", "TOS11000116", 43.8580, 11.1630),
)


class SIRError(RuntimeError):
    """Raised when SIR data cannot be fetched or parsed as expected."""


@dataclass(frozen=True)
class SIRStation:
    name: str
    code: str
    latitude: float
    longitude: float


def _haversine_km(lat1, lon1, lat2, lon2):
    r = 6371.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp = math.radians(lat2 - lat1)
    dl = math.radians(lon2 - lon1)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * r * math.asin(math.sqrt(a))


def nearest_station(stations, latitude, longitude):
    """The :class:`SIRStation` closest (great-circle) to a point."""
    stations = list(stations)
    if not stations:
        raise SIRError("no stations to choose from")
    return min(stations, key=lambda s: _haversine_km(
        latitude, longitude, s.latitude, s.longitude))


# --- pure parsers (unit-tested) ----------------------------------------------

_NUM = re.compile(r"[-+]?\d+(?:[.,]\d+)?")
_DATE_FORMATS = ("%d/%m/%Y", "%Y-%m-%d", "%d-%m-%Y", "%d/%m/%y")


def _parse_date(token):
    token = token.strip()
    for fmt in _DATE_FORMATS:
        try:
            return _dt.datetime.strptime(token, fmt).date()
        except ValueError:
            continue
    return None


def _parse_float(token):
    m = _NUM.search(token.replace("\xa0", " "))
    if not m:
        return None
    val = m.group(0).replace(",", ".")
    try:
        return float(val)
    except ValueError:
        return None


def parse_rain_table(text):
    """Parse SIR daily-rainfall text -> ``{date: mm}`` (robust to CSV/TSV/HTML rows).

    Accepts a two-column ``date, rainfall`` table in any of the SIR/CFR export
    shapes: ``;``/``,``/tab/whitespace-separated CSV, or HTML ``<tr><td>`` rows.
    Italian ``dd/mm/yyyy`` dates and ``,`` decimals are handled; non-data lines
    (headers, totals) are skipped. Raises :class:`SIRError` if nothing parses.
    """
    # flatten simple HTML table rows to "cell|cell" lines
    if "<t" in text.lower():
        rows = re.findall(r"<tr[^>]*>(.*?)</tr>", text, flags=re.I | re.S)
        lines = ["|".join(re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", r, flags=re.I | re.S))
                 for r in rows]
        lines = [re.sub(r"<[^>]+>", " ", ln) for ln in lines]
    else:
        lines = text.splitlines()

    out = {}
    for line in lines:
        line = line.strip()
        if not line:
            continue
        # split on the common separators, keeping tokens
        parts = re.split(r"[;\t|]|\s{2,}|,(?=\s)", line)
        parts = [p for p in (p.strip() for p in parts) if p]
        if len(parts) < 2:
            parts = line.split()
        day = next((d for d in (_parse_date(p) for p in parts) if d is not None), None)
        if day is None:
            continue
        # the rainfall is the first numeric token that is not itself the date
        rain = None
        for p in parts:
            if _parse_date(p) is not None:
                continue
            rain = _parse_float(p)
            if rain is not None:
                break
        if rain is not None:
            out[day] = rain
    if not out:
        raise SIRError("no date/rainfall rows found in SIR response")
    return out


def read_sir_csv(path, *, date_column=0, rain_column=1, delimiter=None):
    """Load a SIR rainfall CSV export -> ``{date: mm}`` (the guaranteed path).

    Point ``date_column`` / ``rain_column`` at the right columns (defaults: date in
    the first, rainfall in the second). ``delimiter`` autodetects ``;``/``,``/tab
    when ``None``. Italian dates and ``,`` decimals are handled.
    """
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        if delimiter is None:
            try:
                delimiter = _csv.Sniffer().sniff(sample, ";,\t").delimiter
            except _csv.Error:
                delimiter = ";"
        out = {}
        for row in _csv.reader(fh, delimiter=delimiter):
            if len(row) <= max(date_column, rain_column):
                continue
            day = _parse_date(row[date_column])
            rain = _parse_float(row[rain_column])
            if day is not None and rain is not None:
                out[day] = rain
    if not out:
        raise SIRError(f"no date/rainfall rows parsed from {path}")
    return out


# --- best-effort HTTP (network; not exercised in unit tests) ------------------

def station_catalog(*, timeout=30):
    """Best-effort pluviometer catalogue; falls back to a built-in local list."""
    try:                                            # pragma: no cover - network
        import requests
        url = f"{SIR_MONITOR_BASE}/archivio/dati.php?D=json_stations"
        r = requests.get(url, timeout=timeout)
        r.raise_for_status()
        data = r.json()
        stations = [SIRStation(name=s.get("nome") or s.get("name") or s["codice"],
                               code=str(s.get("codice") or s.get("code")),
                               latitude=float(s["lat"]), longitude=float(s["lon"]))
                    for s in (data.get("stazioni") or data.get("stations") or data)]
        if stations:
            return stations
    except Exception:                               # pragma: no cover - network
        pass
    return [SIRStation(n, c, la, lo) for (n, c, la, lo) in FALLBACK_STATIONS]


def daily_rainfall(latitude, longitude, start, end, *, station=None,
                   base=SIR_MONITOR_BASE, timeout=30):
    """Daily cumulative rainfall ``{date: mm}`` at the nearest gauge (best-effort).

    Selects the nearest catalogue station (or the given ``station``) and requests
    its daily series for ``[start, end]`` from the SIR monitoring endpoint, parsing
    with :func:`parse_rain_table`. Raises :class:`SIRError` on any network/parse
    problem so the caller can fall back to ERA5 precipitation. ``start``/``end`` are
    :class:`datetime.date`. **Unverified endpoint -- see the module docstring.**
    """
    try:                                            # pragma: no cover - network
        import requests
    except ImportError as exc:                      # pragma: no cover
        raise SIRError("SIR fetch needs `requests`") from exc
    if station is None:                             # pragma: no cover - network
        station = nearest_station(station_catalog(timeout=timeout),
                                  latitude, longitude)
    try:                                            # pragma: no cover - network
        # TODO(confirm): this endpoint + parameter names are a best-effort guess and
        # have so far only returned HTTP 429 (reachable but rate-limited), never a
        # verified 200. Confirm `stazioni.php?type=pluvio_day` (cod/from/to) against a
        # live response and adjust parse_rain_table if the markup differs; until then
        # read_sir_csv is the guaranteed path.
        params = {"type": "pluvio_day", "cod": station.code,
                  "from": start.strftime("%d/%m/%Y"), "to": end.strftime("%d/%m/%Y")}
        r = requests.get(f"{base}/monitoraggio/stazioni.php", params=params,
                         timeout=timeout)
        r.raise_for_status()
        series = parse_rain_table(r.text)
    except SIRError:                                # pragma: no cover - network
        raise
    except Exception as exc:                        # pragma: no cover - network
        raise SIRError(f"SIR request/parse failed: {exc}") from exc
    return {d: v for d, v in series.items() if start <= d <= end}
