# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (c) 2026 Cristiano Foderi <cristiano.foderi@gmail.com>

"""Client for SIR Toscana (Servizio Idrologico Regionale) rain gauges.

Pulls **observed cumulative rainfall** from the Tuscany regional hydrological
network -- the precipitation that drives the FWI System's DMC/DC, which ERA5
captures poorly over complex Tuscan terrain.

    from pyflam import sir_toscana as sir
    obs = sir.recent_cumulative(43.936, 11.096)     # nearest gauge to the point
    obs.cumulative[30], obs.dry_days                # e.g. (9.9 mm, 6 days)

**Confirmed endpoint (2026-07).** The near-real-time monitoring page
``/monitoraggio/stazioni.php?type=pluvio_men`` returns HTTP 200 with an inline
JavaScript ``VALUES`` array -- one row per gauge with the rainfall accumulated over
the **1/2/5/7/10/15/30 days ending "now"**, plus a dry-day count. That is what
:func:`fetch_recent_cumulative` / :func:`recent_cumulative` parse (schema in
:func:`parse_pluvio_men`). It is *near-real-time* (anchored to today), not an
arbitrary historical daily series.

**Still open (TODO).** An arbitrary past **daily** series (needed to spin the FWI
up before a fire that is not recent) lives only in the session-based archive
(``ricerca-dati``), which has no clean API -- so for a historical run export the
series from the SIR site and load it with :func:`read_sir_csv` (the guaranteed
path), or use ERA5 precipitation. The pure parsers and station selection are
unit-tested; the HTTP calls are not.

Data © Servizio Idrologico Regionale della Toscana; see
https://www.sir.toscana.it and https://dati.toscana.it/dataset/pluviometri.
"""

from __future__ import annotations

import csv as _csv
import datetime as _dt
import math
import re
from dataclasses import dataclass, field

# Public monitoring host and the confirmed near-real-time endpoint.
SIR_MONITOR_BASE = "https://www.sir.toscana.it"
PLUVIO_MEN_PATH = "/monitoraggio/stazioni.php?type=pluvio_men"
# Cumulative windows (days) the pluvio_men table reports, in VALUES column order.
CUMULATIVE_WINDOWS = (1, 2, 5, 7, 10, 15, 30)
# Real gauges near the Calvana / Prato ignition, with the confirmed SIR codes (from
# the live pluvio_men table) and approximate coordinates -- used to pick the nearest
# station by lat/lon (the pluvio_men table itself carries no coordinates).
FALLBACK_STATIONS = (
    ("Vaiano acquedotto", "TOS11000503", 43.962, 11.121),
    ("Prato Università", "TOS01001205", 43.880, 11.098),
    ("Fattoria Iavello", "TOS01001273", 43.948, 11.033),
    ("Cantagallo", "TOS01001151", 44.013, 11.065),
    ("Gamberame", "TOS01004779", 44.001, 11.118),
    ("Vernio", "TOS01001171", 44.043, 11.149),
)


class SIRError(RuntimeError):
    """Raised when SIR data cannot be fetched or parsed as expected."""


@dataclass(frozen=True)
class SIRStation:
    name: str
    code: str
    latitude: float
    longitude: float


@dataclass(frozen=True)
class SIRObservation:
    """One gauge's near-real-time rainfall from the pluvio_men table.

    ``cumulative`` maps each window (days in :data:`CUMULATIVE_WINDOWS`) to the mm
    accumulated over that window ending ``as_of``; ``dry_days`` is the reported
    consecutive dry-day count; ``today_mm`` the accumulation since local midnight.
    """
    code: str
    name: str
    comune: str
    province: str
    elevation_m: float | None
    today_mm: float | None
    as_of: str
    dry_days: int | None
    cumulative: dict = field(default_factory=dict)


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

# One VALUES row: ("code","name","comune","prov","zone","today","asof",
#                  1g,2g,5g,7g,10g,15g,30g,"dry","elev","flag")  -- 17 fields.
_VALUES_RE = re.compile(r"VALUES\[\d+\]\s*=\s*new\s+Array\((.*?)\);", re.S)
_FIELD_RE = re.compile(r'"((?:[^"\\]|\\.)*)"')


def _num(token):
    if token is None:
        return None
    t = re.sub(r"<[^>]+>", "", str(token)).replace("\xa0", " ").strip()
    t = t.replace(",", ".")
    m = re.search(r"[-+]?\d+(?:\.\d+)?", t)
    return float(m.group(0)) if m else None


def parse_pluvio_men(html):
    """Parse the SIR ``pluvio_men`` page -> ``{code: SIRObservation}``.

    Extracts the inline ``VALUES[i] = new Array( … )`` rows (the near-real-time
    rain-gauge table) into per-gauge cumulative rainfall over
    :data:`CUMULATIVE_WINDOWS` plus the dry-day count. Raises :class:`SIRError` if
    no rows are found. Pure -- the confirmed schema is unit-tested against a fixture.
    """
    out = {}
    for row in _VALUES_RE.findall(html):
        f = [m.group(1) for m in _FIELD_RE.finditer(row)]
        if len(f) < 16:
            continue
        code = f[0].strip()
        cum = {w: _num(f[7 + i]) for i, w in enumerate(CUMULATIVE_WINDOWS)}
        dry = _num(f[14])
        out[code] = SIRObservation(
            code=code, name=f[1].strip(), comune=f[2].strip(), province=f[3].strip(),
            elevation_m=_num(f[15]) if len(f) > 15 else None,
            today_mm=_num(f[5]), as_of=f[6].strip(),
            dry_days=int(dry) if dry is not None else None, cumulative=cum)
    if not out:
        raise SIRError("no VALUES rows found in pluvio_men page")
    return out

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


def station_catalog():
    """The built-in Calvana/Prato gauge catalogue (real SIR codes + approx coords).

    The pluvio_men table carries no coordinates, so this small list is what
    :func:`recent_cumulative` uses to pick the nearest gauge by lat/lon. Extend it
    (or match on ``comune``) for other areas; the geo catalogue with full
    coordinates is at https://dati.toscana.it/dataset/pluviometri.
    """
    return [SIRStation(n, c, la, lo) for (n, c, la, lo) in FALLBACK_STATIONS]


# --- confirmed HTTP: near-real-time pluvio_men (network; not unit-tested) ------

def fetch_recent_cumulative(*, base=SIR_MONITOR_BASE, timeout=30):
    """Fetch + parse the SIR ``pluvio_men`` table -> ``{code: SIRObservation}``.

    Hits the **confirmed** endpoint ``/monitoraggio/stazioni.php?type=pluvio_men``
    (HTTP 200, verified 2026-07) and parses the near-real-time cumulative rainfall
    for every gauge (:func:`parse_pluvio_men`). Raises :class:`SIRError` on any
    network/parse problem so callers can fall back to ERA5.
    """
    try:                                            # pragma: no cover - network
        import requests
    except ImportError as exc:                      # pragma: no cover
        raise SIRError("SIR fetch needs `requests`") from exc
    try:                                            # pragma: no cover - network
        r = requests.get(f"{base}{PLUVIO_MEN_PATH}", timeout=timeout,
                         headers={"User-Agent": "pyflam/0.2 (fire-weather)"})
        r.raise_for_status()
        return parse_pluvio_men(r.text)
    except SIRError:                                # pragma: no cover - network
        raise
    except Exception as exc:                        # pragma: no cover - network
        raise SIRError(f"SIR pluvio_men request/parse failed: {exc}") from exc


def recent_cumulative(latitude, longitude, *, station=None, catalog=None,
                      base=SIR_MONITOR_BASE, timeout=30):
    """Nearest gauge's near-real-time cumulative rainfall as a :class:`SIRObservation`.

    Picks the nearest catalogue gauge to the point (or the given ``station`` code),
    fetches :func:`fetch_recent_cumulative`, and returns that gauge's observation
    (windowed cumulatives ending "now" + dry-day count). Raises :class:`SIRError`
    if the gauge is absent from the live table. **Near-real-time (anchored to
    today)** -- for an arbitrary past daily series use :func:`read_sir_csv`.
    """
    obs = fetch_recent_cumulative(base=base, timeout=timeout)   # pragma: no cover
    code = station                                              # pragma: no cover
    if code is None:                                            # pragma: no cover
        code = nearest_station(catalog or station_catalog(),
                               latitude, longitude).code
    if code not in obs:                                         # pragma: no cover
        # fall back to the geographically nearest gauge actually present
        present = [s for s in (catalog or station_catalog()) if s.code in obs]
        if not present:
            raise SIRError(f"gauge {code} (and fallbacks) absent from pluvio_men")
        code = nearest_station(present, latitude, longitude).code
    return obs[code]                                            # pragma: no cover
