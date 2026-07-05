# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (c) 2026 Cristiano Foderi <cristiano.foderi@gmail.com>

"""Multi-week fire-weather history -> a spun-up FWI state for the fire day.

The FWI System's slow codes (DMC ~12-day, DC ~52-day lag) only mean something after
weeks of drive, so a fire-day FWI needs a spin-up over the preceding ~3-4 weeks.
This module builds that daily-noon history and runs it through
:func:`pyflam.fwi.spinup_state`.

Following the chosen data strategy, the drivers are **blended**: observed daily
rainfall from the SIR Toscana rain-gauge network (:mod:`pyflam.sir_toscana`) --
which governs DMC/DC and is captured far better by gauges than by reanalysis over
complex terrain -- combined with ERA5 reanalysis for the spatially-complete noon
temperature / relative humidity / wind (:func:`era5_daily_noon_history`). If SIR is
unavailable the ERA5 precipitation is used instead, so the spin-up always runs.

The blend (:func:`blend_records`) is pure and unit-tested; the ERA5 fetch and the
top-level orchestration (:func:`fire_weather_spinup`) touch the network / CDS and
are not exercised offline.
"""

from __future__ import annotations

import datetime as _dt
import warnings

import numpy as np

from . import fwi


def _as_date(d):
    return d.date() if isinstance(d, _dt.datetime) else d


def blend_records(weather, rain_by_date=None):
    """Merge daily-noon weather records with observed rainfall -> FWI records.

    ``weather`` is an iterable of dicts, each with a ``date`` and
    ``temperature`` (C) / ``relative_humidity`` (%) / ``wind_ms`` (and optionally an
    ERA5 ``rain_mm``). ``rain_by_date`` is a ``{date: mm}`` mapping (e.g. from
    :func:`pyflam.sir_toscana.daily_rainfall`); where a day is present it overrides
    that record's rainfall, otherwise the record's own ``rain_mm`` (or 0) is kept.
    Adds the ``month`` each code needs. Returns the list
    :func:`pyflam.fwi.fwi_from_records` / :func:`~pyflam.fwi.spinup_state` consume.
    """
    rain_by_date = rain_by_date or {}
    out = []
    for w in weather:
        r = dict(w)
        day = _as_date(r["date"])
        if day in rain_by_date and rain_by_date[day] is not None:
            r["rain_mm"] = float(rain_by_date[day])
        r.setdefault("rain_mm", 0.0)
        r["month"] = day.month
        out.append(r)
    return out


def spinup_state(weather, rain_by_date=None, *, state=None, hemisphere="north"):
    """Blend weather + observed rain and return the fire-day :class:`~pyflam.fwi.FWIState`."""
    records = blend_records(weather, rain_by_date)
    return fwi.spinup_state(records, state=state, hemisphere=hemisphere)


# --- ERA5 daily-noon history (network / CDS; not run offline) -----------------

ERA5_HISTORY_VARIABLES = [
    "2m_temperature", "2m_dewpoint_temperature",
    "10m_u_component_of_wind", "10m_v_component_of_wind",
    "total_precipitation",
]


def era5_daily_noon_history(latitude, longitude, start, end, *,
                            cache_path=None, hour_utc=12, box_deg=0.25,
                            force=False):  # pragma: no cover - network/CDS
    """Daily-noon ERA5 weather + 24 h precipitation over ``[start, end]`` at a point.

    Downloads hourly ERA5 (T, dewpoint, 10 m wind, total precipitation) for a small
    box around the point and reduces it to one record per day: noon temperature /
    relative humidity / wind speed, and the rainfall accumulated in the 24 h ending
    at that noon (the FWI convention). Needs ``cdsapi`` + CDS credentials and
    ``xarray``. Returns records for :func:`blend_records`.
    """
    import os

    from .atmosphere import (_read_atmosphere_dataset, kelvin_to_celsius,
                             relative_humidity_from_dewpoint)

    start, end = _as_date(start), _as_date(end)
    cache_path = cache_path or f"era5_hist_{start}_{end}.nc"
    if force or not os.path.exists(cache_path):
        import cdsapi
        # one hourly request spanning the window (+1 day for the noon-to-noon tail)
        days = [(start + _dt.timedelta(days=i))
                for i in range((end - start).days + 2)]
        cdsapi.Client().retrieve("reanalysis-era5-single-levels", {
            "product_type": "reanalysis", "format": "netcdf",
            "variable": ERA5_HISTORY_VARIABLES,
            "date": f"{days[0]}/{days[-1]}",
            "time": [f"{h:02d}:00" for h in range(24)],
            "area": [latitude + box_deg, longitude - box_deg,
                     latitude - box_deg, longitude + box_deg],
        }, cache_path)

    ds = _read_atmosphere_dataset(cache_path)       # handles the CDS zip (instant+accum)
    tname = "valid_time" if "valid_time" in ds.coords else "time"
    pt = ds.sel(latitude=latitude, longitude=longitude, method="nearest")
    t = np.asarray(pt[tname].values).astype("datetime64[h]")

    def _series(name):
        for cand in (name, name.upper()):
            if cand in pt:
                return np.asarray(pt[cand].values, dtype=float)
        raise KeyError(name)

    t2m = kelvin_to_celsius(_series("t2m"))
    d2m = kelvin_to_celsius(_series("d2m"))
    u10, v10 = _series("u10"), _series("v10")
    tp_mm = _series("tp") * 1000.0                      # m -> mm (hourly accum)
    wind = np.hypot(u10, v10)
    rh = np.asarray(relative_humidity_from_dewpoint(t2m, d2m), dtype=float)

    hours = t.astype("datetime64[h]")
    records = []
    day = start
    while day <= end:
        noon = np.datetime64(f"{day}T{hour_utc:02d}", "h")
        i = int(np.argmin(np.abs(hours - noon)))
        # rain accumulated in the 24 h ending at noon (FWI's noon-to-noon window)
        lo = noon - np.timedelta64(24, "h")
        window = (hours > lo) & (hours <= noon)
        records.append(dict(
            date=day, temperature=float(t2m[i]), relative_humidity=float(rh[i]),
            wind_ms=float(wind[i]), rain_mm=float(np.nansum(tp_mm[window]))))
        day += _dt.timedelta(days=1)
    return records


def fire_weather_spinup(ignition_date, latitude, longitude, *, days=28,
                        use_sir=True, hemisphere="north", cache_path=None,
                        log=print):  # pragma: no cover - network
    """Fire-day :class:`~pyflam.fwi.FWIState` from ~``days`` of blended history.

    Fetches ERA5 daily-noon weather for the ``days`` before ``ignition_date`` and,
    when ``use_sir``, overrides the daily rainfall with the nearest SIR Toscana
    gauge (falling back to ERA5 precipitation if SIR is unreachable). Runs the FWI
    spin-up and returns the state to carry into the fire-day (gridded) FWI. Degrades
    to the spring-startup :class:`~pyflam.fwi.FWIState` if ERA5 is unavailable.
    """
    ignition_date = _as_date(ignition_date)
    start = ignition_date - _dt.timedelta(days=days)
    end = ignition_date - _dt.timedelta(days=1)
    try:
        weather = era5_daily_noon_history(latitude, longitude, start, end,
                                          cache_path=cache_path)
    except Exception as exc:
        log(f"ERA5 history unavailable ({exc}); using spring-startup FWI state")
        return fwi.FWIState()

    rain_by_date = None
    if use_sir:
        try:
            from . import sir_toscana as sir
            rain_by_date = sir.daily_rainfall(latitude, longitude, start, end)
            log(f"SIR Toscana rain: {len(rain_by_date)} gauge-days "
                f"({sum(rain_by_date.values()):.0f} mm total over the spin-up)")
        except Exception as exc:
            log(f"SIR Toscana unavailable ({exc}); using ERA5 precipitation")

    st = spinup_state(weather, rain_by_date, hemisphere=hemisphere)
    log(f"FWI spin-up ({days} d -> {end}): FFMC {st.ffmc:.1f}  DMC {st.dmc:.1f}  "
        f"DC {st.dc:.1f}")
    return st
