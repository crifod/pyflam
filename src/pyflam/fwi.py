# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (c) 2026 Cristiano Foderi <cristiano.foderi@gmail.com>

"""Canadian Forest Fire Weather Index (FWI) System.

The FWI System (Van Wagner 1987; Van Wagner & Pickett 1985) turns **daily noon**
weather -- temperature, relative humidity, 10 m wind and the previous 24 h rain --
into six standard fire-danger codes, each a physically-based moisture bookkeeping
of a different fuel layer:

* **FFMC** Fine Fuel Moisture Code -- litter / cured fine fuels (fast, ~2/3-day
  time lag). Its moisture equivalent is directly comparable to the 1-h dead fuel
  moisture (:func:`ffmc_to_moisture`), so it is a natural cross-check for the
  weather-driven dead-fuel model.
* **DMC** Duff Moisture Code -- loosely compacted organic layers (~12-day lag).
* **DC** Drought Code -- deep compact organic layers (~52-day lag): the seasonal
  drought memory that a single day's weather cannot show.
* **ISI** Initial Spread Index -- wind + FFMC, an expected rate-of-spread proxy.
* **BUI** Buildup Index -- DMC + DC, the fuel available to a spreading fire.
* **FWI** the final index -- ISI + BUI, a general fire-intensity rating.

Inputs come from any pyflam gridded provider (ERA5 / GFS / ICON) via
``state_at`` -- see :func:`fwi_from_provider`; wind is converted from m/s
(:func:`wind_ms_to_kmh`). Precipitation is *not* carried on
:class:`~pyflam.atmosphere.AtmosphericState`, so the 24 h noon-to-noon rain is
passed explicitly (0 for a dry spell, or summed from ERA5 ``tp`` / GFS ``apcp``).

Because DMC and DC integrate weeks of weather, a meaningful run needs either a
multi-day drive or realistic startup codes; a single day off the spring-startup
defaults (FFMC 85, DMC 6, DC 15) is indicative for the fast FFMC/ISI only. All
functions are scalar- or NumPy-array-safe (a gridded FWI is just array inputs).

References:
    Van Wagner, C.E. 1987. Development and structure of the Canadian Forest Fire
        Weather Index System. Canadian Forestry Service, Forestry Technical
        Report 35. Ottawa.
    Van Wagner, C.E.; Pickett, T.L. 1985. Equations and FORTRAN program for the
        Canadian Forest Fire Weather Index System. Forestry Technical Report 33.
    Wang, Y.; Anderson, K.R.; Suddaby, R.M. 2015. Updated source code for
        calculating fire danger indices in the Canadian Forest Fire Weather Index
        System. Natural Resources Canada, Information Report NOR-X-424.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta

import numpy as np

# Effective day-length factors by month (index 0 = January), northern hemisphere.
# DMC uses an effective day length Le; DC a day-length adjustment Lf.
DMC_DAY_LENGTH = (6.5, 7.5, 9.0, 12.8, 13.9, 13.9, 12.4, 10.9, 9.4, 8.0, 7.0, 6.0)
DC_DAY_LENGTH = (-1.6, -1.6, -1.6, 0.9, 3.8, 5.8, 6.4, 5.0, 2.4, 0.4, -1.6, -1.6)

# Spring-startup codes after snowmelt (Van Wagner 1987); standard when no
# overwintering / prior run is available.
STARTUP_FFMC = 85.0
STARTUP_DMC = 6.0
STARTUP_DC = 15.0


def _scalar(x):
    """0-d array -> Python float; higher-rank arrays pass through."""
    a = np.asarray(x, dtype=float)
    return float(a) if a.ndim == 0 else a


def _month_factor(table, month):
    """Day-length factor for ``month`` (1-12) from a 12-tuple, scalar or array."""
    idx = np.asarray(month, dtype=int) - 1
    vals = np.asarray(table, dtype=float)[np.clip(idx, 0, 11)]
    return _scalar(vals)


def wind_ms_to_kmh(wind_ms):
    """10 m wind m/s -> km/h (the FWI System's wind unit)."""
    return _scalar(np.asarray(wind_ms, dtype=float) * 3.6)


def southern_hemisphere_tables():
    """DMC/DC day-length tables shifted six months for the southern hemisphere."""
    return (DMC_DAY_LENGTH[6:] + DMC_DAY_LENGTH[:6],
            DC_DAY_LENGTH[6:] + DC_DAY_LENGTH[:6])


# --- the six codes ------------------------------------------------------------

def ffmc(temp_c, relative_humidity, wind_kmh, rain_mm, ffmc_prev):
    """Fine Fuel Moisture Code from noon weather and yesterday's FFMC.

    Fast-responding litter moisture code (0-101, higher = drier). ``rain_mm`` is
    the previous 24 h rainfall; ``ffmc_prev`` yesterday's code (see
    :data:`STARTUP_FFMC`). Scalar or array.
    """
    T = np.asarray(temp_c, dtype=float)
    H = np.clip(np.asarray(relative_humidity, dtype=float), 0.0, 100.0)
    W = np.maximum(np.asarray(wind_kmh, dtype=float), 0.0)
    ro = np.maximum(np.asarray(rain_mm, dtype=float), 0.0)
    Fo = np.asarray(ffmc_prev, dtype=float)

    mo = 147.2 * (101.0 - Fo) / (59.5 + Fo)
    # rainfall wetting (only when > 0.5 mm effective)
    rf = np.maximum(ro - 0.5, 0.0)
    safe_rf = np.where(rf > 0.0, rf, 1.0)
    corr = (42.5 * rf * np.exp(-100.0 / (251.0 - mo))
            * (1.0 - np.exp(-6.93 / safe_rf)))
    extra = np.where(mo > 150.0,
                     0.0015 * np.maximum(mo - 150.0, 0.0) ** 2 * np.sqrt(rf), 0.0)
    mr = np.minimum(mo + corr + extra, 250.0)
    mo = np.where(ro > 0.5, mr, mo)

    Ed = (0.942 * H ** 0.679 + 11.0 * np.exp((H - 100.0) / 10.0)
          + 0.18 * (21.1 - T) * (1.0 - np.exp(-0.115 * H)))
    Ew = (0.618 * H ** 0.753 + 10.0 * np.exp((H - 100.0) / 10.0)
          + 0.18 * (21.1 - T) * (1.0 - np.exp(-0.115 * H)))
    ko = 0.424 * (1.0 - (H / 100.0) ** 1.7) + 0.0694 * np.sqrt(W) * (1.0 - (H / 100.0) ** 8)
    kd = ko * 0.581 * np.exp(0.0365 * T)
    m_dry = Ed + (mo - Ed) * 10.0 ** (-kd)
    kl = (0.424 * (1.0 - ((100.0 - H) / 100.0) ** 1.7)
          + 0.0694 * np.sqrt(W) * (1.0 - ((100.0 - H) / 100.0) ** 8))
    kw = kl * 0.581 * np.exp(0.0365 * T)
    m_wet = Ew - (Ew - mo) * 10.0 ** (-kw)
    m = np.where(mo > Ed, m_dry, np.where(mo < Ew, m_wet, mo))

    F = 59.5 * (250.0 - m) / (147.2 + m)
    return _scalar(np.clip(F, 0.0, 101.0))


def ffmc_to_moisture(ffmc_code):
    """FFMC code -> fine fuel moisture content (%, dry-weight scale ~0-250).

    The inverse of the FFMC scale; comparable to the 1-h dead fuel moisture, so an
    FFMC series can be overlaid on a dead-fuel-moisture panel. Scalar or array.
    """
    F = np.asarray(ffmc_code, dtype=float)
    return _scalar(147.2 * (101.0 - F) / (59.5 + F))


def dmc(temp_c, relative_humidity, rain_mm, dmc_prev, month, *,
        day_length=DMC_DAY_LENGTH):
    """Duff Moisture Code from noon weather and yesterday's DMC (0..inf, drier up)."""
    T = np.maximum(np.asarray(temp_c, dtype=float), -1.1)
    H = np.clip(np.asarray(relative_humidity, dtype=float), 0.0, 100.0)
    ro = np.maximum(np.asarray(rain_mm, dtype=float), 0.0)
    Po = np.maximum(np.asarray(dmc_prev, dtype=float), 0.0)
    Le = _month_factor(day_length, month)

    rk = 1.894 * (T + 1.1) * (100.0 - H) * Le * 1e-4
    rw = 0.92 * ro - 1.27
    wmi = 20.0 + np.exp(5.6348 - Po / 43.43)
    safe_Po = np.where(Po > 0.0, Po, 1.0)
    b = np.where(Po <= 33.0, 100.0 / (0.5 + 0.3 * Po),
                 np.where(Po <= 65.0, 14.0 - 1.3 * np.log(safe_Po),
                          6.2 * np.log(safe_Po) - 17.2))
    wmr = wmi + 1000.0 * rw / (48.77 + b * rw)
    pr = 43.43 * (5.6348 - np.log(np.maximum(wmr - 20.0, 1e-6)))
    pr = np.where(ro > 1.5, np.maximum(pr, 0.0), Po)
    return _scalar(np.maximum(pr + rk, 0.0))


def dc(temp_c, rain_mm, dc_prev, month, *, day_length=DC_DAY_LENGTH):
    """Drought Code from noon temperature and yesterday's DC (deep-drought memory)."""
    T = np.maximum(np.asarray(temp_c, dtype=float), -2.8)
    ro = np.maximum(np.asarray(rain_mm, dtype=float), 0.0)
    Do = np.maximum(np.asarray(dc_prev, dtype=float), 0.0)
    Lf = _month_factor(day_length, month)

    pe = np.maximum((0.36 * (T + 2.8) + Lf) / 2.0, 0.0)
    rd = 0.83 * ro - 1.27
    Qo = 800.0 * np.exp(-Do / 400.0)
    Qr = Qo + 3.937 * rd
    Dr = np.maximum(400.0 * np.log(800.0 / np.maximum(Qr, 1e-6)), 0.0)
    Dr = np.where(ro > 2.8, Dr, Do)
    return _scalar(np.maximum(Dr + pe, 0.0))


def isi(wind_kmh, ffmc_code):
    """Initial Spread Index from wind and the current FFMC (spread-rate proxy)."""
    W = np.maximum(np.asarray(wind_kmh, dtype=float), 0.0)
    m = ffmc_to_moisture(ffmc_code)
    m = np.asarray(m, dtype=float)
    f_wind = np.exp(0.05039 * W)
    f_f = 91.9 * np.exp(-0.1386 * m) * (1.0 + m ** 5.31 / 4.93e7)
    return _scalar(0.208 * f_wind * f_f)


def bui(dmc_code, dc_code):
    """Buildup Index from DMC and DC (fuel available to a spreading fire)."""
    P = np.maximum(np.asarray(dmc_code, dtype=float), 0.0)
    D = np.maximum(np.asarray(dc_code, dtype=float), 0.0)
    denom = P + 0.4 * D
    safe = np.where(denom > 0.0, denom, 1.0)
    low = 0.8 * P * D / safe
    high = P - (1.0 - 0.8 * D / safe) * (0.92 + (0.0114 * P) ** 1.7)
    U = np.where(P <= 0.4 * D, low, high)
    U = np.where(denom > 0.0, U, 0.0)
    return _scalar(np.maximum(U, 0.0))


def fwi(isi_val, bui_val):
    """Final Fire Weather Index from ISI and BUI (general fire-intensity rating)."""
    I = np.maximum(np.asarray(isi_val, dtype=float), 0.0)
    U = np.maximum(np.asarray(bui_val, dtype=float), 0.0)
    f_d = np.where(U <= 80.0, 0.626 * U ** 0.809 + 2.0,
                   1000.0 / (25.0 + 108.64 * np.exp(-0.023 * U)))
    B = 0.1 * I * f_d
    # The log-power branch is only taken for B > 1 (where 0.434*ln B > 0); clamp the
    # base so the discarded B <= 1 side of np.where doesn't raise on a negative power.
    base = np.maximum(0.434 * np.log(np.where(B > 1.0, B, np.e)), 0.0)
    S = np.where(B > 1.0, np.exp(2.72 * base ** 0.647), B)
    return _scalar(np.maximum(S, 0.0))


def daily_severity_rating(fwi_val):
    """Daily Severity Rating from FWI (a more linear measure of control effort)."""
    return _scalar(0.0272 * np.asarray(fwi_val, dtype=float) ** 1.77)


# --- stateful daily stepper ---------------------------------------------------

@dataclass
class FWIState:
    """Carried moisture codes (yesterday's values); the FWI System's memory."""
    ffmc: float = STARTUP_FFMC
    dmc: float = STARTUP_DMC
    dc: float = STARTUP_DC


@dataclass
class FWIIndices:
    """All six FWI codes for one day, plus fine fuel moisture and severity."""
    ffmc: float
    dmc: float
    dc: float
    isi: float
    bui: float
    fwi: float
    fine_fuel_moisture: float          # % (dry weight), from FFMC
    daily_severity_rating: float


class FWISystem:
    """Steps the FWI codes day by day, carrying the moisture memory between days.

    ``hemisphere="south"`` shifts the day-length tables six months. Seed
    :class:`FWIState` from a prior run / overwintering when you have it; otherwise
    the spring-startup defaults apply (indicative for FFMC/ISI on day one).
    """

    def __init__(self, state: FWIState | None = None, *, hemisphere: str = "north"):
        self.state = state or FWIState()
        if hemisphere == "south":
            self._le, self._lf = southern_hemisphere_tables()
        else:
            self._le, self._lf = DMC_DAY_LENGTH, DC_DAY_LENGTH

    def step(self, *, temperature, relative_humidity, wind_kmh, rain_mm,
             month) -> FWIIndices:
        """Advance one day from noon weather; updates and returns the codes."""
        F = ffmc(temperature, relative_humidity, wind_kmh, rain_mm, self.state.ffmc)
        P = dmc(temperature, relative_humidity, rain_mm, self.state.dmc, month,
                day_length=self._le)
        D = dc(temperature, rain_mm, self.state.dc, month, day_length=self._lf)
        I = isi(wind_kmh, F)
        U = bui(P, D)
        S = fwi(I, U)
        self.state = FWIState(ffmc=F, dmc=P, dc=D)
        return FWIIndices(ffmc=F, dmc=P, dc=D, isi=I, bui=U, fwi=S,
                          fine_fuel_moisture=ffmc_to_moisture(F),
                          daily_severity_rating=daily_severity_rating(S))


def fwi_from_records(records, *, state: FWIState | None = None,
                     hemisphere: str = "north"):
    """Run the FWI System over a sequence of daily-noon weather records.

    Each record is a mapping with ``temperature`` (C), ``relative_humidity`` (%),
    ``rain_mm``, ``month`` (1-12), and either ``wind_kmh`` or ``wind_ms``. Returns
    a list of :class:`FWIIndices`, one per record (codes carry day to day).
    """
    sys = FWISystem(state=state, hemisphere=hemisphere)
    out = []
    for r in records:
        w = r.get("wind_kmh")
        if w is None:
            w = wind_ms_to_kmh(r["wind_ms"])
        out.append(sys.step(
            temperature=r["temperature"], relative_humidity=r["relative_humidity"],
            wind_kmh=w, rain_mm=r.get("rain_mm", 0.0), month=r["month"]))
    return out


def spinup_state(records, *, state: FWIState | None = None,
                 hemisphere: str = "north") -> FWIState:
    """Run the FWI System over a spin-up history and return the *carried codes*.

    Feed a ~3-4 week series of daily-noon ``records`` (as in :func:`fwi_from_records`)
    ending the day before ignition; the returned :class:`FWIState` holds the
    FFMC/DMC/DC that have integrated those weeks of weather -- the drought memory
    (high DC after a dry spell) that a single day off the spring-startup defaults
    cannot capture. Use it as ``previous`` for :func:`gridded_fwi` on the fire day.
    """
    series = fwi_from_records(records, state=state, hemisphere=hemisphere)
    if not series:
        return state or FWIState()
    last = series[-1]
    return FWIState(ffmc=last.ffmc, dmc=last.dmc, dc=last.dc)


def gridded_fwi(*, temperature, relative_humidity, wind_kmh, month,
                previous: FWIState, rain_mm=0.0, hemisphere: str = "north") -> dict:
    """Per-cell FWI codes for one day from gridded noon weather + a prior state.

    ``temperature`` / ``relative_humidity`` / ``wind_kmh`` (and optional
    ``rain_mm``) are scalars or same-shaped arrays (e.g. an AOI raster from
    :meth:`GriddedAtmosphere.field_on`); ``previous`` is the spun-up scalar
    :class:`FWIState` (:func:`spinup_state`) whose DMC/DC are broadcast across the
    grid -- deep-drought codes are near-uniform over a small domain while FFMC/ISI
    vary with the local weather. Returns a dict of same-shaped code arrays
    (``ffmc``/``dmc``/``dc``/``isi``/``bui``/``fwi``/``fine_fuel_moisture``).
    """
    le, lf = (southern_hemisphere_tables() if hemisphere == "south"
              else (DMC_DAY_LENGTH, DC_DAY_LENGTH))
    F = ffmc(temperature, relative_humidity, wind_kmh, rain_mm, previous.ffmc)
    P = dmc(temperature, relative_humidity, rain_mm, previous.dmc, month, day_length=le)
    D = dc(temperature, rain_mm, previous.dc, month, day_length=lf)
    I = isi(wind_kmh, F)
    U = bui(P, D)
    S = fwi(I, U)
    return {"ffmc": F, "dmc": P, "dc": D, "isi": I, "bui": U, "fwi": S,
            "fine_fuel_moisture": ffmc_to_moisture(F)}


def gridded_fwi_on(field_state, *, month, previous: FWIState, rain_mm=0.0,
                   hemisphere: str = "north") -> dict:
    """:func:`gridded_fwi` from a gridded :class:`~pyflam.atmosphere.AtmosphericState`.

    ``field_state`` is a per-cell state (e.g. ``provider.field_on(ls, noon)``);
    wind is converted from m/s. See :func:`gridded_fwi` for ``previous`` / returns.
    """
    return gridded_fwi(
        temperature=field_state.temperature,
        relative_humidity=field_state.relative_humidity,
        wind_kmh=wind_ms_to_kmh(field_state.wind_speed),
        month=month, previous=previous, rain_mm=rain_mm, hemisphere=hemisphere)


def fwi_from_provider(provider, latitude, longitude, days, *, rain_mm=0.0,
                      hour_utc=12, state: FWIState | None = None,
                      hemisphere: str = "north"):
    """FWI series at a point from an ERA5 / GFS / ICON provider.

    Samples each day's noon (``hour_utc``) state from ``provider.state_at`` for the
    ``days`` (an iterable of :class:`datetime.date`/``datetime``), takes T / RH /
    wind from it, and combines with ``rain_mm`` -- a scalar (same each day) or a
    per-day sequence of 24 h noon-to-noon totals (0 for a dry spell; ERA5 ``tp`` /
    GFS ``apcp`` summed noon-to-noon otherwise). Returns a list of
    :class:`FWIIndices`.
    """
    days = list(days)
    rains = ([rain_mm] * len(days) if np.ndim(rain_mm) == 0 else list(rain_mm))
    records = []
    for d, rr in zip(days, rains):
        when = datetime(d.year, d.month, d.day, hour_utc)
        st = provider.state_at(latitude, longitude, when)
        records.append(dict(
            temperature=st.temperature, relative_humidity=st.relative_humidity,
            wind_ms=st.wind_speed, rain_mm=rr, month=d.month))
    return fwi_from_records(records, state=state, hemisphere=hemisphere)
