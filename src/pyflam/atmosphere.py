# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (c) 2026 Cristiano Foderi <cristiano.foderi@gmail.com>

"""Atmospheric forcing for fire simulations: weather/reanalysis -> fire inputs.

This module lets a pyflam run be driven by real atmospheric data instead of fixed
scenario inputs -- a forecast (GFS, HRRR, WRF output) for **near-real-time** runs,
or a reanalysis (ERA5 from Copernicus) for **re-analysis** runs. It does three
things:

1. **Abstracts the data source.** :class:`AtmosphereProvider` returns an
   :class:`AtmosphericState` for a location and time. Concrete providers wrap a
   constant (testing), a gridded dataset (any xarray/NetCDF/GRIB -- WRF, ERA5,
   GFS) via a variable-name map, or the ERA5/GFS services (lazy, documented).

2. **Carries the fire-relevant variables, with a focus on convection.** Beyond
   the surface state (10 m wind, 2 m temperature/humidity, pressure) the state
   holds the **convective / energy-flux** fields that govern fire-atmosphere
   coupling: surface sensible/latent heat flux, CAPE, CIN, boundary-layer height
   and stability. These set the background buoyancy the fire's own plume develops
   into (see :mod:`pyflam.pyroconvection`).

3. **Derives pyflam inputs from the state** -- physics that is fully testable
   offline: dead fuel moisture from temperature/humidity (NFDRS equilibrium
   moisture content), midflame wind, Monin-Obukhov stability from the surface
   heat flux, an ambient ground heat flux for the buoyant CFD, and a convective
   plume factor (CAPE/stability) that strengthens lofting and spotting in an
   unstable atmosphere.

Network/file providers need optional deps (``xarray``/``cfgrib``/``cdsapi``);
they are imported lazily with a clear message. Everything else is pure NumPy.

References:
    Simard, A.J. 1968. The moisture content of forest fuels. (NFDRS EMC.)
    Stull, R.B. 1988. An Introduction to Boundary Layer Meteorology. (Obukhov L.)
    Hersbach, H. et al. 2020. The ERA5 global reanalysis. QJRMS 146.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from datetime import datetime, timezone

import numpy as np

from .units import (
    celsius_to_fahrenheit, kelvin_to_celsius, m_per_s_to_ft_per_min,
)

_G = 9.81
_RHO_AIR = 1.2
_CP_AIR = 1005.0
_KAPPA = 0.4
# Reference CAPE (J/kg) at which the convective plume factor reaches ~2x; a
# moderately unstable atmosphere. Tunable.
_CAPE_REF = 1500.0


@dataclass
class AtmosphericState:
    """Fire-relevant atmospheric state at a point and time.

    Surface variables drive the standard fire model; the convective / energy-flux
    variables drive the fire-atmosphere (plume) coupling. Missing fields are
    ``None`` and the derivations fall back sensibly.

    Units: wind m/s (``direction`` deg FROM, met); temperature degrees C;
    relative humidity %; pressure hPa; heat fluxes W/m^2 (positive = surface
    heating the air); CAPE/CIN J/kg; boundary-layer height m.
    """

    wind_speed: float                        # 10 m wind speed (m/s)
    wind_direction: float                    # deg FROM (meteorological)
    temperature: float                       # 2 m air temperature (C)
    relative_humidity: float                 # 2 m RH (%)
    pressure: float = 1013.0                 # surface pressure (hPa)
    # --- convective / energy-flux ---
    sensible_heat_flux: float | None = None  # surface sensible heat flux (W/m^2)
    latent_heat_flux: float | None = None    # surface latent heat flux (W/m^2)
    cape: float | None = None                # convective available PE (J/kg)
    cin: float | None = None                 # convective inhibition (J/kg)
    boundary_layer_height: float | None = None   # PBL height (m)
    # --- metadata ---
    time: datetime | None = None
    latitude: float | None = None
    longitude: float | None = None

    @classmethod
    def from_si(cls, *, wind_u=None, wind_v=None, wind_speed=None,
                wind_direction=None, temperature_K=None, temperature_C=None,
                dewpoint_K=None, relative_humidity=None, **kw):
        """Build a state from common forecast/reanalysis SI variables.

        Accepts wind as (u, v) components or speed/direction; temperature in K or
        C; humidity as RH(%) or from a dewpoint (K). Extra keyword fields pass
        straight through (e.g. ``cape=...``, ``sensible_heat_flux=...``).
        """
        def _scalarize(v):
            v = np.asarray(v)
            return float(v) if v.ndim == 0 else v

        if wind_u is not None and wind_v is not None:
            wind_speed = _scalarize(np.hypot(wind_u, wind_v))
            # meteorological FROM-direction of the (u eastward, v northward) wind.
            wind_direction = _scalarize(
                np.degrees(np.arctan2(-np.asarray(wind_u), -np.asarray(wind_v))) % 360.0)
        t_c = (temperature_C if temperature_C is not None
               else kelvin_to_celsius(np.asarray(temperature_K)))
        t_c = _scalarize(t_c)
        if relative_humidity is None and dewpoint_K is not None:
            relative_humidity = relative_humidity_from_dewpoint(
                t_c, kelvin_to_celsius(np.asarray(dewpoint_K)))
        return cls(wind_speed=wind_speed, wind_direction=wind_direction,
                   temperature=t_c, relative_humidity=relative_humidity, **kw)


# --- humidity / fuel moisture -------------------------------------------------

def relative_humidity_from_dewpoint(temp_c, dewpoint_c):
    """Relative humidity (%) from temperature and dewpoint (C), Magnus formula.

    Scalars or arrays.
    """
    def es(t):
        t = np.asarray(t, dtype=float)
        return 6.112 * np.exp(17.67 * t / (t + 243.5))
    rh = np.clip(100.0 * es(dewpoint_c) / es(temp_c), 0.0, 100.0)
    return float(rh) if rh.ndim == 0 else rh


def equilibrium_moisture_content(temp_c, relative_humidity):
    """NFDRS equilibrium moisture content (% , Simard 1968).

    The moisture a fine dead fuel equilibrates to at temperature ``temp_c`` (C)
    and ``relative_humidity`` (%). Piecewise in RH; temperature in F internally.
    Works on scalars or NumPy arrays (returns the matching type).
    """
    h = np.clip(np.asarray(relative_humidity, dtype=float), 0.0, 100.0)
    t = celsius_to_fahrenheit(np.asarray(temp_c, dtype=float))
    low = 0.03229 + 0.281073 * h - 0.000578 * h * t
    mid = 2.22749 + 0.160107 * h - 0.014784 * t
    high = 21.0606 + 0.005565 * h * h - 0.00035 * h * t - 0.483199 * h
    emc = np.where(h < 10.0, low, np.where(h <= 50.0, mid, high))
    emc = np.maximum(emc, 0.0)
    return float(emc) if emc.ndim == 0 else emc


# Standard dead-fuel response time lags (hours): the larger the fuel, the slower
# it tracks the equilibrium moisture (its "memory").
FUEL_TIME_LAGS = {"m_1h": 1.0, "m_10h": 10.0, "m_100h": 100.0}


def time_lag_step(m_prev, emc_fraction, dt_hours, tau_hours):
    """Advance a fuel moisture toward equilibrium over ``dt_hours`` (time-lag law).

    ``m(t+dt) = EMC + (m(t) - EMC) exp(-dt/tau)`` -- exponential approach to the
    equilibrium with the fuel's response time ``tau``. Scalars or arrays.
    """
    return emc_fraction + (np.asarray(m_prev, dtype=float) - emc_fraction) \
        * np.exp(-float(dt_hours) / float(tau_hours))


@dataclass
class DeadFuelMoistureModel:
    """Stateful time-lag dead fuel moisture (the operational Nelson-type model).

    Instantaneous EMC ignores that larger fuels lag the weather; this carries the
    1/10/100-h moistures and steps each toward the current EMC by its own time lag
    (Fosberg & Deeming 1971; cf. Nelson 2000 for the full diffusion model). Use it
    in a time-marched run so fuels remember recent humidity rather than snapping to
    the latest value -- key for diurnal drying/recovery and reanalysis runs.

    Moistures are fractions.
    """

    m_1h: float
    m_10h: float
    m_100h: float

    @classmethod
    def equilibrium(cls, state: "AtmosphericState") -> "DeadFuelMoistureModel":
        """Initialise every class at the current equilibrium moisture content."""
        emc = equilibrium_moisture_content(
            state.temperature, state.relative_humidity) / 100.0
        return cls(m_1h=emc, m_10h=emc, m_100h=emc)

    def update(self, state: "AtmosphericState", dt_minutes: float) -> dict:
        """Step the three classes toward the state's EMC and return the moistures.

        Scalars stay Python floats; a per-cell ``state`` (array temperature /
        relative humidity, e.g. from :meth:`GriddedAtmosphere.field_on`) steps the
        three classes per cell and returns arrays -- so a spatial march can carry
        gridded fuel-moisture memory.
        """
        emc = equilibrium_moisture_content(
            state.temperature, state.relative_humidity) / 100.0
        dt_h = dt_minutes / 60.0

        def _keep(x):                       # 0-d -> float, arrays stay arrays
            a = np.asarray(x)
            return float(a) if a.ndim == 0 else a

        self.m_1h = _keep(time_lag_step(self.m_1h, emc, dt_h, FUEL_TIME_LAGS["m_1h"]))
        self.m_10h = _keep(time_lag_step(self.m_10h, emc, dt_h, FUEL_TIME_LAGS["m_10h"]))
        self.m_100h = _keep(time_lag_step(self.m_100h, emc, dt_h,
                                          FUEL_TIME_LAGS["m_100h"]))
        return {"m_1h": self.m_1h, "m_10h": self.m_10h, "m_100h": self.m_100h}


def dead_fuel_moisture(state: AtmosphericState, *, offsets=(0.0, 1.0, 2.0)) -> dict:
    """Dead 1/10/100-h fuel moisture (fractions) from the atmospheric state.

    The 1-h moisture is the equilibrium moisture content; the slower 10-h and
    100-h classes are offset upward (``offsets`` in %, a standard simple scheme --
    a full Nelson dead-fuel model with lags is a refinement). Returns the
    ``m_1h/m_10h/m_100h`` kwargs that :func:`pyflam.spread` expects.
    """
    emc = equilibrium_moisture_content(state.temperature, state.relative_humidity)
    o1, o10, o100 = offsets
    return {
        "m_1h": (emc + o1) / 100.0,
        "m_10h": (emc + o10) / 100.0,
        "m_100h": (emc + o100) / 100.0,
    }


# --- stability / energy flux --------------------------------------------------

def friction_velocity(wind_speed, *, z0=0.1, height=10.0) -> float:
    """Neutral friction velocity u* (m/s) from a wind at ``height`` (m)."""
    return _KAPPA * float(wind_speed) / np.log((height + z0) / z0)


def obukhov_length(state: AtmosphericState, *, z0=0.1) -> float:
    """Monin-Obukhov length L (m) from the surface sensible heat flux.

    ``L < 0`` unstable (daytime convective), ``L > 0`` stable (nocturnal),
    ``inf`` neutral. Returns ``inf`` when no heat flux is available.
    """
    q = state.sensible_heat_flux
    if q is None or abs(q) < 1e-6:
        return float("inf")
    ustar = friction_velocity(state.wind_speed, z0=z0)
    t_k = state.temperature + 273.15
    return -(ustar ** 3) * _RHO_AIR * _CP_AIR * t_k / (_KAPPA * _G * q)


def stability_class(state: AtmosphericState, *, z0=0.1) -> str:
    """Coarse stability class: ``"unstable"`` / ``"neutral"`` / ``"stable"``.

    Uses the heat flux sign (and CAPE/CIN when present) -- the regime that decides
    whether a fire plume grows freely (unstable) or is capped (stable).
    """
    q = state.sensible_heat_flux
    if state.cape is not None and state.cape > 500.0 and (state.cin or 0.0) < 50.0:
        return "unstable"
    if q is None:
        return "neutral"
    if q > 10.0:
        return "unstable"
    if q < -10.0:
        return "stable"
    return "neutral"


def ambient_surface_heat_flux(state: AtmosphericState) -> float:
    """Ambient ground sensible heat flux (W/m^2) for the buoyant background.

    This is the *atmosphere's own* surface heating into which the fire's plume
    develops (the fire flux from :mod:`pyflam.pyroconvection` adds on top). Uses
    the reported sensible heat flux, or 0 (neutral) if unavailable. Scalar or
    array (gridded state).
    """
    q = state.sensible_heat_flux
    if q is None:
        return 0.0
    q = np.asarray(q, dtype=float)
    return float(q) if q.ndim == 0 else q


def convective_plume_factor(state: AtmosphericState, *,
                            profile: "AtmosphericProfile | None" = None) -> float:
    """Loft enhancement (>= 1) for the fire plume from atmospheric instability.

    An unstable atmosphere lets a fire plume rise higher and entrain less,
    strengthening lofting and spotting; a stable one caps it. Bounded heuristic --
    the convective coupling, not a cloud model.

    Without a ``profile`` this uses surface CAPE + Monin-Obukhov stability (the
    original behaviour). **With a vertical ``profile``** it shifts to the predictors
    the pyroCb literature favours over surface CAPE: a dry mixed layer (high LCL /
    large near-surface dewpoint depression) capped by moisture aloft (the inverted-V
    sounding) and elevated lower-tropospheric instability/dryness (Continuous
    Haines). Surface CAPE is a poor pyroCb predictor -- pyroCb routinely form with
    near-zero CAPE (Castellnou et al. 2022; Peterson et al. 2017; Mills & McCaw 2010).

    **Per cell.** The ``state`` fields may be scalars or 2D arrays (a gridded
    ``field_on`` state); the factor is then an array. The aloft terms (C-Haines,
    mid-level moisture) come from the single column ``profile`` and apply uniformly,
    but the inverted-V boost is realised **per cell via that cell's own surface
    dryness** -- so under a moist-aloft sounding, locally dry (high-LCL) cells get
    the pyroconvective boost and locally moist cells do not.
    """
    cape = np.asarray(state.cape if state.cape is not None else 0.0, dtype=float)
    factor = 1.0 + cape / _CAPE_REF

    # Vectorized stability multiplier (heat-flux sign, with a CAPE/CIN override).
    cin = np.asarray(state.cin if state.cin is not None else 0.0, dtype=float)
    unstable = (cape > 500.0) & (cin < 50.0)
    if state.sensible_heat_flux is not None:
        q = np.asarray(state.sensible_heat_flux, dtype=float)
        unstable = unstable | (q > 10.0)
        stable = (q < -10.0) & ~unstable
    else:
        stable = np.zeros_like(factor, dtype=bool)
    factor = factor * np.where(unstable, 1.2, np.where(stable, 0.7, 1.0))

    if profile is not None:
        ch = continuous_haines(profile)                  # column scalar
        factor = factor * (1.0 + 0.6 * np.clip((ch - 5.0) / 6.0, 0.0, 1.0))
        _, _, mid_rh = inverted_v(profile)               # moisture aloft (column)
        if mid_rh >= 50.0:
            # inverted-V realised per cell: dry surface (large dewpoint depression)
            # under the moist-aloft column. Ramps 1 -> 1.25 over a 10->30 C depression.
            depr = (np.asarray(state.temperature, dtype=float)
                    - dewpoint_from_rh(state.temperature, state.relative_humidity))
            factor = factor * (1.0 + 0.25 * np.clip((depr - 10.0) / 20.0, 0.0, 1.0))

    factor = np.clip(factor, 0.5, 3.0)
    return float(factor) if np.ndim(factor) == 0 else factor


# --- pyroconvection potential (vertical-profile diagnostics) -------------------
#
# High-convective wildfire is a *vertical* problem: whether a fire's plume merely
# rises or breaks through to pyrocumulus / pyrocumulonimbus (pyroCb) is set by the
# boundary-layer -> LCL -> free-convection geometry and by mid-level moisture, NOT
# by surface CAPE (pyroCb often form with little or no surface-based CAPE). The
# canonical pyroCb environment is a deep, dry, well-mixed boundary layer capped by
# a moist layer aloft -- the "inverted-V" sounding. These diagnostics distil that
# into numbers the fire model can act on.
#
# References:
#   Castellnou, M., et al. 2022. Pyroconvection classification based on atmospheric
#       vertical profiling. J. Geophys. Res. Atmos. 127, e2022JD036920.
#   Mills, G.A.; McCaw, W.L. 2010. Atmospheric stability environments and fire
#       weather in Australia -- the Continuous Haines index. CAWCR Tech. Rep. 20.
#   Peterson, D.A., et al. 2017. PyroCb climatology (mid-troposphere humidity as a
#       pyroCb discriminator).

def _es_hpa(temp_c):
    """Saturation vapour pressure (hPa) over water, Magnus form."""
    t = np.asarray(temp_c, dtype=float)
    return 6.112 * np.exp(17.67 * t / (t + 243.5))


def dewpoint_from_rh(temp_c, relative_humidity):
    """Dewpoint (C) from temperature (C) and RH (%), inverse Magnus. Scalars/arrays."""
    rh = np.clip(np.asarray(relative_humidity, dtype=float), 1e-3, 100.0)
    t = np.asarray(temp_c, dtype=float)
    gamma = np.log(rh / 100.0) + 17.67 * t / (t + 243.5)
    dp = 243.5 * gamma / (17.67 - gamma)
    return float(dp) if dp.ndim == 0 else dp


def lcl_height_m(temp_c, relative_humidity):
    """Lifting condensation level height above ground (m) -- Espy / Lawrence (2005).

    ``LCL ~ 125 * (T - Td)`` metres, with ``T - Td`` the dewpoint depression (C); a
    hot, dry mixed layer (large depression) pushes the LCL high, the inverted-V
    setup. Scalars or arrays.
    """
    td = dewpoint_from_rh(temp_c, relative_humidity)
    depression = np.maximum(np.asarray(temp_c, dtype=float) - td, 0.0)
    h = 125.0 * depression
    return float(h) if np.ndim(h) == 0 else h


@dataclass
class AtmosphericProfile:
    """A coarse vertical sounding: temperature + dewpoint at pressure levels.

    ``pressure`` (hPa, descending or ascending), ``temperature`` and ``dewpoint``
    (C) are equal-length 1-D arrays. This is the minimal profile the pyroconvection
    diagnostics (Continuous Haines, mid-level moisture) need beyond the surface
    state; build one from a sounding, a model column, or the standard pressure
    levels of a reanalysis. ``from_rh`` builds it from RH instead of dewpoint.
    """

    pressure: np.ndarray            # hPa
    temperature: np.ndarray         # C
    dewpoint: np.ndarray            # C

    @classmethod
    def from_rh(cls, pressure, temperature, relative_humidity):
        return cls(pressure=np.asarray(pressure, float),
                   temperature=np.asarray(temperature, float),
                   dewpoint=dewpoint_from_rh(temperature, relative_humidity))

    def _at(self, level_hpa):
        """Linear-interpolate (temperature, dewpoint) to a pressure level."""
        p = np.asarray(self.pressure, float)
        order = np.argsort(p)
        t = float(np.interp(level_hpa, p[order], np.asarray(self.temperature)[order]))
        d = float(np.interp(level_hpa, p[order], np.asarray(self.dewpoint)[order]))
        return t, d


def continuous_haines(profile: AtmosphericProfile, *,
                      low_hpa: float = 850.0, high_hpa: float = 700.0) -> float:
    """Continuous Haines index (C-Haines) -- lower-tropospheric stability + dryness.

    ``CA = (T_low - T_high)/2 - 2``  (a stability term: the 850->700 hPa lapse),
    ``CB = min((T_high - Td_high)/3 - 1, 5)`` capped, plus the >5 reduction
    (Mills & McCaw 2010). ``C-Haines = CA + CB``; higher = drier & more unstable
    aloft = greater pyroconvective / blow-up potential (typically 0-13). The level
    pair defaults to the Australian 850/700 hPa convention.
    """
    t_lo, _ = profile._at(low_hpa)
    t_hi, td_hi = profile._at(high_hpa)
    ca = (t_lo - t_hi) / 2.0 - 2.0
    cb = (t_hi - td_hi) / 3.0 - 1.0
    if cb > 5.0:                      # Mills & McCaw cap on the moisture term
        cb = 5.0 + (cb - 5.0) / 2.0
    return float(ca + cb)


def inverted_v(profile: AtmosphericProfile, *, mid_hpa: float = 600.0,
               surface_depression_c: float = 10.0, mid_rh_pct: float = 50.0):
    """Detect the inverted-V (dry mixed layer, moist aloft) pyroCb-prone sounding.

    Returns ``(is_inverted_v, surface_dewpoint_depression_c, mid_level_rh_pct)``.
    The signature is a **large near-surface dewpoint depression** (dry, deep mixed
    layer -> high LCL) together with **moister air aloft** (mid-level RH above
    ``mid_rh_pct``) -- the environment in which a fire plume that reaches the free
    convection level taps mid-level moisture and instability, even with little
    surface CAPE (Castellnou et al. 2022; Peterson et al. 2017).
    """
    p = np.asarray(profile.pressure, float)
    i_sfc = int(np.argmax(p))                       # highest pressure = surface
    sfc_depr = float(profile.temperature[i_sfc] - profile.dewpoint[i_sfc])
    t_mid, td_mid = profile._at(mid_hpa)
    mid_rh = float(np.clip(100.0 * _es_hpa(td_mid) / _es_hpa(t_mid), 0.0, 100.0))
    flag = (sfc_depr >= surface_depression_c) and (mid_rh >= mid_rh_pct)
    return bool(flag), sfc_depr, mid_rh


@dataclass
class PyroconvectionPotential:
    """Vertical-profile pyroconvection diagnostics for a fire-weather column."""

    lcl_height_m: float
    continuous_haines: float | None      # needs a profile
    inverted_v: bool | None
    mid_level_rh_pct: float | None
    mixed_layer_depth_m: float | None
    plume_dominated_favorable: bool      # deep dry ABL + (moist aloft | high C-Haines)
    notes: str


def pyroconvection_potential(state: AtmosphericState, *,
                             profile: AtmosphericProfile | None = None,
                             chaines_threshold: float = 9.0) -> PyroconvectionPotential:
    """Assess a column's potential for plume-dominated / pyroconvective fire.

    Combines what the surface ``state`` gives (LCL from T/RH; PBL depth) with the
    vertical ``profile`` when available (Continuous Haines, the inverted-V test).
    The headline flag ``plume_dominated_favorable`` fires when the boundary layer is
    deep and dry (high LCL) **and** there is moisture/instability aloft (inverted-V,
    or C-Haines past ``chaines_threshold``) -- the geometry, not surface CAPE, that
    the pyroCb literature identifies. Surface-only (no profile) returns the LCL and
    a coarser flag from PBL depth + LCL alone.
    """
    lcl = lcl_height_m(state.temperature, state.relative_humidity)
    if np.ndim(lcl) != 0:
        raise ValueError("pyroconvection_potential takes a single (scalar) state")
    lcl = float(lcl)
    pbl = state.boundary_layer_height
    deep_dry_abl = lcl >= 1500.0 or (pbl is not None and float(pbl) >= 2000.0)

    if profile is None:
        flag = bool(deep_dry_abl and lcl >= 2000.0)
        return PyroconvectionPotential(
            lcl_height_m=lcl, continuous_haines=None, inverted_v=None,
            mid_level_rh_pct=None,
            mixed_layer_depth_m=(float(pbl) if pbl is not None else None),
            plume_dominated_favorable=flag,
            notes="surface-only: no profile; flag from LCL + PBL depth")

    ch = continuous_haines(profile)
    iv, sfc_depr, mid_rh = inverted_v(profile)
    aloft = iv or ch >= chaines_threshold
    flag = bool(deep_dry_abl and aloft)
    notes = (f"C-Haines {ch:.1f}; inverted-V {iv}; sfc dewpoint depression "
             f"{sfc_depr:.1f}C; mid RH {mid_rh:.0f}%")
    return PyroconvectionPotential(
        lcl_height_m=lcl, continuous_haines=ch, inverted_v=iv,
        mid_level_rh_pct=mid_rh,
        mixed_layer_depth_m=(float(pbl) if pbl is not None else None),
        plume_dominated_favorable=flag, notes=notes)


# --- Briggs bent-over plume rise + pyroCb firepower threshold ------------------
#
# How high a fire's plume rises (and whether it reaches the level where free moist
# convection can develop) follows the Briggs (1969) buoyant-plume framework: the
# observed wildfire plume is the archetypal bent-over plume in a crosswind (Lareau
# & Clements 2017). Inverting that plume against a single sounding gives the minimum
# firepower for pyroCb -- the PyroCb Firepower Threshold (Tory & Kepert 2021). The
# forms below are the standard analytical Briggs solutions; the PFT here is a
# simplified, documented version (reach the LCL against the capping stability), not
# the full Tory & Kepert plume-condensation integration.

# Reference air properties for the buoyancy-flux conversion.
_T_REF_K = 293.0                 # K
_PI = math.pi


def briggs_buoyancy_flux(heat_flux_w: float, *, temperature_k: float = _T_REF_K) -> float:
    """Plume buoyancy flux F (m^4 s^-3) from total convective heat flux (W).

    ``F = g * Q_c / (pi * rho * cp * T)`` -- the Briggs buoyancy flux for a heat
    source of power ``Q_c`` (the convective fraction of the firepower).
    """
    return (_G * float(heat_flux_w)) / (_PI * _RHO_AIR * _CP_AIR * float(temperature_k))


def brunt_vaisala_squared(profile: AtmosphericProfile,
                          low_hpa: float = 700.0, high_hpa: float = 500.0) -> float:
    """Static-stability parameter s = (g/theta) dtheta/dz (s^-2) over a layer.

    Positive in a stable layer (the moist cap a pyroCb plume must punch through),
    ~0 in a well-mixed layer. Heights come from the hypsometric approximation.
    """
    t_lo, _ = profile._at(low_hpa)
    t_hi, _ = profile._at(high_hpa)
    th_lo = (t_lo + 273.15) * (1000.0 / low_hpa) ** 0.286
    th_hi = (t_hi + 273.15) * (1000.0 / high_hpa) ** 0.286
    # Approx layer thickness (m): hypsometric, mean T.
    tbar = 0.5 * (t_lo + t_hi) + 273.15
    dz = 287.0 * tbar / _G * math.log(low_hpa / high_hpa)
    if dz <= 0.0:
        return 0.0
    dtheta_dz = (th_hi - th_lo) / dz
    return float(max(_G / (0.5 * (th_lo + th_hi)) * dtheta_dz, 0.0))


def briggs_plume_rise(heat_flux_w: float, wind_ms: float, *,
                      stability_s2: float = 0.0, distance_m: float | None = None,
                      temperature_k: float = _T_REF_K) -> float:
    """Bent-over plume rise (m) above the source -- Briggs (1969).

    Stable layer (``stability_s2`` > 0): the **final** rise
    ``dh = 2.6 (F / (U s))^(1/3)``. Neutral (s = 0) with a downwind ``distance_m``:
    the transitional rise ``dh = 1.6 F^(1/3) x^(2/3) / U`` (a neutral bent-over
    plume has no final height, so a distance is required). ``wind_ms`` is the
    cross-plume wind; very light winds are floored to keep the bent-over scaling
    valid.
    """
    f = briggs_buoyancy_flux(heat_flux_w, temperature_k=temperature_k)
    u = max(float(wind_ms), 0.5)
    if stability_s2 > 1e-9:
        return float(2.6 * (f / (u * stability_s2)) ** (1.0 / 3.0))
    if distance_m is None:
        raise ValueError("neutral plume rise needs distance_m (no final height)")
    return float(1.6 * f ** (1.0 / 3.0) * float(distance_m) ** (2.0 / 3.0) / u)


def pyrocb_firepower_threshold(state: AtmosphericState, profile: AtmosphericProfile,
                               *, convective_fraction: float = 0.6,
                               target_height_m: float | None = None) -> float:
    """Minimum firepower (W) for the plume to reach free moist convection -- a PFT.

    The PyroCb Firepower Threshold (Tory & Kepert 2021): the least fire power that,
    in *this* atmosphere, lifts the plume to where free moist convection can begin.
    Here that target height defaults to the **LCL** (from the surface state) -- the
    plume must reach it to condense and tap the moist instability aloft -- and the
    plume must rise against the static stability of the capping layer
    (:func:`brunt_vaisala_squared`). Inverting the stable bent-over Briggs rise
    ``dh = 2.6 (F/(U s))^(1/3)`` for the heat flux gives
    ``Q_c = (dh/2.6)^3 * U * s * pi rho cp T / g`` and ``firepower = Q_c /
    convective_fraction``. Lower threshold = more easily pyroconvective.

    The threshold **rises with the LCL height and with the cap stability** (a
    drier surface or a stronger inversion needs a more powerful fire). It is
    **simplified** vs the full Tory & Kepert plume-condensation integration, and is
    only meaningful **paired with** :func:`pyroconvection_potential`: a *low*
    threshold in a moist, stable column does not mean pyroCb -- you also need the
    dry, unstable, moist-aloft (inverted-V) environment for deep convection to
    follow. Returns ``inf`` only for a degenerate (non-positive) target height.
    """
    h = float(target_height_m if target_height_m is not None
              else lcl_height_m(state.temperature, state.relative_humidity))
    s = brunt_vaisala_squared(profile)
    if h <= 0.0:
        return math.inf
    u = max(float(state.wind_speed), 0.5)
    f_req = (h / 2.6) ** 3 * u * max(s, 0.0)                # required buoyancy flux
    q_c = f_req * _PI * _RHO_AIR * _CP_AIR * _T_REF_K / _G  # convective heat flux (W)
    return float(q_c / max(convective_fraction, 1e-6))


# --- pyroconvection TYPE classification (Castellnou et al. 2022) ---------------
#
# Beyond a binary "pyroCb-favourable" flag, Castellnou et al. (2022, JGR-Atmos
# 127, e2022JD036920) classify pyroconvection into an ordered ladder of prototypes
# (their Fig. 7) set by vertical-profile diagnostics. This reproduces that ladder
# so a column (or a grid of columns) can be labelled by *type*, the basis of the
# operational "tipus de piroconvecció" forecast product. The conditioning
# variables and thresholds are taken directly from the paper:
#
#   * mixed-layer (ABL) stability via the potential-temperature gradient dtheta/dz:
#       unstable < 1.0e-4 K/m ; neutral 1.0e-4..1.1e-3 ; stable > 1.1e-3 K/m
#       (after Liu & Liang 2010). Stable ML -> non-pyroCu plume only.
#   * the LCL/ABL height ratio: > 1 -> overshooting (brief) ; < 1 -> resilient.
#   * the upper-layer stability gamma-theta (free-troposphere dtheta/dz): a weak
#       cap (low gamma-theta) lets a resilient/overshooting pyroCu deepen to
#       pyroCb. The paper's cases bracket the boundary: M11 gamma-theta=4.2e-3
#       (resilient, not deep) vs SCQ51 gamma-theta=3.9e-3 (deep pyroCb).
#   * a fire-power gate: pyroCu coincided with fireline intensity > 1e4 kW/m
#       (Tedim et al. 2018) -- below that the column only supports a surface plume.

# Ordered low -> high pyroconvective activity, with an ordinal level and a
# fire-weather colour ramp (white -> green -> yellow -> orange -> dark red).
PYROCONVECTION_TYPES = (
    "surface_plume", "convection_plume", "overshooting_pyrocu",
    "resilient_pyrocu", "deep_pyrocu_pyrocb",
)
PYROCONVECTION_TYPE_LEVEL = {t: i for i, t in enumerate(PYROCONVECTION_TYPES)}
PYROCONVECTION_TYPE_COLOR = {
    "surface_plume": "#ffffff", "convection_plume": "#a6d96a",
    "overshooting_pyrocu": "#fee08b", "resilient_pyrocu": "#f46d43",
    "deep_pyrocu_pyrocb": "#7f0000",
}
# The Catalan Bombers (GRAF) operational palette, sampled from their published
# "Tipus de piroconvecció" legend bar (ICON-EU 26 Jul 2026 00Z) and rounded to the CSS named
# colours it is built from: darkred / gold / yellow / lightgreen. Offered so a forecaster can
# read the two products side by side without re-learning the key -- the class *names*
# correspond one-to-one, only the integers run opposite (GRAF numbers pyroCb 1, pyflam 4).
#
# Their scale has no surface-plume class, so white is retained for it.
#
# Trade-off worth knowing before choosing this: gold and yellow are adjacent hues, so the
# resilient/overshooting pair is markedly harder to tell apart than in the default
# (ColorBrewer RdYlGn) palette, where they are orange-red against buff. That cost is real --
# it is why the visual reading of GRAF's own maps cannot reliably separate those two classes.
PYROCONVECTION_TYPE_COLOR_GRAF = {
    "surface_plume": "#ffffff", "convection_plume": "#90ee90",
    "overshooting_pyrocu": "#ffff00", "resilient_pyrocu": "#ffd700",
    "deep_pyrocu_pyrocb": "#8b0000",
}
PYROCONVECTION_PALETTES = {"pyflam": PYROCONVECTION_TYPE_COLOR,
                           "graf": PYROCONVECTION_TYPE_COLOR_GRAF}


def pyroconvection_colors(palette: str = "pyflam") -> dict:
    """Class -> hex colour map. ``palette`` is ``"pyflam"`` (default) or ``"graf"``."""
    try:
        return PYROCONVECTION_PALETTES[palette]
    except KeyError:
        raise ValueError(f"unknown palette {palette!r}; "
                         f"expected one of {sorted(PYROCONVECTION_PALETTES)}") from None
# Human-readable English labels for plots/reports (keys = PYROCONVECTION_TYPES).
PYROCONVECTION_TYPE_LABEL = {
    "surface_plume": "Surface plume",
    "convection_plume": "Convection plume",
    "overshooting_pyrocu": "Overshooting pyroCu",
    "resilient_pyrocu": "Resilient pyroCu",
    "deep_pyrocu_pyrocb": "Deep pyroCu / pyroCb",
}

# Thresholds (K/m) from Castellnou et al. 2022, §2.4.1 / §3.3 / Table 1.
_ML_UNSTABLE = 1.0e-4
_ML_STABLE = 1.1e-3
_GAMMA_DEEP = 4.0e-3            # gamma-theta below this -> deep pyroCu/pyroCb
_FLI_PYROCU_KW = 1.0e4         # fireline intensity gate for pyroCu (Tedim 2018)

# Dry-air thermodynamics for the profile diagnostics. (_KAPPA above is von Karman,
# hence _KAPPA_DRY for the Poisson exponent Rd/cp.)
_RD = 287.05                   # J kg-1 K-1, gas constant for dry air
_CP_DRY = 1004.0               # J kg-1 K-1
_KAPPA_DRY = _RD / _CP_DRY     # ~0.286
_EPSILON = 0.622               # Rd/Rv

# Bulk-Richardson ABL + profile diagnostics.
_RIB_CRITICAL = 0.33           # Rib crossing that marks the ABL top (Zhang et al. 2014)
_PARCEL_EXCESS_K = 0.5         # theta excess defining the parcel mixing depth
# Start height (AGL) for the Rib search. Zhang et al. (2014) recommend 200 m, to clear the
# surface layer without discarding a shallow ABL. Castellnou et al. (2022) sec.2.4.2 raise it
# to ~400 m for their *in-plume* sondes, which are launched into the indraft converging on the
# plume base (Charland & Clements 2013) -- a fire-specific correction. Ambient forecast columns
# carry no indraft, so this path keeps Zhang's 200 m. The choice is not cosmetic: on the
# 2026-07-26 00Z run it moves the median 18Z Tuscany ABL from 225 m (200 m start) to 404 m
# (400 m start), though it barely changes how many columns are classifiable.
_RIB_START_M = 200.0
_ABL_MIN_M = 150.0             # plausibility clip on the diagnosed ABL depth
_ABL_MAX_M = 4500.0
# A shear *height* cannot be located from a handful of pressure levels: the
# moving-window maximum then just reflects the interpolation, not the flow. Below
# this level count the adaptive classifier drops the shear test rather than
# fabricate one. ICON-2I open data publishes 6 levels, so it takes the no-shear
# path; ERA5 (19+ levels) and radiosondes take the full path.
_SHEAR_MIN_LEVELS = 10


@dataclass(frozen=True)
class PyroconvThresholds:
    """Thresholds for the profile-based pyroconvection ladder.

    **Every field is annotated with its provenance**, in one of three grades, because the
    defensibility of a classification is the defensibility of its weakest threshold and that
    should be visible without archaeology:

    * ``[PRIMARY]``  traceable to a peer-reviewed source, with the citation.
    * ``[DERIVED]``  interpolated from a small number of published cases -- directionally
      supported, but the exact value is a choice.
    * ``[UNSOURCED]`` inherited from the reference (MARI) implementation with no citation
      found in Castellnou et al. (2022) or its references. Retained because removing them
      would change the ladder's structure, but they are **not** evidence and should not be
      argued from.

    No field is fitted to another group's product. One was -- ``rh_top_moist`` at 60 %, tuned
    so the deep-pyroCb fraction matched the Catalan Bombers ICON-EU output for 2026-07-15 --
    and it has been removed, because a threshold calibrated to a third party's forecast cannot
    then be used as independent corroboration of it. See ``rh_top_moist`` below.
    """

    # [PRIMARY] Castellnou et al. (2022) sec. 2.4.1, after Liu & Liang (2010): unstable
    # < 0.1e-3, stable > 1.1e-3, neutral between. Neutral and unstable are both pyroCu-capable,
    # so 1.1e-3 is the gate. (The 1.0e-3 in the reference implementation is a transcription
    # slip -- the paper says 1.1e-3.)
    ml_stable: float = 1.1e-3

    # [UNSOURCED] No counterpart in Castellnou et al. (2022), which admits overshooting on a
    # non-stable mixed layer without an upper bound on the gradient. Inherited from the MARI
    # port.
    ml_overshoot_max: float = 3.0e-3

    # [DERIVED] The paper brackets the deepening cap with two cases: M11 at 4.2e-3 (resilient,
    # did not deepen) against SCQ51 at 3.9e-3 (deep pyroCb). Two points, so the boundary is
    # located to within about that interval and no better.
    gamma_weak_cap: float = 4.2e-3
    # [DERIVED] Upper bracket from SCQ41 at 5.1e-3 (strong cap, inhibited). The 4.8e-3 used
    # here sits inside the gap between that case and gamma_weak_cap; the paper does not give it.
    gamma_strong_cap: float = 4.8e-3

    # [PRIMARY] Castellnou et al. (2022) sec. 2.4.1 conditions on LCL/ABL *above or below 1*:
    # "Values >1 or <1 point to turbulence above or below ABL."
    lcl_ratio_resilient_max: float = 1.00
    # [UNSOURCED] The paper gives no upper bound on the overshooting band, nor a separate
    # ratio ceiling for the deep class. Both are MARI-port values. They matter: together they
    # decide how far above the condensation level a column may sit and still be classified.
    lcl_ratio_overshoot_max: float = 1.60
    lcl_ratio_deep_max: float = 1.10

    # [UNSOURCED] The paper requires the shear maximum to sit *near* the ABL/LCL for the top
    # class but gives no numeric distance; 0.30 is the MARI port's.
    shear_distance_deep: float = 0.30

    # [UNSOURCED, was FITTED] RH at the ABL top required for classes 3-4.
    #
    # Restored to the reference (MARI) port's 80 % on 2026-07-28. It had been lowered to 60 %,
    # and the justification recorded in this file was that 60 % "brings the deep-pyroCb
    # fraction into agreement" with the Catalan Bombers ICON-EU product for a single day
    # (2026-07-15, 3.0 % against their 3.4 % at 15Z). That is a fit to a third party's
    # forecast on one day, and it sits on the gate for the two classes that matter, so any
    # later agreement with that product was partly manufactured. Removed for that reason, not
    # because 80 % is better supported -- it is not, it is equally uncited.
    #
    # The consequence is real and should be seen rather than tuned away: at 80 % only a small
    # fraction of well-mixed Mediterranean afternoons pass, so classes 3-4 become rare in
    # exactly the regime where pyroconvection occurs. If that is wrong, the fix is a moisture
    # criterion with a physical basis -- not a number chosen to match someone's map.
    rh_top_moist: float = 80.0

    # [UNSOURCED] Used by the scalar :func:`pyroconvection_score`, not by the ladder.
    residual_score: float = 35.0

    # [PYFLAM] Plausibility clips on the diagnosed depth, not classification thresholds.
    min_abl_m: float = _ABL_MIN_M
    max_abl_m: float = _ABL_MAX_M


DEFAULT_PYROCONV_THRESHOLDS = PyroconvThresholds()


def theta_kelvin(temp_c, pressure_hpa):
    """Potential temperature (K) from temperature (C) and pressure (hPa)."""
    return (np.asarray(temp_c, float) + 273.15) * (1000.0 / np.asarray(pressure_hpa, float)) ** 0.286


def theta_gradient(profile: "AtmosphericProfile", low_hpa: float, high_hpa: float) -> float:
    """Potential-temperature gradient dtheta/dz (K/m) across a pressure layer.

    Positive = stable. Heights from the hypsometric approximation. Used for both
    the mixed-layer stability (e.g. surface->850 hPa) and the upper-layer cap
    gamma-theta (e.g. 700->500 hPa) in the pyroconvection classification.
    """
    t_lo, _ = profile._at(low_hpa)
    t_hi, _ = profile._at(high_hpa)
    th_lo = theta_kelvin(t_lo, low_hpa)
    th_hi = theta_kelvin(t_hi, high_hpa)
    tbar = 0.5 * (t_lo + t_hi) + 273.15
    dz = 287.0 * tbar / _G * math.log(low_hpa / high_hpa)
    return float((th_hi - th_lo) / dz) if dz > 0 else 0.0


def bulk_richardson_abl_height(height_m, temperature_c, wind_u, wind_v, *,
                               pressure_hpa=None, surface_start_m: float = 400.0,
                               rib_crit: float = _RIB_CRITICAL) -> float:
    """ABL (mixing-layer) height (m) from the bulk Richardson number profile.

    Implements the Castellnou et al. (2022) procedure: compute the bulk Richardson
    number ``Rib(z) = (g/theta_s)(theta(z)-theta_s)(z-z_s) / (u^2+v^2)`` starting
    from ``surface_start_m`` and return the height where ``Rib`` first reaches
    ``rib_crit``, linearly interpolated.

    The defaults are the paper's own (sec.2.4.2): ``rib_crit`` **0.33** after Zhang et al.
    (2014), and a 400 m start. The 400 m is a *fire-sonde* correction -- Zhang recommends
    200 m, and the paper raises it because its sondes are launched into the plume indraft,
    which depresses the estimate (Charland & Clements 2013). This function is the in-plume
    counterpart, so it keeps 400 m; :func:`bulk_richardson_abl_grid` runs on ambient forecast
    columns with no indraft and keeps Zhang's 200 m. ``rib_crit`` was 0.25 here until
    2026-07-26, which matched neither the paper nor the gridded path.

    ``height_m`` (AGL), ``temperature_c`` and the wind components are 1-D arrays at
    the same levels; pass ``pressure_hpa`` to use exact potential temperature,
    otherwise theta is approximated by temperature. Falls back to the top height if
    the criterion is never reached.
    """
    z = np.asarray(height_m, float)
    order = np.argsort(z)
    z = z[order]
    if pressure_hpa is not None:
        th = theta_kelvin(np.asarray(temperature_c, float)[order],
                          np.asarray(pressure_hpa, float)[order])
    else:
        th = np.asarray(temperature_c, float)[order] + 273.15
    u = np.asarray(wind_u, float)[order]
    v = np.asarray(wind_v, float)[order]

    th_s = float(np.interp(surface_start_m, z, th))
    u_s = float(np.interp(surface_start_m, z, u))
    v_s = float(np.interp(surface_start_m, z, v))
    above = z > surface_start_m
    z2, th2, u2, v2 = z[above], th[above], u[above], v[above]
    if z2.size == 0:
        return float(z[-1])
    shear2 = np.maximum((u2 - u_s) ** 2 + (v2 - v_s) ** 2, 1e-6)
    rib = (_G / th_s) * (th2 - th_s) * (z2 - surface_start_m) / shear2

    prev_z, prev_r = surface_start_m, 0.0
    for zi, ri in zip(z2, rib):
        if ri >= rib_crit:
            if ri == prev_r:
                return float(zi)
            frac = (rib_crit - prev_r) / (ri - prev_r)
            return float(prev_z + frac * (zi - prev_z))
        prev_z, prev_r = zi, ri
    return float(z2[-1])


# --- moist thermodynamics for the profile diagnostics --------------------------

def saturation_vapour_pressure_pa(temp_k):
    """Saturation vapour pressure (Pa) over liquid water (Bolton 1980). Array-safe."""
    tc = np.asarray(temp_k, float) - 273.15
    return 611.2 * np.exp(17.67 * tc / (tc + 243.5))


def specific_humidity_from_rh(relative_humidity, temp_k, pressure_pa):
    """Specific humidity (kg/kg) from RH (%), temperature (K) and pressure (Pa)."""
    e = np.clip(np.asarray(relative_humidity, float) / 100.0, 0.0, 1.0) \
        * saturation_vapour_pressure_pa(temp_k)
    p = np.asarray(pressure_pa, float)
    return _EPSILON * e / np.maximum(p - (1.0 - _EPSILON) * e, 1.0)


def virtual_potential_temperature(theta_k, specific_humidity):
    """Virtual potential temperature (K) -- theta corrected for water-vapour buoyancy."""
    return np.asarray(theta_k, float) * (
        1.0 + 0.61 * np.asarray(specific_humidity, float))


def lcl_height_bolton_m(temp_k, dewpoint_k, pressure_pa):
    """LCL height above ground (m), exact form (Bolton 1980, eq. 21 + hypsometric).

    More faithful than the Espy rule in :func:`lcl_height_m`, which is kept for the
    surface-only diagnostics. Array-safe; returns 0 where the parcel is saturated.
    """
    t = np.asarray(temp_k, float)
    td = np.minimum(np.asarray(dewpoint_k, float), t)
    p = np.asarray(pressure_pa, float)
    # Temperature at the LCL (K).
    t_lcl = 1.0 / (1.0 / (td - 56.0) + np.log(t / td) / 800.0) + 56.0
    # Dry adiabat from the surface to the LCL, then hypsometric to a height.
    p_lcl = p * (t_lcl / t) ** (1.0 / _KAPPA_DRY)
    tbar = 0.5 * (t + t_lcl)
    z = _RD * tbar / _G * np.log(np.maximum(p, 1.0) / np.maximum(p_lcl, 1.0))
    return np.maximum(z, 0.0)


def bulk_richardson_abl_grid(height_agl_m, theta_v, wind_u, wind_v, *,
                             theta_v_surface, wind_u_surface, wind_v_surface,
                             fallback_m=None, rib_crit: float = _RIB_CRITICAL,
                             start_m: float = _RIB_START_M,
                             min_abl_m: float = _ABL_MIN_M,
                             max_abl_m: float = _ABL_MAX_M):
    """ABL depth (m AGL) over a whole grid from the bulk Richardson profile.

    The gridded counterpart of :func:`bulk_richardson_abl_height`: same
    ``Rib(z) = (g/theta_v_s)(theta_v(z)-theta_v_s) z / |U(z)-U_s|^2`` criterion, but
    vectorised over ``(nlev, ny, nx)`` stacks so a forecast grid is one array op
    instead of a per-cell Python loop.

    ``height_agl_m``, ``theta_v``, ``wind_u``, ``wind_v`` are level-major stacks
    (levels ascending in height); the ``*_surface`` arrays are the 2 m / 10 m
    reference fields. The ABL is the first height above ``start_m`` where Rib reaches
    ``rib_crit``, linearly interpolated between the bracketing levels. Cells that
    never cross, or that fall outside ``[min_abl_m, max_abl_m]``, take ``fallback_m``
    (the model's own boundary-layer height, if you have it) and are then clipped.

    The surface reference is the model's own 2 m / 10 m state, not a 400 m parcel:
    the paper's 400 m start applies to soundings *inside* the plume, where the fire's
    indraft contaminates the lowest levels. Here the profile is ambient, so the
    search still begins at ``start_m`` (200 m) to skip the surface layer.
    """
    z = np.asarray(height_agl_m, float)
    thv = np.asarray(theta_v, float)
    u = np.asarray(wind_u, float)
    v = np.asarray(wind_v, float)
    thv_s = np.asarray(theta_v_surface, float)[None, ...]
    u_s = np.asarray(wind_u_surface, float)[None, ...]
    v_s = np.asarray(wind_v_surface, float)[None, ...]

    shear2 = np.maximum((u - u_s) ** 2 + (v - v_s) ** 2, 0.5)
    rib = (_G / thv_s) * (thv - thv_s) * z / shear2

    shape = z.shape[1:]
    abl = np.full(shape, np.nan)
    prev_z = np.full(shape, start_m)
    prev_r = np.zeros(shape)
    for k in range(z.shape[0]):
        zk, rk = z[k], rib[k]
        # Levels that are underground (a 1000 hPa level over the Apennines) or
        # otherwise unusable arrive as nan and must not advance the bracket -- a nan
        # `prev_r` would poison the interpolation for every level above them.
        usable = np.isfinite(rk) & np.isfinite(zk) & (zk > start_m)
        hit = np.isnan(abl) & usable & (rk >= rib_crit)
        if hit.any():
            denom = np.where(rk != prev_r, rk - prev_r, np.inf)
            frac = np.clip((rib_crit - prev_r) / denom, 0.0, 1.0)
            abl = np.where(hit, prev_z + frac * (zk - prev_z), abl)
        prev_z = np.where(usable, zk, prev_z)
        prev_r = np.where(usable, rk, prev_r)

    if fallback_m is not None:
        fb = np.asarray(fallback_m, float)
        bad = ~np.isfinite(abl) | (abl < min_abl_m) | (abl > max_abl_m)
        abl = np.where(bad & np.isfinite(fb), fb, abl)

    # A column that never crosses Rib_crit is *well mixed through the profile* (a
    # superadiabatic afternoon surface layer gives Rib < 0 all the way up), not a
    # column we know nothing about. Returning nan there would silently drop exactly
    # the deeply-mixed cells that matter. Fall back to the top of the usable profile,
    # as the single-column bulk_richardson_abl_height does, and let the clip below cap
    # it. Only a column with no usable level at all stays nan.
    usable_z = np.where(np.isfinite(rib) & np.isfinite(z), z, np.nan)
    with np.errstate(invalid="ignore"):
        z_top = np.nanmax(np.where(np.isnan(usable_z), -np.inf, usable_z), axis=0)
    z_top = np.where(np.isfinite(z_top) & (z_top > -np.inf), z_top, np.nan)
    abl = np.where(~np.isfinite(abl) & np.isfinite(z_top), z_top, abl)

    return np.clip(abl, min_abl_m, max_abl_m)


def parcel_mixing_depth_grid(height_agl_m, theta, theta_surface, *,
                             excess_k: float = _PARCEL_EXCESS_K,
                             min_m: float = _ABL_MIN_M, max_m: float = _ABL_MAX_M):
    """Well-mixed-layer depth (m AGL) by the parcel method (Holzworth 1964), gridded.

    The first height at which ``theta`` exceeds the surface potential temperature by
    ``excess_k``, linearly interpolated -- i.e. the top of the layer a surface parcel
    mixes through, **below the entrainment jump**.

    This is *not* the same height as :func:`bulk_richardson_abl_grid`, and the
    difference is the point. The Rib height is where turbulence dies, which is at the
    top of the entrainment zone; the parcel height is the top of the well-mixed layer,
    at its base. Castellnou et al. (2022, sec.2.1.1) treat the mixed-layer gradient
    dtheta/dz and the entrainment jump Delta-theta as **separate** state variables of
    the mixed-layer slab model (Vila-Guerau de Arellano et al. 2015), and the
    classification ladder conditions on the former only.

    So: use the Rib height for the LCL/ABL and shear/ABL *ratios* (that is what the
    paper uses it for), and this height as the upper bound when measuring the
    mixed-layer dtheta/dz. Evaluating theta at the Rib top instead folds Delta-theta
    into the gradient and reports the capping inversion rather than the mixing --
    which, on a coarse profile, is enough to push a well-mixed summer CBL over the
    "stable" threshold and gate it out of the pyroCu branch entirely.

    ``height_agl_m`` and ``theta`` are level-major ``(nlev, ny, nx)`` stacks with
    levels ascending in height; ``theta_surface`` is the 2 m field. Unusable levels
    carry ``nan``. Columns that never cross take the top of the usable profile.
    """
    z = np.asarray(height_agl_m, float)
    th = np.asarray(theta, float)
    target = np.asarray(theta_surface, float) + excess_k

    out = np.full(target.shape, np.nan)
    for k in range(1, z.shape[0]):
        z0, z1, t0, t1 = z[k - 1], z[k], th[k - 1], th[k]
        good = np.isfinite(z0) & np.isfinite(z1) & np.isfinite(t0) & np.isfinite(t1)
        hit = np.isnan(out) & good & (t1 >= target)
        if hit.any():
            denom = np.where(t1 != t0, t1 - t0, np.inf)
            frac = np.clip((target - t0) / denom, 0.0, 1.0)
            out = np.where(hit, z0 + frac * (z1 - z0), out)

    usable = np.where(np.isfinite(z) & np.isfinite(th), z, -np.inf)
    z_top = np.max(usable, axis=0)
    out = np.where(~np.isfinite(out) & np.isfinite(z_top) & (z_top > -np.inf), z_top, out)
    return np.clip(out, min_m, max_m)


# --- fire-induced boundary layer (dry-pyrocloud encroachment) ------------------
#
# The moist diagnostics above (LCL, cap, shear, the Castellnou ladder) describe
# whether the *plume* reaches free moist convection. They do not describe the
# other route to fire-driven boundary-layer decoupling: the fire's *sensible*
# heat flux alone growing a fire-induced boundary layer (a "fireABL") that
# encroaches into the free troposphere with no condensation involved -- the
# dry-pyrocloud mechanism (Castellnou et al. 2022; Castellnou Ribau et al. 2024).
#
# This is the fire-forced analogue of ordinary daytime convective-boundary-layer
# growth by encroachment (Stull 1988, convective-mixed-layer chapter: "thermals
# rise until they hit the stable layer capping the ML"), with the fire's
# convective heat flux substituted for solar insolation. It ports stage 4 of the
# GRAF/WUR pipeline, with one physics fix: the fireABL top is found by
# intersecting the fire-heated parcel with the *actual* theta(z) profile (as in
# :func:`parcel_mixing_depth_grid`), not by extrapolating a short free-tropo
# linear fit -- the extrapolation runs the fireABL to unphysical heights (13 km
# from a 2 km fit) on the very hours it matters most.

# GRAF/WUR fire-plume scale height (m) in the convective-velocity scale. A named
# constant because it is a calibration target, not a first-principles value: it
# sets only the theta-excess magnitude (theta' ~ H^(-1/3)).
#
# CALIBRATION STATUS (magnitude is NOT calibrated; the ranking is):
#
#  * Linear-geometry fit vs the 90-fire sonde set: rescaling theta' by 1.08
#    (H ~ 56 m) removes the linear model's low bias. But that geometry is not the
#    one this module uses.
#  * Profile-intersection fit vs the SCQ sonde (the one fire with a full sounding
#    AND observed heights, 19 h): the profile fireABL OVER-predicts by ~3x
#    (median 6.1 km vs sonde 2.0 km), and the scale height that would fix it runs
#    off the top of the physical range (best H > 6 km, still +0.4 km biased). So
#    the H knob CANNOT calibrate the magnitude -- the theta' FORM is the problem:
#    F here is the fire-*front* flux ((Ib/2)/front_depth), a locally intense value,
#    used to force a mixed-layer-scale encroachment; a mixed-layer-averaged flux is
#    what that intersection should see. (Tellingly, GRAF's own published
#    ``expected_fireABL`` matches the sondes to within ~-370 m, i.e. their refined
#    model runs ~3x lower than the Demo stage-4 code ported here.)
#
# FIX (:func:`mixed_layer_fire_flux`): the over-prediction is a flux problem, not a
# scale-height one. Feed the encroachment a mixed-layer-averaged flux -- the fire's
# convective power over the ABL depth, not the ~10-20 m flaming-front depth. Confirmed
# across the 4 fires with ERA5 soundings + sonde fireABLs (36 h; zenodo_6433389):
# bias -15 m, MAE 400 m, r 0.91 (one global coefficient for the table's unit
# ambiguity), and it correctly captures the shallow-ABL strong-decoupling case
# (Torroella) the front flux misses by ~4x. So the recommended forcing for a real
# fire is ``mixed_layer_fire_flux(I, abl)``, not the front flux.
# Default scale height left at the reference 70 m; pass ``scale_height_m`` to override.
_FIRE_PLUME_SCALE_M = 70.0


def fire_parcel_theta_excess(heat_flux_w_m2, theta_mean_below_k, *,
                             scale_height_m: float = _FIRE_PLUME_SCALE_M):
    """Convective velocity and temperature scales for a surface heat flux (Deardorff).

    Mixed-layer similarity, array-safe:

        ``Q0     = F / (rho cp)``                    kinematic heat flux (K m/s)
        ``w*     = (3 g Q0 H / (2 theta_v))^(1/3)``  convective velocity (m/s)
        ``theta* = Q0 / w*``                         convective temperature scale (K)

    ``heat_flux_w_m2`` is the cell-averaged convective heat flux (W/m^2), e.g. from
    :func:`pyflam.pyroconvection.fire_heat_flux`; ``theta_mean_below_k`` is the mean
    potential temperature of the sub-ABL mixed layer (K). Returns ``(theta_star_k, w_star_ms)``.
    Zero flux gives zero excess (and w* -> 0 is handled, not a divide-by-zero).

    **Corrected 2026-07-26 -- the specific heat capacity was missing.** The previous form
    computed ``F / (rho w0)``, which carries units of J/kg, not kelvin, and was nonetheless
    added straight to a potential temperature in :func:`fire_induced_abl_grid`. It overstated
    the scale by ~88x: for F = 200 W/m^2 it returned 20.1 K where the dimensionally correct
    value is 0.23 K.

    **Do not use this as a fire plume's temperature excess.** theta* is the *turbulence*
    scale of a convectively mixed layer -- a few tenths of a kelvin for any realistic
    cell-averaged flux. A plume core is a coherent buoyant structure, not a turbulent
    fluctuation: the GRAF campaign measured in-plume anomalies of 0.1-13.1 K (Castellnou
    Ribau et al. 2025). The two differ by more than an order of magnitude, and no choice of
    ``heat_flux_w_m2`` reconciles them, because a cell-averaged flux is the fire's heat spread
    over ground that is mostly not burning. Prescribe the excess directly instead -- see
    ``theta_excess`` on :func:`fire_induced_abl_grid`, which is also what the campaign ingest
    (``scripts/ingest_inplume_sondes.py``) does deliberately, "so no fire intensity or flux
    parameterisation is needed".
    """
    f = np.asarray(heat_flux_w_m2, float)
    thv = np.asarray(theta_mean_below_k, float)
    pos = f > 0.0
    q0 = f / (_RHO_AIR * _CP_DRY)
    w0 = np.where(pos, (3.0 * _G * q0 * scale_height_m
                        / (2.0 * np.maximum(thv, 1.0))) ** (1.0 / 3.0), 0.0)
    excess = np.where(pos & (w0 > 0.0), q0 / np.maximum(w0, 1e-9), 0.0)
    return excess, w0


def mixed_layer_fire_flux(fireline_intensity_w_m, abl_m, *,
                          convective_fraction: float = 0.5):
    """Mixed-layer-averaged fire heat flux (W/m^2) for the fireABL encroachment.

    The fire's convective power per unit front (``convective_fraction * I``, I in
    W/m) spread over the **boundary-layer depth** -- the scale over which the heat
    is actually mixed -- rather than over the narrow flaming-front depth. Feed the
    result as ``heat_flux`` to :func:`fire_induced_abl_grid`.

    This fixes the fireABL *magnitude* bug (see the calibration note on
    ``_FIRE_PLUME_SCALE_M``): the reference stage-4 forcing uses the fire-front flux
    ``(I/2)/front_depth`` (front_depth ~ 10-20 m, a locally intense value), which
    drives the encroachment too high. Substituting the ABL depth for the front depth
    couples the forcing to the mixing scale.

    Validation. On the SCQ sonde alone (self-consistent units) the swap is
    parameter-free and collapses the bias from +3.6 km to +0.06 km (r 0.84 -> 0.94).
    Across the 4 fires with ERA5 soundings AND sonde-observed fireABLs (SCQ,
    Martorell, PoblaMassaluca, Torroella; 36 h; zenodo_6433389), with one global
    coefficient absorbing the aggregated table's ambiguous I units, the ML form
    gives bias -15 m, MAE 400 m, r 0.91. Its decisive advantage is the shallow-ABL,
    strong-decoupling regime (Torroella, ABL 315 m: ML +27 m vs the front flux's
    -444 m) -- exactly where the dry mechanism matters; on deep-ABL fires a rescaled
    front flux does nearly as well. In pyflam the intensity is in known units
    (W/m from :func:`pyflam.pyroconvection.fire_heat_flux`), so no coefficient is
    needed. ``convective_fraction`` defaults to 0.5 to match the reference ``I/2``.
    """
    i = np.asarray(fireline_intensity_w_m, float)
    return convective_fraction * i / np.maximum(np.asarray(abl_m, float), 1.0)


def fire_induced_abl_grid(height_agl_m, theta, *, theta_mean_below, blh,
                          heat_flux=None, theta_excess=None,
                          scale_height_m: float = _FIRE_PLUME_SCALE_M,
                          min_m: float = _ABL_MIN_M, max_m: float | None = None):
    """Fire-induced boundary-layer top (m AGL) by fire-forced encroachment, gridded.

    The height at which the fire-heated parcel
    ``theta_mean_below + theta'(heat_flux)`` reaches neutral buoyancy in the
    ambient sounding -- i.e. the first level *at or above* the ambient ABL where
    ambient ``theta`` catches the heated parcel, linearly interpolated. Same
    interpolation pattern as :func:`parcel_mixing_depth_grid`, but with the
    fire-forced target of :func:`fire_parcel_theta_excess` instead of a fixed
    surface excess, and the search anchored above ``blh`` (the mixed layer below
    is already accounted for by ``theta_mean_below``).

    Unlike the GRAF/WUR stage-4 original, the intersection is against the real
    ``theta(z)`` stack, so the result is bounded by the sounding: a column whose
    parcel is warmer than the entire profile returns the profile top (the fire
    penetrates the modelled domain) rather than an extrapolated height.

    ``height_agl_m`` and ``theta`` are level-major ``(nlev, ny, nx)`` stacks with
    levels ascending in height (K for ``theta``); ``theta_mean_below``,
    ``heat_flux`` (W/m^2) and ``blh`` (m) are ``(ny, nx)`` fields. Unusable levels
    carry ``nan``. Clipped to ``[min_m, max_m]`` (``max_m=None`` -> only the
    sounding bounds it) and never returned below ``blh``.

    Magnitude caveat: this reproduces the reference stage-4 mechanism and ranks
    columns well (r ~ 0.8 vs sondes) but its absolute heights run ~3x high vs
    sonde-observed fireABLs, and no scale height fixes it (see the calibration note
    on ``_FIRE_PLUME_SCALE_M``). Read the height and the decoupling ratio
    qualitatively until the forcing flux is recalibrated.

    Note the fireABL is a *height* diagnostic only; whether the decoupled fire
    also goes moist is a separate question -- compare a plume height against
    ``LCL + offset`` for that, where the in-plume ``offset`` is a configurable
    parameter (default 0 = ambient LCL). The +1 km cloud-base offset reported in
    the literature (Kablick 2018; Lareau & Clements 2016) is *not* supported by
    the GRAF prototype labels (a fit prefers ~ -0.5 km), so it is not baked in.
    """
    z = np.asarray(height_agl_m, float)
    th = np.asarray(theta, float)
    thm = np.asarray(theta_mean_below, float)
    blh = np.asarray(blh, float)
    if (heat_flux is None) == (theta_excess is None):
        raise ValueError("pass exactly one of heat_flux= or theta_excess=")
    if theta_excess is not None:
        # Prescribed plume excess (K), the preferred route: the flux -> theta' similarity
        # scaling returns a mixed-layer *turbulence* scale of a few tenths of a kelvin, an
        # order of magnitude below the 0.1-13.1 K anomalies the GRAF campaign measured inside
        # real plumes, and no cell-averaged flux closes that gap. The campaign ingest forces
        # the encroachment with the measured excess for the same reason.
        excess = np.asarray(theta_excess, float)
    else:
        excess, _ = fire_parcel_theta_excess(heat_flux, thm, scale_height_m=scale_height_m)
    target = thm + excess

    out = np.full(target.shape, np.nan)
    for k in range(1, z.shape[0]):
        z0, z1, t0, t1 = z[k - 1], z[k], th[k - 1], th[k]
        good = np.isfinite(z0) & np.isfinite(z1) & np.isfinite(t0) & np.isfinite(t1)
        # only bracket intervals reaching above the ambient ABL top
        hit = np.isnan(out) & good & (z1 >= blh) & (t1 >= target)
        if hit.any():
            denom = np.where(t1 != t0, t1 - t0, np.inf)
            frac = np.clip((target - t0) / denom, 0.0, 1.0)
            out = np.where(hit, z0 + frac * (z1 - z0), out)

    # Parcel warmer than the whole profile -> fire penetrates the domain top.
    usable = np.where(np.isfinite(z) & np.isfinite(th), z, -np.inf)
    z_top = np.max(usable, axis=0)
    out = np.where(~np.isfinite(out) & np.isfinite(z_top) & (z_top > -np.inf), z_top, out)
    hi = z_top if max_m is None else np.minimum(z_top, max_m)
    out = np.clip(out, min_m, hi)
    # The fireABL cannot sit below the ambient ABL it grows out of. Applied last so
    # the sounding clip cannot push it back under the ABL on a thin/degenerate column
    # (profile top below the ABL), where ``hi`` would otherwise win.
    return np.maximum(out, blh)


def entrainment_jump_grid(height_agl_m, theta, abl_m, *,
                          ml_lo_frac: float = 0.10, ml_hi_frac: float = 0.80,
                          ez_frac: float = 0.15, ez_min_m: float = 100.0,
                          samples: int = 5):
    """Entrainment-zone potential-temperature jump ``delta-theta`` (K), gridded.

    The step in ``theta`` across the entrainment zone, between the well-mixed layer and the
    free atmosphere above the ABL top -- the ``d(theta)`` of mixed-layer (slab) theory
    (Vila-Guerau de Arellano et al. 2015; Stull 1988), drawn as the ``theta jump`` in
    Castellnou et al. (2022) Fig. 1c.

    This is the quantity the source method names as the control on penetration:

        *"The jumps at the entrainment zone between ABL and free atmosphere (dq, d-theta)
        and lapse rates on the free atmosphere assess the ability of a parcel to penetrate
        above ABL and achieve free convection."* -- Castellnou et al. (2022) sec.2.1.2

    pyflam previously carried only the free-atmosphere lapse rate ``gamma_theta`` and no jump
    at all, so nothing in the ladder expressed how hard it is to *get out of* the mixed layer
    -- only how the plume fares once out. On a well-mixed summer afternoon that omission lets
    essentially every column qualify, which is the mechanism behind the over-extent documented
    in ``docs/graf_vs_pyflam_2026-07-26.md``.

    In a real (non-slab) profile the jump is smeared over an entrainment zone of order
    0.1-0.2 of the ABL depth rather than being a discontinuity, so it is measured as
    ``theta`` at ``abl + max(ez_frac*abl, ez_min_m)`` minus the mixed-layer mean over
    ``[ml_lo_frac, ml_hi_frac] * abl``. Both bounds are deliberately inside their layers:
    starting the mean above the superadiabatic surface layer and ending it below the
    entrainment zone keeps the two terms from sampling each other.

    Returns a 2-D array (K), ``nan`` where the column cannot support the estimate. Values are
    clipped at 0 from below: a negative jump is an unstable top, which mixed-layer theory does
    not admit and which here signals a degenerate column rather than a physical state.
    """
    z = np.asarray(height_agl_m, float)
    th = np.asarray(theta, float)
    abl = np.asarray(abl_m, float)

    lo = np.maximum(ml_lo_frac * abl, 20.0)
    hi = np.maximum(ml_hi_frac * abl, lo + 50.0)
    acc = np.zeros(abl.shape)
    for i in range(samples):
        acc = acc + _interp_profile_at(z, th, lo + (hi - lo) * (i / (samples - 1.0)))
    th_ml = acc / samples

    z_top = abl + np.maximum(ez_frac * abl, ez_min_m)
    th_top = _interp_profile_at(z, th, z_top)

    jump = th_top - th_ml
    ok = np.isfinite(jump) & np.isfinite(abl) & (abl > 0)
    return np.where(ok, np.maximum(jump, 0.0), np.nan)


def _interp_profile_at(z, x, ztarget):
    """Per-column linear interpolation of ``x(z)`` to a 2-D target height (nan-aware).

    Level-major ``(nlev, ny, nx)`` stacks with ``z`` ascending on axis 0; targets outside the
    usable span clamp to the nearest usable level. Shared by the entrainment-jump and
    residual-layer diagnostics so they cannot drift apart in how they read a profile.
    """
    z = np.asarray(z, float)
    x = np.asarray(x, float)
    out = np.full(np.shape(ztarget), np.nan)
    for k in range(z.shape[0] - 1):
        z0, z1, x0, x1 = z[k], z[k + 1], x[k], x[k + 1]
        span = z1 - z0
        good = (np.isfinite(z0) & np.isfinite(z1) & np.isfinite(x0) & np.isfinite(x1)
                & (np.abs(span) > 1e-6))
        f = np.divide(ztarget - z0, span, out=np.zeros_like(out), where=good)
        seg = np.isnan(out) & good & (ztarget >= z0) & (ztarget <= z1)
        out = np.where(seg, x0 + np.clip(f, 0.0, 1.0) * (x1 - x0), out)

    ok = np.isfinite(z) & np.isfinite(x)
    any_ok = ok.any(axis=0)
    i_lo = np.argmax(ok, axis=0)
    i_hi = (ok.shape[0] - 1) - np.argmax(ok[::-1], axis=0)
    take = lambda a, i: np.take_along_axis(a, i[None, ...], axis=0)[0]
    z_lo, x_lo = take(z, i_lo), take(x, i_lo)
    z_hi, x_hi = take(z, i_hi), take(x, i_hi)
    out = np.where(np.isnan(out) & any_ok & (ztarget <= z_lo), x_lo, out)
    out = np.where(np.isnan(out) & any_ok & (ztarget >= z_hi), x_hi, out)
    return out


def fire_cape_grid(height_agl_m, theta_v, *, theta_excess, top_m=None):
    """FireCAPE (J/kg) -- convective available potential energy of a fire-heated parcel.

    Potter (2005), as used by Castellnou et al. (2022) Eq. 2: the environmental CAPE
    recomputed for a surface parcel carrying the fire's potential-temperature perturbation,

        ``FireCAPE = g * integral[ (theta_parcel - theta_v(z)) / theta_v(z) ] dz``

    over the layer where the heated parcel is positively buoyant. The parcel conserves
    ``theta_v(surface) + theta_excess`` (dry ascent, no condensation -- the moist contribution
    is deliberately excluded, matching the dry framing of the fireABL diagnostic), so this is
    the *dry* fire CAPE, and only positive buoyancy is accumulated, as for ordinary CAPE.

    ``theta_excess`` is the fire's perturbation (K), e.g. from
    :func:`fire_parcel_theta_excess`. Returns a 2-D array; 0 where the parcel is nowhere
    buoyant, ``nan`` where the column is unusable.

    This closes the second of the two variables the source method computes and pyflam did not
    (the other being the entrainment jump). Note it is *not* surface CAPE, which this product
    deliberately omits -- the fire perturbation is what makes it meaningful for a plume.
    """
    z = np.asarray(height_agl_m, float)
    thv = np.asarray(theta_v, float)
    exc = np.asarray(theta_excess, float)

    ok = np.isfinite(z) & np.isfinite(thv)
    any_ok = ok.any(axis=0)
    i_lo = np.argmax(ok, axis=0)
    thv_sfc = np.take_along_axis(thv, i_lo[None, ...], axis=0)[0]
    parcel = thv_sfc + exc

    cape = np.zeros(np.shape(parcel))
    for k in range(z.shape[0] - 1):
        z0, z1, t0, t1 = z[k], z[k + 1], thv[k], thv[k + 1]
        dz = z1 - z0
        good = (np.isfinite(z0) & np.isfinite(z1) & np.isfinite(t0) & np.isfinite(t1)
                & (dz > 0))
        if top_m is not None:
            good = good & (z1 <= top_m)
        b0 = (parcel - t0) / np.where(good, np.maximum(t0, 1.0), 1.0)
        b1 = (parcel - t1) / np.where(good, np.maximum(t1, 1.0), 1.0)
        # trapezoid over the positive-buoyancy part of the layer only
        both = good & (b0 > 0) & (b1 > 0)
        cape = np.where(both, cape + _G * 0.5 * (b0 + b1) * dz, cape)
        part = good & ((b0 > 0) ^ (b1 > 0))
        frac = np.abs(b0) / np.maximum(np.abs(b0) + np.abs(b1), 1e-12)
        cape = np.where(part, cape + _G * 0.5 * np.maximum(np.maximum(b0, b1), 0.0)
                        * dz * np.where(b0 > 0, frac, 1.0 - frac), cape)
    return np.where(any_ok, cape, np.nan)


def plume_entrainment_fraction(lcl_m, z_fc_m, *, beta: float = 0.6):
    """Mass fraction of environment entrained between cloud base and free convection.

    A bent-over Briggs plume has radius ``R = beta * z`` and mass flux proportional to
    ``R^2 U``, so rising from ``z1`` to ``z2`` its mass grows by ``(z2/z1)^2`` and the
    entrained share of the final mass is ``1 - (z1/z2)^2``. (``beta`` cancels; it is kept in
    the signature because the *form* is Briggs's and the reader should see which plume model
    is being assumed.)

    The interval is not a free choice. It runs from the **LCL**, below which there is no
    condensate to evaporate, to the **free-convection height** ``z_fc``
    (:func:`pyrocb_firepower_threshold_grid`), above which the parcel is buoyant on its own
    and no longer depends on surviving entrainment. Both ends are already diagnosed, so this
    closes :func:`critical_rh_for_cloud_persistence` with no additional parameter.

    Assumes a bent-over plume. For an upright plume in weak wind the geometry differs, and
    this is the least-settled part of the moisture criterion.
    """
    z1 = np.asarray(lcl_m, float)
    z2 = np.asarray(z_fc_m, float)
    with np.errstate(invalid="ignore", divide="ignore"):
        chi = 1.0 - (z1 / np.where(z2 > 0, z2, np.nan)) ** 2
    return np.clip(chi, 0.05, 1.0)


def critical_rh_for_cloud_persistence(theta_excess_k, temp_k, pressure_pa, *,
                                      mixing_fraction: float = 1.0):
    """Minimum ambient RH (%) for a plume cloud to survive entrainment -- derived, not assumed.

    A pyrocloud that has condensed keeps its buoyancy only if the liquid water it carries is
    not evaporated away by the sub-saturated air it entrains. Saturating a mass fraction
    ``mixing_fraction`` of environmental air costs the parcel

        ``dq   = q_s(T, p) * (1 - RH/100)``      moisture deficit of the entrained air
        ``dth  = (Lv/cp) * mixing_fraction * dq``  latent cooling paid to close it

    and the cloud persists while the plume's own buoyancy excess covers that cost. Setting
    ``dth = theta_excess`` and solving for RH gives the criterion returned here:

        ``RH_crit = 100 * (1 - theta_excess * cp / (Lv * mixing_fraction * q_s(T, p)))``

    ``mixing_fraction = 1`` is equal-mass mixing, the conservative limit and the usual choice
    in the cloud-top entrainment instability literature (Randall 1980; Deardorff 1980).

    **Why this replaces a fixed threshold.** The ladder inherited a constant 80 % RH at the ABL
    top for classes 3-4, with no citation in Castellnou et al. (2022); it had at one point been
    lowered to 60 % to match a third party's forecast. Neither is a physical statement. Tory &
    Kepert (2021) point at this derivation themselves -- their footnote 10 notes the buoyancy
    buffer "would ideally vary with evaporation potential of entrained air (e.g., the difference
    in specific humidity between the plume parcel and the environment)", after Peterson et al.
    (2017) on mid-tropospheric humidity as a pyroCb control.

    The criterion is **not** a constant: it falls with height as ``q_s`` falls, and falls as the
    fire gets stronger. Evaluated against the plume excesses the GRAF campaign measured
    (0.1-13.1 K, median 4.0), the inherited constants correspond to particular fire strengths
    rather than to atmospheric physics -- 80 % is roughly a 6-8 K plume at 1-2 km, and 60 % is
    a ~13 K plume, at the very top of the observed range. That is why 60 % over-produced the
    deep class.

    Returns RH in percent, clipped to [0, 100]. A cloud in air already at 100 % RH loses nothing
    to evaporation, so the criterion is always satisfiable in a saturated environment.
    """
    thx = np.asarray(theta_excess_k, float)
    t = np.asarray(temp_k, float)
    p = np.asarray(pressure_pa, float)
    es = saturation_vapour_pressure_pa(t)
    qs = _EPSILON * es / np.maximum(p - (1.0 - _EPSILON) * es, 1.0)
    denom = _LV * max(mixing_fraction, 1e-6) * np.maximum(qs, 1e-9)
    return np.clip(100.0 * (1.0 - thx * _CP_DRY / denom), 0.0, 100.0)


def max_rh_abl_grid(height_agl_m, relative_humidity, *, zmin: float = 150.0,
                    zmax: float = 4000.0, min_drop_pct: float = 2.0, smooth: int = 1):
    """ABL height (m AGL) as the height of maximum relative humidity -- the GRAF criterion.

    Castellnou Ribau et al. (2025, AMT 18, 7805-7831) sec. 2.6 define the boundary-layer top
    this way rather than by a Richardson crossing:

        *"The height of the maximum RH value is used as a criterion to estimate the height of
        the atmospheric boundary layer. This criterion is based on the observation that
        specific humidity tends to be well mixed in the convective boundary layer... temperature
        decreases with height, leading to an increase in relative humidity with altitude,
        reaching a peak at the inversion level. Above this inversion, the air becomes drier and
        warmer, resulting in a decrease in RH."*

    It is a *moisture* criterion where :func:`bulk_richardson_abl_grid` is a *dynamic* one, and
    on the campaign's own sondes the two disagree by 0.5x to 4.7x -- enough to move a column
    across class boundaries. Provided so the two definitions can be compared rather than
    conflated; the classifier still uses the Rib depth by default.

    The peak is sought between ``zmin`` (above the surface layer, where a shallow morning
    moisture maximum would otherwise win) and ``zmax``. A column qualifies only if RH falls by
    at least ``min_drop_pct`` somewhere above the peak -- without that guard a profile whose RH
    rises monotonically to the top of the data returns its top level, which is not an inversion.
    ``smooth`` is a running-mean half-width in levels, for noisy sonde data; model profiles need
    none. Returns ``nan`` where no qualifying peak exists.

    Note the source applies this *visually* on a plotted profile, with a human rejecting
    spurious maxima. This is the automated analogue and will not always agree with a hand
    reading.
    """
    z = np.asarray(height_agl_m, float)
    rh = np.asarray(relative_humidity, float)
    if smooth > 1:
        k = np.ones(smooth) / smooth
        rh = np.apply_along_axis(lambda a: np.convolve(a, k, mode="same"), 0, rh)

    win = (z >= zmin) & (z <= zmax) & np.isfinite(rh) & np.isfinite(z)
    cand = np.where(win, rh, -np.inf)
    kmax = np.argmax(cand, axis=0)
    peak_rh = np.take_along_axis(cand, kmax[None, ...], axis=0)[0]
    peak_z = np.take_along_axis(np.where(np.isfinite(z), z, np.nan), kmax[None, ...], axis=0)[0]

    # the peak must be followed by a real decrease aloft, else it is not an inversion
    above = z > peak_z[None, ...]
    rh_above = np.where(above & np.isfinite(rh), rh, np.inf)
    drop = peak_rh - np.min(rh_above, axis=0)
    ok = win.any(axis=0) & np.isfinite(peak_z) & (drop >= min_drop_pct)
    return np.where(ok, peak_z, np.nan)


def inversion_height_grid(height_agl_m, theta, *, zmin: float = 150.0,
                          zmax: float = 4000.0, min_grad: float = 2.0e-3):
    """Height (m AGL) of the strongest potential-temperature gradient -- the capping inversion.

    A structural marker independent of both the Richardson and maximum-RH criteria, used to
    corroborate them (:func:`corroborated_abl_grid`). Returns ``nan`` where no layer in the
    window reaches ``min_grad``, i.e. where there is no inversion worth calling a cap.
    """
    z = np.asarray(height_agl_m, float)
    th = np.asarray(theta, float)
    dz = np.diff(z, axis=0)
    grad = np.divide(np.diff(th, axis=0), dz, out=np.full_like(dz, np.nan),
                     where=np.abs(dz) > 1e-6)
    mid = 0.5 * (z[:-1] + z[1:])
    win = (mid >= zmin) & (mid <= zmax) & np.isfinite(grad)
    cand = np.where(win, grad, -np.inf)
    k = np.argmax(cand, axis=0)
    best = np.take_along_axis(cand, k[None, ...], axis=0)[0]
    hgt = np.take_along_axis(np.where(np.isfinite(mid), mid, np.nan), k[None, ...], axis=0)[0]
    return np.where(win.any(axis=0) & (best >= min_grad), hgt, np.nan)


def corroborated_abl_grid(height_agl_m, theta, relative_humidity, rib_abl_m, *,
                          tol_frac: float = 0.35, **maxrh_kwargs):
    """ABL depth using the maximum-RH criterion **where an inversion corroborates it**.

    Castellnou Ribau et al. (2025) sec. 2.6 identify the boundary-layer top as the height of
    maximum relative humidity, and state that they *"supplement this visually-based approach
    with numerical calculations using the bulk Richardson number"* -- the two are meant to be
    used together, not chosen between. That matters because the criterion is applied by a human
    reading a plotted profile, who discards spurious maxima; an automated maximum does not, and
    picks near-surface moisture peaks in a minority of columns.

    Here the maximum-RH height is accepted only when it sits within ``tol_frac`` of the
    strongest theta gradient in the column (:func:`inversion_height_grid`) -- the inversion it
    is supposed to be detecting. Where the two disagree, or where no inversion is resolvable,
    the bulk-Richardson depth is used instead. Returns ``(abl_m, from_max_rh)`` so a caller can
    report what fraction of the domain each criterion supplied.

    The tolerance is deliberately loose: the two markers measure different things (a moisture
    maximum sits at the inversion *base*, the gradient maximum inside it), so this tests
    agreement of structure, not of value.
    """
    z = np.asarray(height_agl_m, float)
    rib = np.asarray(rib_abl_m, float)
    z_rh = max_rh_abl_grid(z, relative_humidity, **maxrh_kwargs)
    z_inv = inversion_height_grid(z, theta)
    with np.errstate(invalid="ignore"):
        agree = (np.isfinite(z_rh) & np.isfinite(z_inv)
                 & (np.abs(z_rh - z_inv) <= tol_frac * np.maximum(z_inv, 1.0)))
    return np.where(agree, z_rh, rib), agree


def residual_layer_grid(height_agl_m, theta, *, search_max_m: float = 3000.0,
                        excess_k: float = _PARCEL_EXCESS_K, min_top_m: float = 200.0):
    """Residual-layer top (m AGL) -- the well-mixed depth left over after the CBL decays.

    After sunset the convective mixed layer collapses into a shallow stable layer, but the
    air above it keeps the day's near-neutral profile: the *residual layer*
    (Stull 1988, ch. 12). A surface-referenced parcel method finds only the shallow stable
    layer and reports a depth of tens of metres; a least-squares mixed-layer gradient fitted
    inside that depth has no levels to work with and returns ``nan``. That is what leaves the
    evening hours unclassifiable even after the ABL floor is removed -- the ladder's stability
    diagnostic is ill-posed in a regime the fire literature cares about, since the campaign of
    Castellnou Ribau et al. (2025) sampled late-afternoon and evening fires.

    Here the parcel is referenced not to the surface but to the **coldest level in the lower
    column** -- the top of the nocturnal stable layer -- and the residual top is the first
    height above it where ``theta`` exceeds that reference by ``excess_k``. By day the theta
    minimum sits at the surface and this reduces to the ordinary parcel mixing depth, so the
    same call is valid around the clock; only in the evening do the two separate.

    Returns a 2-D array (m AGL); ``nan`` where no usable column exists.
    """
    z = np.asarray(height_agl_m, float)
    th = np.asarray(theta, float)
    usable = np.isfinite(z) & np.isfinite(th) & (z <= search_max_m)

    big = np.where(usable, th, np.inf)
    k_ref = np.argmin(big, axis=0)
    take = lambda a: np.take_along_axis(a, k_ref[None, ...], axis=0)[0]
    z_ref, th_ref = take(z), take(th)
    target = th_ref + excess_k

    out = np.full(z_ref.shape, np.nan)
    for k in range(z.shape[0] - 1):
        z0, z1, t0, t1 = z[k], z[k + 1], th[k], th[k + 1]
        good = (np.isfinite(z0) & np.isfinite(z1) & np.isfinite(t0) & np.isfinite(t1)
                & (z0 >= z_ref))
        hit = np.isnan(out) & good & (t1 >= target)
        denom = np.where(t1 != t0, t1 - t0, np.inf)
        frac = np.clip((target - t0) / denom, 0.0, 1.0)
        out = np.where(hit, z0 + frac * (z1 - z0), out)

    z_top = np.max(np.where(usable, z, -np.inf), axis=0)
    out = np.where(np.isnan(out) & (z_top > -np.inf), z_top, out)
    return np.where(usable.any(axis=0), np.maximum(out, min_top_m), np.nan)


_LV = 2.501e6                  # J kg-1, latent heat of vaporisation at 0 C
# Tory & Kepert (2021) PFT constants. eq 31 is their operational reduction of eq 25 with
# rho_0 = 0.755 kg/m3 folded in: PFT[GW] = 0.3 * z_fc[km]^2 * U_ML[m/s] * dtheta_fc[K].
_PFT_C = 0.3
_PFT_CLOUD_TOP_C = -20.0       # conservative electrification level (their assumption 5)
_PFT_BUOY_BUFFER_K = 0.5       # dtheta_b, for evaporative loss to entrained dry air
_PFT_HEAT_TO_MOISTURE = 15.0   # K per g/kg, the fire's heat:moisture ratio (assumption 4)


def _moist_lapse_rate(temp_k, pressure_pa):
    """Saturated (pseudo)adiabatic lapse rate, K/m. Array-safe."""
    t = np.asarray(temp_k, float)
    p = np.asarray(pressure_pa, float)
    es = saturation_vapour_pressure_pa(t)
    qs = _EPSILON * es / np.maximum(p - (1.0 - _EPSILON) * es, 1.0)
    num = 1.0 + _LV * qs / (_RD * t)
    den = 1.0 + _LV ** 2 * qs * _EPSILON / (_CP_DRY * _RD * t ** 2)
    return _G / _CP_DRY * num / den


def pyrocb_firepower_threshold_grid(height_agl_m, temperature_k, pressure_pa,
                                    spec_humidity, wind_u, wind_v, *,
                                    surface_pressure_pa,
                                    buoyancy_buffer_k: float = _PFT_BUOY_BUFFER_K,
                                    cloud_top_temp_c: float = _PFT_CLOUD_TOP_C,
                                    heat_to_moisture: float = _PFT_HEAT_TO_MOISTURE,
                                    beta_max: float = 0.12, beta_steps: int = 60):
    """PyroCb Firepower Threshold (GW), gridded -- Tory & Kepert (2021) eq 31.

    The minimum **total** firepower that, in this atmosphere, lifts a fire plume to where
    free moist convection can carry it to the electrification level. Smaller values favour
    pyroCb. Returns a dict with ``pft_gw`` and its three drivers ``z_fc_m``,
    ``delta_theta_fc_k`` and ``u_ml_ms``, plus the mixed-layer state used.

        PFT[GW] = 0.3 * z_fc[km]^2 * U_ML[m/s] * dtheta_fc[K]

    The procedure follows their section 4b-4c rather than the earlier inversion in
    :func:`pyrocb_firepower_threshold`, which targeted the **LCL** and used a bulk
    Brunt-Vaisala stability. Those are not the paper's quantities: ``z_fc`` is the
    *free-convection* height, 3.5-4.8 km in their worked cases against a ~1.5 km LCL, and
    ``dtheta_fc`` is a parcel-vs-mixed-layer potential-temperature difference, not a
    stability. Scoring the LCL understates PFT by orders of magnitude.

    Steps, per column:

    1. **Mixed layer.** ``theta_ML``/``q_ML`` are height-weighted means up to the ML-LCL
       (weighted because entrained mass flux grows linearly with height, their eqs 15-16),
       found by iterating the ML depth until the ML-LCL sits inside it.
    2. **Saturation-point curve** (Tory et al. 2018). For a fire buoyancy increment
       ``beta``, the plume parcel carries ``theta_SP = (1 + beta) theta_ML`` and
       ``q_SP = q_ML + beta theta_ML / heat_to_moisture`` (their eqs 29-30, with the
       assumed 15 K per 1 g/kg fire heat-to-moisture ratio). Its saturation point is the
       height at which that parcel condenses.
    3. **Free-convection height.** The smallest ``beta`` whose moist adiabat from the
       saturation point stays warmer than the environment -- by ``buoyancy_buffer_k``, the
       allowance for evaporative cooling by entrained dry air -- all the way to the
       ``cloud_top_temp_c`` level. ``z_fc`` is that saturation-point height and
       ``dtheta_fc = beta theta_ML``.
    4. **Wind.** ``U_ML`` is the magnitude of the *vector* mean wind between the surface and
       ``z_fc`` (their step 5) -- not the 10 m wind.

    ``nan`` where no ``beta <= beta_max`` produces a buoyant cloud: the column cannot make a
    pyroCb at any firepower, which is a stronger statement than a large threshold.

    Scale check against their published cases: Black Saturday 1000 LST (z_fc 4.8 km,
    U_ML 20 m/s, dtheta_fc 9 K) gives 1244 GW against their stated 1240. Real events run
    ~100 GW (Chisholm afternoon) to ~1240 GW (Black Saturday morning), with Sir Ivan at
    ~300 GW described as near the upper limit for most wildfires. Compare against a fire's
    **total** power -- Byram intensity times head-fire length, as they do (100 MW/m over a
    5 km front = 500 GW) -- not against an intensity.
    """
    z = np.asarray(height_agl_m, float)
    T = np.asarray(temperature_k, float)
    p = np.asarray(pressure_pa, float)
    q = np.asarray(spec_humidity, float)
    u, v = np.asarray(wind_u, float), np.asarray(wind_v, float)
    ps = np.asarray(surface_pressure_pa, float)
    shape = ps.shape

    # --- 1. mixed layer: iterate depth until the ML-LCL sits inside it.
    #
    # The window is measured from each column's own lowest usable level, not from 0 m. A
    # gridded profile starts within ~10 m of the ground and the distinction is immaterial,
    # but a radiosonde ascent can begin well above it -- the GRAF campaign's Guissona sonde
    # starts at 840 m AGL -- and an absolute `z <= 500` window then selects *no* levels,
    # making theta_ML nan and every subsequent beta fail. That produced a nan PFT which read
    # as "no firepower suffices" when it actually meant "the mixed layer was never sampled".
    theta = theta_kelvin(T - 273.15, p / 100.0)
    usable = np.isfinite(theta) & np.isfinite(q) & np.isfinite(z)
    z_base = np.min(np.where(usable, z, np.inf), axis=0)
    z_base = np.where(np.isfinite(z_base), z_base, 0.0)
    depth = np.full(shape, 500.0)
    for _ in range(6):
        w = usable & (z <= (z_base + depth)[None, ...])
        # Empty window (a coarse profile whose first level already clears the ML-LCL): fall
        # back to the lowest usable level so the column still yields a mixed-layer state.
        if not w.any():
            w = usable & (z <= (z_base + 1e-6)[None, ...])
        else:
            empty = ~w.any(axis=0)
            if empty.any():
                w = w | (usable & empty[None, ...] & (z <= (z_base + 1e-6)[None, ...]))
        # linear-in-height weighting: entrained mass flux grows with height (eqs 15-16)
        wt = np.where(w, np.maximum(z, 1.0), 0.0)
        tot = np.maximum(wt.sum(0), 1e-9)
        th_ml = (wt * np.where(w, theta, 0.0)).sum(0) / tot
        q_ml = (wt * np.where(w, q, 0.0)).sum(0) / tot
        t_sfc = th_ml * (ps / 100000.0) ** _KAPPA_DRY
        e = q_ml * ps / (_EPSILON + (1.0 - _EPSILON) * q_ml)
        td = 243.5 * np.log(np.maximum(e, 1.0) / 611.2) / (17.67 - np.log(np.maximum(e, 1.0) / 611.2))
        depth = np.clip(lcl_height_bolton_m(t_sfc, td + 273.15, ps), 100.0, 5000.0)

    # --- 2-3. walk the saturation-point curve for the critical beta
    z_fc = np.full(shape, np.nan)
    dth_fc = np.full(shape, np.nan)
    for beta in np.linspace(0.0, beta_max, beta_steps)[1:]:
        todo = ~np.isfinite(z_fc)
        if not todo.any():
            break
        th_p = (1.0 + beta) * th_ml
        q_p = q_ml + beta * th_ml / heat_to_moisture * 1e-3      # g/kg -> kg/kg
        t_p = th_p * (ps / 100000.0) ** _KAPPA_DRY
        e = q_p * ps / (_EPSILON + (1.0 - _EPSILON) * q_p)
        le = np.log(np.maximum(e, 1.0) / 611.2)
        td = 243.5 * le / (17.67 - le)
        z_sp = lcl_height_bolton_m(t_p, td + 273.15, ps)
        t_sp = t_p - _G / _CP_DRY * z_sp                          # dry ascent to the SP
        p_sp = ps * np.maximum(t_sp / np.maximum(t_p, 1.0), 1e-6) ** (1.0 / _KAPPA_DRY)

        # moist ascent from the SP; must beat the environment by the buffer up to the
        # cloud-top level, else this beta is not enough
        tp, pp, zp = t_sp.copy(), p_sp.copy(), z_sp.copy()
        ok = np.ones(shape, bool)
        reached = np.zeros(shape, bool)
        for k in range(z.shape[0]):
            zk = z[k]
            step = zk - zp
            adv = step > 0
            if not adv.any():
                continue
            tp = np.where(adv, tp - _moist_lapse_rate(tp, pp) * np.maximum(step, 0.0), tp)
            pp = np.where(adv, pp * np.exp(-_G * np.maximum(step, 0.0) / (_RD * np.maximum(tp, 1.0))), pp)
            zp = np.where(adv, zk, zp)
            env = T[k]
            above_sp = adv & np.isfinite(env)
            ok &= ~(above_sp & ~reached & (tp < env + buoyancy_buffer_k))
            reached |= above_sp & (env <= cloud_top_temp_c + 273.15)
        good = todo & ok & reached
        z_fc = np.where(good, z_sp, z_fc)
        dth_fc = np.where(good, beta * th_ml, dth_fc)

    # --- 4. vector-mean wind between the surface and z_fc
    inlay = (z <= np.where(np.isfinite(z_fc), z_fc, 0.0)[None, ...]) & np.isfinite(u) & np.isfinite(v)
    n = np.maximum(inlay.sum(0), 1)
    u_ml = np.hypot((np.where(inlay, u, 0.0)).sum(0) / n, (np.where(inlay, v, 0.0)).sum(0) / n)

    pft = _PFT_C * (z_fc / 1000.0) ** 2 * u_ml * dth_fc
    return dict(pft_gw=pft, z_fc_m=z_fc, delta_theta_fc_k=dth_fc, u_ml_ms=u_ml,
                theta_ml_k=th_ml, q_ml=q_ml, ml_depth_m=depth)


def shear_height_grid(height_agl_m, wind_u, wind_v, *, zmin: float = 200.0,
                      zmax: float = 6000.0):
    """Height (m AGL) of maximum vector wind shear over a grid, from consecutive levels.

    The gridded, native-model-level analogue of :func:`shear_height_window`. On a resolved
    model-level column (dz of tens to a few hundred metres) the consecutive-level shear
    ``|dU/dz|`` already resolves the shear maximum, so no 50 m re-interpolation is needed and
    the whole grid is one vectorised op instead of a per-cell loop. This is what makes the
    fifth (shear) diagnostic affordable on the ICON-EU model-level path, which the coarse
    5-pressure-level path cannot supply (hence :func:`shear_height_none` there).

    ``height_agl_m``/``wind_u``/``wind_v`` are ``(nlev, ny, nx)`` stacks with levels ascending
    in height. Returns the 2-D height of the strongest shear between adjacent levels whose
    midpoint lies in ``[zmin, zmax]``; ``nan`` where no level pair qualifies.
    """
    z = np.asarray(height_agl_m, float)
    u = np.asarray(wind_u, float)
    v = np.asarray(wind_v, float)
    dz = z[1:] - z[:-1]
    zmid = 0.5 * (z[1:] + z[:-1])
    with np.errstate(invalid="ignore", divide="ignore"):
        shear = np.hypot(u[1:] - u[:-1], v[1:] - v[:-1]) / np.where(dz > 1.0, dz, np.nan)
    ok = (zmid >= zmin) & (zmid <= zmax) & np.isfinite(shear)
    masked = np.where(ok, shear, -np.inf)
    k = np.argmax(masked, axis=0)
    zsel = np.take_along_axis(zmid, k[None, ...], axis=0)[0]
    return np.where(ok.any(axis=0), zsel, np.nan)


def shear_distance_grid(shear_height_m, abl_m, lcl_m):
    """Gridded ``min(|z_shear - ABL|, |z_shear - LCL|) / ABL`` (the ladder's shear test)."""
    zs = np.asarray(shear_height_m, float)
    abl = np.asarray(abl_m, float)
    lcl = np.asarray(lcl_m, float)
    with np.errstate(invalid="ignore", divide="ignore"):
        sd = np.minimum(np.abs(zs - abl), np.abs(zs - lcl)) / np.where(abl > 0, abl, np.nan)
    return np.where(np.isfinite(zs) & (abl > 0), sd, np.nan)


# --- shear height: three implementations (see pyroconvection_type) -------------
#
# The distance from the ABL/LCL to the height of maximum wind shear is the fifth
# diagnostic of the operational ladder. On a coarse pressure-level model it cannot be
# honestly located (three per-column implementations below, one picked by data
# richness); on the ICON-EU model levels it can, via shear_height_grid above.

def shear_height_window(height_m, wind_u, wind_v, *, zmin: float = 200.0,
                        zmax: float = 6000.0, half_window_m: float = 100.0,
                        dz: float = 50.0) -> float:
    """VARIANT 1 -- height (m AGL) of maximum vector shear, moving-window method.

    Interpolates one column onto a uniform ``dz`` grid and returns the height where
    ``|dU/dz|`` over a ``2 * half_window_m`` window is largest. This is the reference
    method and reproduces the operational product's numbers.

    **It requires a genuinely resolved profile.** Given only a handful of pressure
    levels the interpolation invents the sub-window structure, and the "maximum"
    then reflects the interpolant rather than the flow; see :data:`_SHEAR_MIN_LEVELS`
    and prefer :func:`shear_height_adaptive` for model data of unknown richness.
    """
    z = np.asarray(height_m, float)
    u = np.asarray(wind_u, float)
    v = np.asarray(wind_v, float)
    good = np.isfinite(z) & np.isfinite(u) & np.isfinite(v)
    if good.sum() < 2:
        return float("nan")
    z, u, v = z[good], u[good], v[good]
    order = np.argsort(z)
    z, u, v = z[order], u[order], v[order]

    grid = np.arange(max(z[0], 0.0), min(z[-1], zmax) + dz, dz)
    if grid.size < 3:
        return float("nan")
    ug = np.interp(grid, z, u)
    vg = np.interp(grid, z, v)
    half = max(1, int(round(half_window_m / dz)))
    if grid.size <= 2 * half:
        return float("nan")
    span = grid[2 * half:] - grid[:-2 * half]
    du = ug[2 * half:] - ug[:-2 * half]
    dv = vg[2 * half:] - vg[:-2 * half]
    shear = np.hypot(du, dv) / np.maximum(span, 1e-6)
    centre = grid[half:-half]
    ok = (centre >= zmin) & (centre <= zmax) & np.isfinite(shear)
    if not ok.any():
        return float("nan")
    return float(centre[ok][np.nanargmax(shear[ok])])


def shear_height_none(*_args, **_kwargs) -> float:
    """VARIANT 2 -- decline to estimate a shear height; always ``nan``.

    The honest choice when the vertical resolution cannot locate a shear maximum.
    Feeding ``nan`` to :func:`pyroconvection_type` drops the shear test from the
    ladder (see :func:`pyroconvection_type_noshear`) rather than conditioning the
    classification on a number the data do not contain.
    """
    return float("nan")


def shear_height_adaptive(height_m, wind_u, wind_v, *,
                          min_levels: int = _SHEAR_MIN_LEVELS, **kwargs) -> float:
    """VARIANT 3 (pipeline default) -- window method when resolved, else ``nan``.

    Counts the finite levels in the column and calls :func:`shear_height_window`
    only when there are at least ``min_levels`` of them; otherwise returns ``nan``
    so the classifier degrades to the no-shear ladder. This lets one code path serve
    both a 6-level ICON-2I column (no shear test) and a 19-level ERA5 or sounding
    column (full ladder), and makes which path ran an inspectable property of the
    output rather than a silent assumption.
    """
    z = np.asarray(height_m, float)
    u = np.asarray(wind_u, float)
    v = np.asarray(wind_v, float)
    n = int((np.isfinite(z) & np.isfinite(u) & np.isfinite(v)).sum())
    if n < min_levels:
        return float("nan")
    return shear_height_window(z, u, v, **kwargs)


def shear_distance_ratio(shear_height_m, abl_m, lcl_m) -> float:
    """Distance from the shear maximum to the nearer of the ABL top / LCL, /ABL.

    The ladder's ``shear_close`` test: a shear layer sitting on the ABL top or the
    cloud base couples to the plume, one far above it does not.
    """
    if not np.isfinite([shear_height_m, abl_m, lcl_m]).all() or abl_m <= 0:
        return float("nan")
    return float(min(abs(shear_height_m - abl_m),
                     abs(shear_height_m - lcl_m)) / abl_m)


# --- the ladder ----------------------------------------------------------------

def pyroconvection_score(*, lcl_abl_ratio: float, ml_theta_gradient: float,
                         gamma_theta: float, rh_top_abl: float,
                         shear_distance: float | None = None,
                         thresholds: PyroconvThresholds = DEFAULT_PYROCONV_THRESHOLDS
                         ) -> float:
    """Continuous 0-100 pyroconvective-favourability score for a column.

    A weighted blend of the ladder's diagnostics: LCL/ABL proximity (30), mixed-layer
    instability (25), cap weakness (20), shear proximity (15) and moisture at the ABL
    top (10). When ``shear_distance`` is ``nan``/``None`` its 15 points are
    redistributed across the other four terms in proportion, so scores stay on the
    same 0-100 scale whichever ladder ran.

    This is a *diagnostic* ordering, not a calibrated probability of occurrence --
    it ranks columns, it does not tell you how likely a pyroCb is.
    """
    th = thresholds
    if not np.isfinite([lcl_abl_ratio, ml_theta_gradient, gamma_theta, rh_top_abl]).all():
        return float("nan")

    lcl_term = float(np.clip(1.0 - abs(lcl_abl_ratio - 1.0), 0.0, 1.0))
    stability_term = float(np.clip(
        (th.ml_overshoot_max - ml_theta_gradient) / (th.ml_overshoot_max + 1.0e-3),
        0.0, 1.0))
    cap_term = float(np.clip(
        (th.gamma_strong_cap - gamma_theta) / (th.gamma_strong_cap - 2.5e-3), 0.0, 1.0))
    moisture_term = float(np.clip((rh_top_abl - 55.0) / 40.0, 0.0, 1.0))

    terms = [(30.0, lcl_term), (25.0, stability_term), (20.0, cap_term),
             (10.0, moisture_term)]
    has_shear = shear_distance is not None and np.isfinite(shear_distance)
    if has_shear:
        terms.append((15.0, float(np.clip(1.0 - shear_distance, 0.0, 1.0))))
    total_w = sum(w for w, _ in terms)
    score = sum(w * t for w, t in terms) * (100.0 / total_w)
    return float(np.clip(score, 0.0, 100.0))


def pyroconvection_type_castellnou(*, lcl_abl_ratio: float, ml_theta_gradient: float,
                                   gamma_theta: float,
                                   fireline_intensity_kw: float | None = None,
                                   fli_threshold_kw: float = _FLI_PYROCU_KW) -> str:
    """THREE-DIAGNOSTIC ladder -- Castellnou et al. (2022) Fig. 7 as published.

    The original pyflam classifier, kept as the floor of :func:`pyroconvection_type`:
    it is the only ladder that runs on the variables the paper's own Table 1 reports,
    so it is what the reference cases are validated against. It uses no moisture or
    shear information, and (unlike the profile ladders) lets a weak cap promote to
    deep pyroCb from any LCL/ABL ratio, which is why it saturates on a well-mixed
    summer afternoon.

    1. **Fire-power gate.** ``fireline_intensity_kw`` below ``fli_threshold_kw``
       (1e4 kW/m, Tedim et al. 2018) -> ``surface_plume`` whatever the atmosphere.
       ``None`` reports the *potential* type, assuming a pyroCu-capable fire.
    2. **Mixed-layer stability.** Stable ML (> 1.1e-3 K/m) -> ``convection_plume`` (T21).
    3. **LCL/ABL ratio.** ratio > 1 -> ``overshooting_pyrocu`` (SCQ32);
       ratio < 1 -> ``resilient_pyrocu`` (M11).
    4. **Cap.** ``gamma_theta`` < 4.0e-3 -> ``deep_pyrocu_pyrocb`` (SCQ51 3.9e-3 deep
       vs M11 4.2e-3 resilient; SCQ41 5.1e-3 strong cap inhibits).
    """
    if fireline_intensity_kw is not None and fireline_intensity_kw < fli_threshold_kw:
        return "surface_plume"
    if ml_theta_gradient > _ML_STABLE:
        return "convection_plume"
    base = "overshooting_pyrocu" if lcl_abl_ratio > 1.0 else "resilient_pyrocu"
    if gamma_theta < _GAMMA_DEEP:
        return "deep_pyrocu_pyrocb"
    return base


def pyroconvection_type_shear(*, lcl_abl_ratio: float, ml_theta_gradient: float,
                              gamma_theta: float, rh_top_abl: float,
                              shear_distance: float,
                              fireline_intensity_kw: float | None = None,
                              fli_threshold_kw: float = _FLI_PYROCU_KW,
                              thresholds: PyroconvThresholds = DEFAULT_PYROCONV_THRESHOLDS
                              ) -> str:
    """VARIANT 1 -- FIVE-DIAGNOSTIC ladder (LCL/ABL, ML theta, cap, RH-top, shear).

    The full operational ladder. Relative to :func:`pyroconvection_type_castellnou`
    it is markedly more conservative about the top of the ladder: reaching
    ``deep_pyrocu_pyrocb`` needs the LCL *near* the ABL top (ratio <= 1.10), a
    neutral/unstable ML, a weak cap, moist air at the ABL top, **and** a shear
    maximum close to the plume -- not a weak cap alone.

    Requires a shear height, so it requires a resolved profile (>= ~10 levels);
    on 6-level model output use :func:`pyroconvection_type` (variant 3), which
    falls back rather than trusting an interpolated shear maximum.
    """
    th = thresholds
    if fireline_intensity_kw is not None and fireline_intensity_kw < fli_threshold_kw:
        return "surface_plume"
    required = [lcl_abl_ratio, ml_theta_gradient, gamma_theta, rh_top_abl, shear_distance]
    if not np.isfinite(required).all():
        return "surface_plume"

    score = pyroconvection_score(
        lcl_abl_ratio=lcl_abl_ratio, ml_theta_gradient=ml_theta_gradient,
        gamma_theta=gamma_theta, rh_top_abl=rh_top_abl,
        shear_distance=shear_distance, thresholds=th)

    unstable = ml_theta_gradient <= th.ml_stable
    slightly_stable_or_better = ml_theta_gradient <= th.ml_overshoot_max
    weak_cap = gamma_theta <= th.gamma_weak_cap
    strong_cap = gamma_theta >= th.gamma_strong_cap
    shear_close = shear_distance <= th.shear_distance_deep
    moist_top = rh_top_abl >= th.rh_top_moist

    if (lcl_abl_ratio <= th.lcl_ratio_deep_max and unstable and weak_cap
            and shear_close and moist_top):
        return "deep_pyrocu_pyrocb"
    if (lcl_abl_ratio < th.lcl_ratio_resilient_max and unstable and moist_top
            and (strong_cap or not shear_close or not weak_cap)):
        return "resilient_pyrocu"
    if 1.0 <= lcl_abl_ratio <= th.lcl_ratio_overshoot_max and slightly_stable_or_better:
        return "overshooting_pyrocu"
    if lcl_abl_ratio > 1.0:
        return "convection_plume"
    # LCL below the ABL but the column is not coherent/moist enough for a
    # persistent cloud: keep a cautious class only if it is at least well mixed.
    if unstable and np.isfinite(score) and score >= th.residual_score:
        return "overshooting_pyrocu"
    return "surface_plume"


def pyroconvection_type_noshear(*, lcl_abl_ratio: float, ml_theta_gradient: float,
                                gamma_theta: float, rh_top_abl: float,
                                fireline_intensity_kw: float | None = None,
                                fli_threshold_kw: float = _FLI_PYROCU_KW,
                                thresholds: PyroconvThresholds = DEFAULT_PYROCONV_THRESHOLDS
                                ) -> str:
    """VARIANT 2 -- FOUR-DIAGNOSTIC ladder: variant 1 with the shear test removed.

    For data whose vertical resolution cannot locate a shear maximum. The ladder is
    otherwise identical, with the shear clauses dropped: ``deep_pyrocu_pyrocb`` needs
    LCL/ABL <= 1.10, a neutral/unstable ML, a weak cap and a moist ABL top; the
    ``resilient_pyrocu`` fallback keys on the cap alone.

    Dropping the shear requirement makes class 4 *easier* to reach than in variant 1
    (one necessary condition fewer), so this ladder is the more permissive of the
    two at the top end. It remains far stricter than the Castellnou ladder, which
    has neither the moisture nor the LCL-proximity requirement.
    """
    th = thresholds
    if fireline_intensity_kw is not None and fireline_intensity_kw < fli_threshold_kw:
        return "surface_plume"
    if not np.isfinite([lcl_abl_ratio, ml_theta_gradient, gamma_theta, rh_top_abl]).all():
        return "surface_plume"

    score = pyroconvection_score(
        lcl_abl_ratio=lcl_abl_ratio, ml_theta_gradient=ml_theta_gradient,
        gamma_theta=gamma_theta, rh_top_abl=rh_top_abl, thresholds=th)

    unstable = ml_theta_gradient <= th.ml_stable
    slightly_stable_or_better = ml_theta_gradient <= th.ml_overshoot_max
    weak_cap = gamma_theta <= th.gamma_weak_cap
    strong_cap = gamma_theta >= th.gamma_strong_cap
    moist_top = rh_top_abl >= th.rh_top_moist

    if lcl_abl_ratio <= th.lcl_ratio_deep_max and unstable and weak_cap and moist_top:
        return "deep_pyrocu_pyrocb"
    if (lcl_abl_ratio < th.lcl_ratio_resilient_max and unstable and moist_top
            and (strong_cap or not weak_cap)):
        return "resilient_pyrocu"
    if 1.0 <= lcl_abl_ratio <= th.lcl_ratio_overshoot_max and slightly_stable_or_better:
        return "overshooting_pyrocu"
    if lcl_abl_ratio > 1.0:
        return "convection_plume"
    if unstable and np.isfinite(score) and score >= th.residual_score:
        return "overshooting_pyrocu"
    return "surface_plume"


PYROCONVECTION_LADDERS = ("castellnou", "noshear", "shear", "adaptive")


def pyroconvection_type(*, lcl_abl_ratio: float, ml_theta_gradient: float,
                        gamma_theta: float, ladder: str = "castellnou",
                        rh_top_abl: float | None = None,
                        shear_distance: float | None = None,
                        fireline_intensity_kw: float | None = None,
                        fli_threshold_kw: float = _FLI_PYROCU_KW,
                        thresholds: PyroconvThresholds = DEFAULT_PYROCONV_THRESHOLDS
                        ) -> str:
    """Classify a column into a pyroconvection prototype. Returns a :data:`PYROCONVECTION_TYPES` member.

    ``ladder`` selects the decision ladder. They are *not* interchangeable -- they
    condition on different variables and differ in strictness -- so the choice is
    explicit rather than inferred:

    * ``"castellnou"`` (**default**) -- :func:`pyroconvection_type_castellnou`, the
      three-diagnostic ladder exactly as published (LCL/ABL, mixed-layer dtheta/dz,
      cap gamma-theta). The only ladder validated against the paper's Table 1 cases,
      and the only one that runs without moisture or shear input.
    * ``"noshear"`` -- :func:`pyroconvection_type_noshear`, four diagnostics; also
      needs ``rh_top_abl``.
    * ``"shear"`` -- :func:`pyroconvection_type_shear`, the full five; also needs
      ``shear_distance``, hence a profile with enough levels to locate a shear
      maximum (>= ~10; see :func:`shear_height_adaptive`).
    * ``"adaptive"`` -- run the richest ladder the supplied diagnostics support:
      five if ``shear_distance`` is finite, else four if ``rh_top_abl`` is, else
      fall back to ``"castellnou"``. Intended for pipelines that must serve both a
      6-level model column and a full sounding from one call site; record which tier
      ran, since it changes what the class means.

    The profile ladders demand a moist ABL top and an LCL close to it before allowing
    deep pyroCb, where the Castellnou ladder promotes on a weak cap alone -- so they
    are markedly more conservative at the top of the ladder.

    All ladders share the fire-power gate: ``fireline_intensity_kw`` below
    ``fli_threshold_kw`` (1e4 kW/m, Tedim et al. 2018) returns ``surface_plume``
    whatever the atmosphere; ``None`` reports the *potential* type, assuming a
    pyroCu-capable fire.

    Raises ``ValueError`` on an unknown ladder, or on one whose required diagnostics
    are missing -- a silent downgrade would misreport what the class means. Use
    ``"adaptive"`` if you want the fallback.
    """
    if ladder not in PYROCONVECTION_LADDERS:
        raise ValueError(f"unknown ladder {ladder!r}; "
                         f"expected one of {PYROCONVECTION_LADDERS}")
    has_rh = rh_top_abl is not None and np.isfinite(rh_top_abl)
    has_shear = shear_distance is not None and np.isfinite(shear_distance)

    if ladder == "adaptive":
        ladder = "shear" if (has_rh and has_shear) else "noshear" if has_rh else "castellnou"
    if ladder in ("noshear", "shear") and not has_rh:
        raise ValueError(f"ladder={ladder!r} needs a finite rh_top_abl")
    if ladder == "shear" and not has_shear:
        raise ValueError("ladder='shear' needs a finite shear_distance")

    if ladder == "shear":
        return pyroconvection_type_shear(
            lcl_abl_ratio=lcl_abl_ratio, ml_theta_gradient=ml_theta_gradient,
            gamma_theta=gamma_theta, rh_top_abl=rh_top_abl,
            shear_distance=shear_distance, fireline_intensity_kw=fireline_intensity_kw,
            fli_threshold_kw=fli_threshold_kw, thresholds=thresholds)
    if ladder == "noshear":
        return pyroconvection_type_noshear(
            lcl_abl_ratio=lcl_abl_ratio, ml_theta_gradient=ml_theta_gradient,
            gamma_theta=gamma_theta, rh_top_abl=rh_top_abl,
            fireline_intensity_kw=fireline_intensity_kw,
            fli_threshold_kw=fli_threshold_kw, thresholds=thresholds)
    return pyroconvection_type_castellnou(
        lcl_abl_ratio=lcl_abl_ratio, ml_theta_gradient=ml_theta_gradient,
        gamma_theta=gamma_theta, fireline_intensity_kw=fireline_intensity_kw,
        fli_threshold_kw=fli_threshold_kw)


# --- integration with the fire model ------------------------------------------

def midflame_wind_ft_per_min(state: AtmosphericState, *,
                             wind_reduction_factor: float = 0.4) -> float:
    """Midflame wind (ft/min) from the 10 m wind via a wind reduction factor."""
    # 10 m wind ~ 20-ft wind (6.1 m); the WRF reduces it to midflame.
    return m_per_s_to_ft_per_min(state.wind_speed) * wind_reduction_factor


def spread_inputs_from_state(state: AtmosphericState, *,
                             wind_reduction_factor: float = 0.4,
                             moisture_offsets=(0.0, 1.0, 2.0)) -> dict:
    """All scenario inputs for :func:`pyflam.spread` / ``spread_field`` from a state.

    Bundles dead fuel moisture (from T/RH), midflame wind (ft/min) and wind
    direction (deg FROM). Add your own live-fuel moistures and ``load_factor``.
    """
    out = dead_fuel_moisture(state, offsets=moisture_offsets)
    out["wind_midflame"] = midflame_wind_ft_per_min(
        state, wind_reduction_factor=wind_reduction_factor)
    out["wind_direction"] = state.wind_direction
    return out


def atmospheric_firebrand_physics(state: AtmosphericState, base=None):
    """A :class:`pyflam.spotting.FirebrandPhysics` modulated by the convection.

    Scales the plume length scale by :func:`convective_plume_factor`, so spotting
    reaches farther in an unstable, high-CAPE atmosphere and is suppressed in a
    stable one.
    """
    from .spotting import FirebrandPhysics
    base = base or FirebrandPhysics()
    return replace(base, front_length=base.front_length
                   * convective_plume_factor(state))


# --- providers ----------------------------------------------------------------

def latlon_grid(ls):
    """Cell-centre latitude/longitude (2D arrays) for a landscape, or ``None``.

    Returns ``(lat2d, lon2d)`` if the landscape's CRS can be resolved to
    geographic coordinates (directly when already geographic, else via ``pyproj``
    when projected), or ``None`` if there is no usable CRS -- the caller then
    falls back to a single representative point.
    """
    nrows, ncols = ls.shape
    cols = ls.west + (np.arange(ncols) + 0.5) * ls.cellsize_x
    rows = ls.north - (np.arange(nrows) + 0.5) * ls.cellsize_y
    xx, yy = np.meshgrid(cols, rows)
    crs = ls.crs
    if crs is None:
        return None
    try:
        from pyproj import CRS, Transformer
        crs = CRS.from_user_input(crs)
        if crs.is_geographic:
            return yy, xx                         # already lon/lat in x/y
        tr = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
        lon2d, lat2d = tr.transform(xx, yy)
        return np.asarray(lat2d), np.asarray(lon2d)
    except Exception:
        return None


def _broadcast_state(state: AtmosphericState, shape) -> AtmosphericState:
    """Broadcast a scalar state's numeric fields to 2D arrays of ``shape``."""
    kw = {}
    for f in ("wind_speed", "wind_direction", "temperature", "relative_humidity",
              "pressure", "sensible_heat_flux", "latent_heat_flux", "cape",
              "cin", "boundary_layer_height"):
        v = getattr(state, f)
        kw[f] = None if v is None else np.full(shape, float(v))
    return replace(state, **kw)


def wind_field_from_state(state: AtmosphericState, ls):
    """A :class:`~pyflam.wind.WindField` (m/s) from a (gridded or scalar) state."""
    from .wind import WindField
    spd = np.broadcast_to(np.asarray(state.wind_speed, float), ls.shape).copy()
    dirn = np.broadcast_to(np.asarray(state.wind_direction, float), ls.shape).copy()
    return WindField(speed=spd, direction=dirn, cellsize=ls.cellsize_x,
                     west=ls.west, north=ls.north, speed_units="m/s", crs=ls.crs)


class AtmosphereProvider:
    """Returns an :class:`AtmosphericState` for a location and time."""

    def state_at(self, latitude: float, longitude: float,
                 time: datetime | None = None) -> AtmosphericState:
        raise NotImplementedError

    def field_on(self, ls, time=None, *, latlon=None) -> AtmosphericState:
        """Sample the state onto a landscape grid (a *field*: array-valued state).

        Default: sample at the landscape centre and broadcast (uniform over the
        domain). Gridded providers override this to vary per cell. ``latlon`` may
        supply precomputed ``(lat2d, lon2d)`` arrays.
        """
        if latlon is None:
            latlon = latlon_grid(ls)
        if latlon is None:
            lat = lon = None
        else:
            lat = float(np.mean(latlon[0]))
            lon = float(np.mean(latlon[1]))
        return _broadcast_state(self.state_at(lat, lon, time), ls.shape)


class ConstantAtmosphere(AtmosphereProvider):
    """A provider that returns the same state everywhere (testing / idealized)."""

    def __init__(self, state: AtmosphericState):
        self.state = state

    def state_at(self, latitude=None, longitude=None, time=None):
        return replace(self.state, latitude=latitude, longitude=longitude,
                       time=time or self.state.time)

    def field_on(self, ls, time=None, *, latlon=None):
        return _broadcast_state(replace(self.state, time=time or self.state.time),
                                ls.shape)


# Canonical variable -> dataset-variable-name maps for common sources. Override
# per dataset as needed; only the names present are read.
ERA5_VARS = {
    "wind_u": "u10", "wind_v": "v10", "temperature_K": "t2m",
    "dewpoint_K": "d2m", "pressure": "sp", "sensible_heat_flux": "sshf",
    "latent_heat_flux": "slhf", "cape": "cape", "boundary_layer_height": "blh",
}
GFS_VARS = {
    "wind_u": "u10", "wind_v": "v10", "temperature_K": "t2m",
    "relative_humidity": "r2", "pressure": "sp", "cape": "cape", "cin": "cin",
    # HPBL is undecodable by some eccodes builds (shortName 'unknown'); fetch_gfs
    # fetches it alone and renames the sole field to 'hpbl' (see fetch_gfs).
    "boundary_layer_height": "hpbl",
    # GFS surface heat fluxes decode as time-mean fields and are positive upward
    # (surface heating the air) -- pyflam's convention -- so no transform is needed.
    # Only present in forecast files (fxx >= 3); absent from the f000 analysis.
    "sensible_heat_flux": "avg_ishf", "latent_heat_flux": "avg_slhtf",
}

# ERA5 surface heat fluxes are *accumulated* (J/m^2) over the product step and
# positive **downward** (into the surface). pyflam wants an instantaneous flux in
# W/m^2 positive **upward** (surface heating the air), so divide by the
# accumulation period and flip the sign. Hourly ERA5 accumulates over 1 hour.
ERA5_ACCUMULATION_SECONDS = 3600.0


def era5_flux_to_watts(accumulated_j_m2, accumulation_seconds=ERA5_ACCUMULATION_SECONDS):
    """ERA5 accumulated flux (J/m^2, down) -> instantaneous W/m^2 (up)."""
    return -np.asarray(accumulated_j_m2, dtype=float) / float(accumulation_seconds)


# Per-variable post-processing applied after sampling (canonical name -> fn).
ERA5_TRANSFORMS = {
    "sensible_heat_flux": era5_flux_to_watts,
    "latent_heat_flux": era5_flux_to_watts,
}


def _naive_utc(time):
    """A tz-aware datetime -> naive UTC, so it compares with naive ``datetime64``.

    Gridded datasets (ERA5/GFS/ICON NetCDF) carry tz-naive time coordinates;
    xarray's ``.sel`` cannot compare those with a tz-aware label. Callers may pass
    tz-aware datetimes, so normalise to naive UTC before selecting.
    """
    if time is not None and getattr(time, "tzinfo", None) is not None:
        return time.astimezone(timezone.utc).replace(tzinfo=None)
    return time


class GriddedAtmosphere(AtmosphereProvider):
    """Provider backed by a gridded dataset (xarray): WRF / ERA5 / GFS / NetCDF.

    ``dataset`` is an ``xarray.Dataset`` with latitude/longitude (and optional
    time) coordinates; ``var_map`` maps canonical names (see :data:`ERA5_VARS`)
    to its variables. ``state_at`` selects the nearest grid point (and time).
    Build one from a downloaded file with :func:`open_atmosphere`.
    """

    def __init__(self, dataset, var_map: dict, *, transforms=None,
                 lat_name="latitude", lon_name="longitude", time_name="time"):
        self.ds = dataset
        self.var_map = var_map
        self.transforms = transforms or {}
        self.lat_name, self.lon_name, self.time_name = lat_name, lon_name, time_name

    def _apply(self, canon, value):
        fn = self.transforms.get(canon)
        return fn(value) if fn is not None else value

    def _wrap_lon(self, lon):
        """Match the dataset's longitude convention (GFS is 0-360, ERA5 -180-180)."""
        lons = np.asarray(self.ds[self.lon_name].values)
        if lons.size and float(np.nanmax(lons)) > 180.0:
            return np.asarray(lon, dtype=float) % 360.0
        return lon

    def state_at(self, latitude, longitude, time=None):
        sel = {self.lat_name: latitude, self.lon_name: self._wrap_lon(longitude)}
        point = self.ds.sel(**sel, method="nearest")
        # Only select on time when it is an indexable dimension; a single-time
        # field carries time as a scalar coordinate that cannot be `.sel`-ed.
        if time is not None and self.time_name in self.ds.dims:
            point = point.sel({self.time_name: _naive_utc(time)}, method="nearest")
        kw = {}
        for canon, var in self.var_map.items():
            if var in point:
                v = float(np.asarray(point[var].values).reshape(-1)[0])
                kw[canon] = float(np.asarray(self._apply(canon, v)).reshape(-1)[0])
        return AtmosphericState.from_si(
            latitude=latitude, longitude=longitude, time=time, **kw)

    def field_on(self, ls, time=None, *, latlon=None):
        """Sample every landscape cell from the gridded dataset (per-cell state)."""
        import xarray as xr
        if latlon is None:
            latlon = latlon_grid(ls)
        if latlon is None:                        # no georeferencing -> uniform
            return super().field_on(ls, time)
        lat2d, lon2d = latlon
        lat_da = xr.DataArray(np.asarray(lat2d).ravel(), dims="p")
        lon_da = xr.DataArray(self._wrap_lon(np.asarray(lon2d).ravel()), dims="p")
        pts = self.ds.sel({self.lat_name: lat_da, self.lon_name: lon_da},
                          method="nearest")
        if time is not None and self.time_name in self.ds.dims:
            pts = pts.sel({self.time_name: _naive_utc(time)}, method="nearest")
        shape = ls.shape
        kw = {}
        for canon, var in self.var_map.items():
            if var in pts:
                v = np.asarray(pts[var].values, dtype=float).reshape(shape)
                kw[canon] = np.asarray(self._apply(canon, v), dtype=float)
        return AtmosphericState.from_si(time=time, **kw)


def _read_atmosphere_dataset(path, **kwargs):
    """Open a forecast/reanalysis file as a single xarray ``Dataset``.

    Handles the **zip** the new Copernicus CDS returns for ERA5 (instantaneous and
    accumulated variables in separate NetCDFs): extracts the members, merges them
    and loads the result into memory so the archive can be cleaned up.
    """
    import xarray as xr
    import zipfile

    if not zipfile.is_zipfile(path):
        return xr.open_dataset(path, **kwargs)

    import glob
    import os
    import shutil
    import tempfile

    tmp = tempfile.mkdtemp(prefix="pyflam_era5_")
    try:
        with zipfile.ZipFile(path) as z:
            z.extractall(tmp)
        members = sorted(glob.glob(os.path.join(tmp, "*.nc")))
        if not members:
            raise ValueError(f"no NetCDF member found inside archive {path!r}")
        parts = [xr.open_dataset(m, **kwargs) for m in members]
        try:
            ds = xr.merge(parts, compat="override", combine_attrs="override").load()
        finally:
            for p in parts:
                p.close()
        return ds
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def open_atmosphere(path: str, source: str = "era5", **kwargs) -> GriddedAtmosphere:
    """Open a downloaded forecast/reanalysis file as a provider (needs xarray).

    ``source`` selects the variable map (``"era5"``, ``"gfs"``, or pass
    ``var_map=`` for WRF/other). GRIB needs the ``cfgrib`` engine. This is the
    offline path -- download once (ERA5 via the Copernicus CDS, GFS via NOMADS)
    and point pyflam at the file. The new CDS returns ERA5 as a **zip** of
    NetCDFs; that is detected and the members are merged transparently.
    """
    try:
        import xarray as xr  # noqa: F401 -- ensures a clear error if xarray is absent
    except ImportError as exc:  # pragma: no cover - exercised only without xarray
        raise ImportError(
            "Atmospheric file reading needs xarray (and cfgrib for GRIB): "
            "pip install xarray cfgrib netcdf4"
        ) from exc
    var_map = kwargs.pop("var_map", None) or {
        "era5": ERA5_VARS, "gfs": GFS_VARS}.get(source, ERA5_VARS)
    transforms = kwargs.pop("transforms", None)
    if transforms is None and source == "era5":
        transforms = ERA5_TRANSFORMS          # accumulated J/m^2 down -> W/m^2 up
    ds = _read_atmosphere_dataset(path, **kwargs)
    # The new CDS ERA5 NetCDFs use 'valid_time' as the time coordinate; older
    # files use 'time'. Point the provider at whichever is present.
    gkw = {}
    if "time" not in ds.coords and "time" not in ds.dims and "valid_time" in ds.coords:
        gkw["time_name"] = "valid_time"
    return GriddedAtmosphere(ds, var_map, transforms=transforms, **gkw)


# Default ERA5 single-level variables for fire + convection (CDS long names).
ERA5_FIRE_VARIABLES = [
    "10m_u_component_of_wind", "10m_v_component_of_wind", "2m_temperature",
    "2m_dewpoint_temperature", "surface_pressure",
    "surface_sensible_heat_flux", "surface_latent_heat_flux",
    "convective_available_potential_energy", "boundary_layer_height",
]


def era5_request(*, date: str, time, area, variables=None) -> dict:
    """Build the Copernicus CDS retrieval dict for an ERA5 fire/convection query.

    ``date`` ``"YYYY-MM-DD"``; ``time`` an hour string ``"13:00"`` or a list;
    ``area`` ``(north, west, south, east)`` in degrees. Pure -- no network -- so
    it is unit-testable; :func:`fetch_era5` submits it.
    """
    return {
        "product_type": "reanalysis",
        "format": "netcdf",
        "variable": list(variables or ERA5_FIRE_VARIABLES),
        "date": date,
        "time": [time] if isinstance(time, str) else list(time),
        "area": [area[0], area[1], area[2], area[3]],
    }


def fetch_era5(cache_path: str, *, date: str, time, area, variables=None,
               force: bool = False) -> GriddedAtmosphere:
    """Download (and cache) an ERA5 slice from Copernicus and open it as a provider.

    Returns a :class:`GriddedAtmosphere`. If ``cache_path`` already exists and not
    ``force``, the download is skipped (caching for repeat / reanalysis runs).
    Needs the ``cdsapi`` package and CDS credentials (``~/.cdsapirc``).

    Note: ERA5 surface heat fluxes are time-accumulated (J/m^2) and positive
    downward; convert to W/m^2 upward for :class:`AtmosphericState` if you use
    them quantitatively (a documented post-processing step).
    """
    import os
    if force or not os.path.exists(cache_path):
        try:
            import cdsapi
        except ImportError as exc:  # pragma: no cover
            raise ImportError("ERA5 fetch needs cdsapi: pip install cdsapi "
                              "(and configure ~/.cdsapirc)") from exc
        cdsapi.Client().retrieve(
            "reanalysis-era5-single-levels",
            era5_request(date=date, time=time, area=area, variables=variables),
            cache_path)
    return open_atmosphere(cache_path, source="era5")


def fetch_gfs(*, run, fxx: int = 0, cache_dir: str | None = None,
              product: str = "pgrb2.0p25"):
    """Fetch a GFS forecast field and open it as a provider (near-real-time).

    ``run`` is the model run time (datetime or ``"YYYY-MM-DD HH:MM"``), ``fxx``
    the forecast hour. Uses Herbie (which caches downloads under ``cache_dir``).
    Needs ``herbie-data`` (and ``cfgrib``). Returns a :class:`GriddedAtmosphere`.

    Retrieves the surface + 2 m + 10 m fire-weather and convection / energy-flux
    fields: 10 m wind, 2 m T/RH, surface pressure, CAPE, CIN, boundary-layer
    height and the surface sensible/latent heat fluxes. The **heat fluxes are
    time-mean forecast fields and only exist for ``fxx >= 3``** (not in the f000
    analysis); request a forecast hour to drive the plume / energy-flux diagnostics.
    Boundary-layer height (HPBL) is fetched separately because some eccodes builds
    cannot decode its short name from a merged request.
    """
    try:
        from herbie import Herbie
    except ImportError as exc:  # pragma: no cover
        raise ImportError("GFS fetch needs Herbie: pip install herbie-data "
                          "cfgrib") from exc
    import warnings

    import xarray as xr
    herbie_kw = {} if cache_dir is None else {"save_dir": cache_dir}
    h = Herbie(run, model="gfs", product=product, fxx=fxx, **herbie_kw)

    def _merge(ds):
        if isinstance(ds, list):   # cfgrib returns one dataset per "hypercube"
            merged = ds[0]
            for d in ds[1:]:
                merged = xr.merge([merged, d], compat="override",
                                  combine_attrs="override")
            return merged
        return ds

    # Surface + 2 m + 10 m + convective / energy-flux fields. Heat fluxes
    # (SHTFL/LHTFL) and CIN/PRES are at the surface; the f000 analysis omits the
    # (time-mean) fluxes -- they simply won't be present and read back as None.
    search = (r":(UGRD|VGRD):10 m above ground:|:(TMP|RH):2 m above ground:"
              r"|:(CAPE|CIN|PRES):surface:|:(SHTFL|LHTFL):surface:")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")            # cfgrib / Herbie merge chatter
        ds = _merge(h.xarray(search, remove_grib=False))
        # HPBL alone -> the sole field is the boundary-layer height; rename it to
        # the canonical 'hpbl' even when eccodes leaves it as 'unknown'.
        try:
            hpbl = _merge(h.xarray(r":HPBL:surface:", remove_grib=False))
            name = list(hpbl.data_vars)[0]
            ds = xr.merge([ds, hpbl[name].rename("hpbl")],
                          compat="override", combine_attrs="override")
        except Exception:                          # HPBL optional -- keep going
            pass
    return GriddedAtmosphere(ds, GFS_VARS, lat_name="latitude",
                             lon_name="longitude")


# --- ICON-2I 2.2 km (Italy) from the MISTRAL / AgenziaItaliaMeteo open archive ---
#
# The Italian convection-permitting model ICON-2I (2.2 km, full national domain) is
# published as open data (CC-BY) on MISTRAL / AgenziaItaliaMeteo MeteoHub, one GRIB2
# file per variable and per level. This is the convection-permitting source for fire
# weather over Italy (ICON-D2 does not cover Italy; GFS/ICON-EU are coarser). The
# files carry the full forecast (out to ~72 h) on a 0.025 deg (~2.2 km) lat/lon grid.

ICON2I_MISTRAL_BASE = ("https://meteohub.agenziaitaliameteo.it/nwp/"
                       "ICON-2I_SURFACE_PRESSURE_LEVELS")

# (local_name, variable_subdir, typeOfLevel, level) -- the fields the pyroconvection
# classifier + fuel gate need: pressure-level T, 2 m T/dewpoint, 10 m wind.
ICON2I_FIRE_FIELDS = [
    ("T850", "T", "isobaricInhPa", 850), ("T700", "T", "isobaricInhPa", 700),
    ("T500", "T", "isobaricInhPa", 500),
    ("T2M", "T_2M", "heightAboveGround", 2), ("TD2M", "TD_2M", "heightAboveGround", 2),
    ("U10", "U_10M", "heightAboveGround", 10), ("V10", "V_10M", "heightAboveGround", 10),
]

# Pressure levels the ICON-2I open-data archive actually publishes. 250 hPa (~10 km)
# is above anything the ABL/cap diagnostics look at, so the profile stops at 500.
ICON2I_PROFILE_LEVELS = (1000, 925, 850, 700, 500)

# The full profile the Rib-based ABL + moisture ladder needs: geopotential (FI ->
# real per-cell heights, replacing the standard-atmosphere assumption), temperature,
# RH and wind on each level, plus the surface reference state, the surface pressure
# and orography (to get heights AGL) and the land fraction. RELHUM is taken straight
# from the archive rather than derived from QV -- one variable fewer to download.
ICON2I_PROFILE_FIELDS = (
    [(f"{var}{p}", var, "isobaricInhPa", p)
     for var in ("FI", "T", "RELHUM", "U", "V") for p in ICON2I_PROFILE_LEVELS]
    + [("T2M", "T_2M", "heightAboveGround", 2), ("TD2M", "TD_2M", "heightAboveGround", 2),
       ("U10", "U_10M", "heightAboveGround", 10), ("V10", "V_10M", "heightAboveGround", 10),
       ("PS", "PS", "surface", 0), ("HSURF", "HSURF", "surface", 0),
       ("FRLAND", "FR_LAND", "surface", 0)]
)


def fetch_icon2i_mistral(date, run: int = 0, *, cache_dir: str = ".",
                         fields=None, base_url: str = ICON2I_MISTRAL_BASE,
                         force: bool = False, timeout: int = 600) -> dict:
    """Download ICON-2I 2.2 km GRIB files from the MISTRAL/AgenziaItaliaMeteo archive.

    ``date`` is a ``datetime``/``date`` (the run day); ``run`` the run hour (0 or 12).
    Downloads one GRIB per field in ``fields`` (default :data:`ICON2I_FIRE_FIELDS` --
    the fields the pyroconvection map needs) into ``cache_dir``, skipping any already
    present unless ``force``. Returns ``{local_name: path}``; open each with
    ``xarray.open_dataset(path, engine="cfgrib")`` (one variable per file).

    The archive directory is open (no token) for this NWP product; the broader
    MeteoHub extraction API is separate. File naming follows
    ``ICON_2I_SURFACE_PRESSURE_LEVELS_{YYYYMMDDHH}_{typeOfLevel}-{level}.grib`` under
    a per-variable subdirectory and a ``{YYYYMMDDHH}`` run directory.
    """
    import os
    import urllib.request

    stamp = f"{date:%Y%m%d}{int(run):02d}"
    fields = fields or ICON2I_FIRE_FIELDS
    os.makedirs(cache_dir, exist_ok=True)
    out = {}
    for name, subdir, leveltype, level in fields:
        fn = f"ICON_2I_SURFACE_PRESSURE_LEVELS_{stamp}_{leveltype}-{level}.grib"
        url = f"{base_url}/{stamp}/{subdir}/{fn}"
        dst = os.path.join(cache_dir, f"{name}.grib")
        if force or not os.path.exists(dst) or os.path.getsize(dst) == 0:
            urllib.request.urlretrieve(url, dst)            # raises on HTTP error
        out[name] = dst
    return out


# --- ICON-EU (DWD open data): coarser horizontally (6.5 km) but with the native
# vertical grid, which the Italian 2.2 km product does not publish. Used for the
# hybrid pyroconvection product: profile diagnostics (ABL, mixed-layer stability,
# cap, LCL) from the model levels, fuel gate from ICON-2I's 2.2 km surface fields.
# See scripts/validation/: on the model levels the mixed-layer dtheta/dz becomes a
# genuine measurement (~10 levels inside the mixed layer, vs ~2 on ICON-2I).

ICON_EU_BASE = "https://opendata.dwd.de/weather/nwp/icon-eu/grib"

# ICON-EU has 74 model levels; level 74 is the lowest (~10 m AGL). HHL is a *half*-level
# height field, so the full level k needs half levels k and k+1 -- hence the +1 fetched below.
#
# Measured half-level heights over Tuscany (m MSL, terrain ~300 m): 51 = 4407, 48 = 5291,
# 45 = 6177, 42 = 7066, 39 = 7958, 36 = 8851 -- roughly 295 m per level in this range.
#
# The set ran 74..51 (surface to ~4 km) until 2026-07-27, which is ample for the ABL
# diagnostics but stops at about +2 C over Tuscany in summer. The PyroCb Firepower Threshold
# (:func:`pyrocb_firepower_threshold_grid`) needs the profile up to the **-20 C** cloud-top
# level that Tory & Kepert (2021) use as the electrification criterion -- near 7.4 km here --
# so it returned nan in every column on the shorter stack.
#
# 74..39 (~8 km MSL) proved too tight: over Tuscany on 2026-07-27 the stack top sat at
# -20.3 C at the median, so 29 % of land columns could not be tested at all -- and that 29 %
# matched the nan fraction exactly, i.e. every failure was "cannot see high enough" rather
# than "cannot make pyroCb". Extended again to 74..36 (~8.9 km MSL, ~-26 C) for real margin.
#
# Cost: 39 levels rather than 24, so ~235 MB per step against ~148 MB, or ~5.6 GB for a
# 3-day 8-hourly run. Callers that only need the ABL diagnostics can pass the shorter range
# to :func:`fetch_icon_eu` and skip the upper fifteen.
ICON_EU_MODEL_LEVELS = tuple(range(74, 35, -1))     # ascending in height (74 -> 36)
ICON_EU_MODEL_VARS = ("T", "QV", "U", "V", "P")
ICON_EU_SURFACE_VARS = ("T_2M", "TD_2M", "PS", "U_10M", "V_10M")


def fetch_icon_eu(date, run: int = 0, step: int = 0, *, cache_dir: str = ".",
                  levels=ICON_EU_MODEL_LEVELS, model_vars=ICON_EU_MODEL_VARS,
                  surface_vars=ICON_EU_SURFACE_VARS, base_url: str = ICON_EU_BASE,
                  force: bool = False, timeout: int = 600) -> dict:
    """Download one ICON-EU forecast step (model levels + surface) from DWD open data.

    ``date`` is the run day, ``run`` the run hour (00/03/.../21 -- ICON-EU has 8 runs
    a day), ``step`` the forecast lead in hours. Downloads, for that step, each
    ``model_vars`` variable on every model level in ``levels`` (default the lowest 24,
    surface to ~4 km), the ``surface_vars`` single-level fields, the ``HHL`` half-level
    heights and ``FR_LAND`` (both time-invariant). Files are bz2 on the server and are
    stored decompressed in ``cache_dir``; anything already present is skipped unless
    ``force``. Returns ``{local_name: path}`` -- ``"{VAR}{level}"`` for model levels,
    ``"HHL{level}"`` for the half-level heights, and the bare variable name otherwise.

    One step is ~222 MB (36 levels x 5 vars). The daily hybrid product fetches 8 steps
    per run (~1.6 GB), against ~2.6 GB for the ICON-2I product; ICON-EU serves one small
    file per level/step rather than whole-domain blobs.
    """
    import bz2
    import os
    import urllib.request

    stamp = f"{date:%Y%m%d}{int(run):02d}"
    os.makedirs(cache_dir, exist_ok=True)

    def grab(url, dst):
        if force or not os.path.exists(dst) or os.path.getsize(dst) == 0:
            tmp, _ = urllib.request.urlretrieve(url)        # raises on HTTP error
            with open(tmp, "rb") as f:
                raw = bz2.decompress(f.read())
            os.remove(tmp)
            with open(dst, "wb") as f:
                f.write(raw)

    out = {}
    ml = f"{base_url}/{int(run):02d}"
    for var in model_vars:
        for lev in levels:
            fn = (f"icon-eu_europe_regular-lat-lon_model-level_{stamp}_"
                  f"{int(step):03d}_{lev}_{var}.grib2")
            dst = os.path.join(cache_dir, f"{var}{lev}.grib2")
            grab(f"{ml}/{var.lower()}/{fn}.bz2", dst)
            out[f"{var}{lev}"] = dst
    # HHL half levels: full level k lies between half levels k and k+1.
    for lev in tuple(levels) + (max(levels) + 1,):
        fn = f"icon-eu_europe_regular-lat-lon_time-invariant_{stamp}_{lev}_HHL.grib2"
        dst = os.path.join(cache_dir, f"HHL{lev}.grib2")
        grab(f"{ml}/hhl/{fn}.bz2", dst)
        out[f"HHL{lev}"] = dst
    for var in surface_vars:
        fn = (f"icon-eu_europe_regular-lat-lon_single-level_{stamp}_"
              f"{int(step):03d}_{var}.grib2")
        dst = os.path.join(cache_dir, f"{var}.grib2")
        grab(f"{ml}/{var.lower()}/{fn}.bz2", dst)
        out[var] = dst
    fn = f"icon-eu_europe_regular-lat-lon_time-invariant_{stamp}_FR_LAND.grib2"
    dst = os.path.join(cache_dir, "FR_LAND.grib2")
    grab(f"{ml}/fr_land/{fn}.bz2", dst)
    out["FR_LAND"] = dst
    return out
