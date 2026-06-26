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
from datetime import datetime

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
        """Step the three classes toward the state's EMC and return the moistures."""
        emc = equilibrium_moisture_content(
            state.temperature, state.relative_humidity) / 100.0
        dt_h = dt_minutes / 60.0
        self.m_1h = float(time_lag_step(self.m_1h, emc, dt_h, FUEL_TIME_LAGS["m_1h"]))
        self.m_10h = float(time_lag_step(self.m_10h, emc, dt_h, FUEL_TIME_LAGS["m_10h"]))
        self.m_100h = float(time_lag_step(self.m_100h, emc, dt_h,
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
    the pyroCb literature favours over surface CAPE: a deep, dry mixed layer (high
    LCL) capped by moisture aloft (the inverted-V sounding) and elevated lower-
    tropospheric instability/dryness (Continuous Haines). Surface CAPE is a poor
    pyroCb predictor -- pyroCb routinely form with near-zero CAPE -- so when the
    profile shows a pyroconvection-favourable column the factor is boosted toward
    its cap (Castellnou et al. 2022; Peterson et al. 2017; Mills & McCaw 2010).
    """
    cape = state.cape or 0.0
    factor = 1.0 + cape / _CAPE_REF
    cls = stability_class(state)
    if cls == "stable":
        factor *= 0.7
    elif cls == "unstable":
        factor *= 1.2

    if profile is not None:
        ch = continuous_haines(profile)
        iv, _, _ = inverted_v(profile)
        # C-Haines ~13 max: ramp a multiplicative boost up to ~1.6x by CH=11,
        # plus an extra kick for a confirmed inverted-V (the canonical pyroCb set-up).
        factor *= 1.0 + 0.6 * np.clip((ch - 5.0) / 6.0, 0.0, 1.0)
        if iv:
            factor *= 1.25
    return float(np.clip(factor, 0.5, 3.0))


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
    "boundary_layer_height": "hpbl", "sensible_heat_flux": "shtfl",
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
            point = point.sel({self.time_name: time}, method="nearest")
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
            pts = pts.sel({self.time_name: time}, method="nearest")
        shape = ls.shape
        kw = {}
        for canon, var in self.var_map.items():
            if var in pts:
                v = np.asarray(pts[var].values, dtype=float).reshape(shape)
                kw[canon] = np.asarray(self._apply(canon, v), dtype=float)
        return AtmosphericState.from_si(time=time, **kw)


def open_atmosphere(path: str, source: str = "era5", **kwargs) -> GriddedAtmosphere:
    """Open a downloaded forecast/reanalysis file as a provider (needs xarray).

    ``source`` selects the variable map (``"era5"``, ``"gfs"``, or pass
    ``var_map=`` for WRF/other). GRIB needs the ``cfgrib`` engine. This is the
    offline path -- download once (ERA5 via the Copernicus CDS, GFS via NOMADS)
    and point pyflam at the file.
    """
    try:
        import xarray as xr
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
    ds = xr.open_dataset(path, **kwargs)
    return GriddedAtmosphere(ds, var_map, transforms=transforms)


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
    # Surface + 2 m + 10 m + convective fields relevant to fire behaviour.
    search = (r":(UGRD|VGRD):10 m above ground:|:(TMP|RH):2 m above ground:"
              r"|:CAPE:surface:|:HPBL:")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")            # cfgrib / Herbie merge chatter
        ds = h.xarray(search, remove_grib=False)
        if isinstance(ds, list):   # cfgrib returns one dataset per "hypercube"
            merged = ds[0]
            for d in ds[1:]:
                merged = xr.merge([merged, d], compat="override",
                                  combine_attrs="override")
            ds = merged
    return GriddedAtmosphere(ds, GFS_VARS, lat_name="latitude",
                             lon_name="longitude")
