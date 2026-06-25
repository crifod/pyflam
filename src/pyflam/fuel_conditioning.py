"""Per-cell dead fuel moisture conditioning (terrain insolation + canopy shading).

Turns gridded *or* scalar weather (air temperature, relative humidity) plus a
landscape's slope / aspect / elevation / canopy-cover bands into **per-cell**
1/10/100-h dead fuel moisture rasters -- the pyflam analog of FlamMap's "dead
fuel moisture conditioning" step. The output is the ``m_1h`` / ``m_10h`` /
``m_100h`` dict that :func:`pyflam.spread_field` consumes (it already accepts
per-cell moisture arrays, grouping cells by moisture internally), so conditioned
moisture flows straight into MTT growth and burn probability.

The dominant control is **solar exposure**. Sun-exposed fuels (south-facing in the
northern hemisphere, on steep slopes tilted toward the sun, in the open) absorb
more shortwave, run warmer than the 2 m air, and therefore sit at a *lower*
equilibrium moisture content (EMC) -- they are drier. Shaded fuels (north-facing,
under canopy, or at night) stay near the ambient air EMC and are moister. This is
the mechanism documented for complex terrain by Holden et al. (2011): aspect and
canopy cover produce substantial, fire-relevant variation in dead fuel moisture
that a flat-terrain weather model misses entirely.

The scheme, per cell:

1. A dimensionless **sun-exposure index** ``S in [0, 1]`` = (beam factor on the
   slope facet) x (fraction of canopy openness). ``S = 1`` is sun perpendicular to
   fully open fuel (maximum heating); ``S = 0`` is no direct beam (night, deep
   shade, or a slope facing away from a low sun).
2. A near-fuel temperature rise ``dT = dt_sun * S`` above the air temperature.
3. Holding the near-surface water-vapour pressure fixed, the warmer fuel sees a
   lower local relative humidity, and a moisture submodel is evaluated at that
   warmer, drier microclimate -- either the NFDRS EMC
   (:func:`pyflam.atmosphere.equilibrium_moisture_content`, ``model="emc"``) or the
   semi-mechanistic VPD model (:func:`dead_fuel_moisture_vpd`, ``model="vpd"``;
   Resco de Dios 2015 / Nolan 2016), the more accurate point estimator for fine
   dead fuels. The slower 10-/100-h classes are offset upward by a small fixed
   amount (a simple, standard scheme; a full time-lag run with
   :class:`pyflam.atmosphere.DeadFuelMoistureModel` is the refinement).

This is a quasi-steady "conditioning" of the moisture field for a single weather
snapshot, not a time-marched fuel-moisture forecast. Units: temperature degrees C,
relative humidity %, elevation metres, moisture fractions.

References:
    Holden, Z.A.; Jolly, W.M. 2011. Modeling topographic influences on fuel
        moisture and fire danger in complex terrain to improve wildland fire
        management decision support. Forest Ecology and Management 262: 2033-2041.
    Rothermel, R.C. 1983. How to predict the spread and intensity of forest and
        range fires. USDA Forest Service GTR INT-143 (dead-fuel-moisture
        corrections for aspect, slope, shading, time of day).
    Nelson, R.M. 2000. Prediction of diurnal change in 10-h fuel stick moisture
        content. Canadian Journal of Forest Research 30: 1071-1087.
"""

from __future__ import annotations

from datetime import datetime

import numpy as np

from .atmosphere import equilibrium_moisture_content, latlon_grid

# Standard environmental lapse rate (degrees C per km) for optional elevation
# adjustment of a single reference air temperature to each cell's elevation.
DEFAULT_LAPSE_RATE = 6.5

# Default maximum near-fuel heating (degrees C) for a fully sun-exposed, open
# fuel relative to the 2 m air temperature. Sunlit surface fuels commonly run
# 10-20 C above air; 17 C is a representative open-sun maximum.
DEFAULT_DT_SUN = 17.0


def _es(temp_c):
    """Saturation vapour pressure (hPa) over water, Magnus form."""
    t = np.asarray(temp_c, dtype=float)
    return 6.112 * np.exp(17.67 * t / (t + 243.5))


def equation_of_time(day_of_year):
    """Equation of time (minutes) -- Spencer (1971) Fourier approximation.

    The signed difference between apparent (sundial) and mean (clock) solar time
    caused by Earth's orbital eccentricity and axial tilt; ranges about -14 to +16
    minutes over the year. Add it (with the longitude correction) to convert local
    standard clock time to true solar time.
    """
    b = 2.0 * np.pi * (day_of_year - 1) / 365.0
    return 229.18 * (0.000075 + 0.001868 * np.cos(b) - 0.032077 * np.sin(b)
                     - 0.014615 * np.cos(2 * b) - 0.040849 * np.sin(2 * b))


def solar_position(latitude, day_of_year, hour, *, longitude=None, timezone=None):
    """Sun zenith and azimuth (degrees) for a latitude, day of year and hour.

    ``latitude`` in degrees (N positive), ``day_of_year`` 1-365. By default ``hour``
    (0-24) is the **local solar hour** (noon = sun on the meridian). Pass both
    ``longitude`` (degrees, E positive) and ``timezone`` (hours from UTC, e.g. +1
    for CET) to instead treat ``hour`` as **local standard clock time**: it is then
    converted to true solar time with the longitude correction and the equation of
    time. Azimuth is measured clockwise from true north (a southern sun at solar
    noon in the northern hemisphere returns ~180). Returns ``(zenith_deg, azimuth_deg)``.
    """
    if longitude is not None and timezone is not None:
        # Local standard meridian at 15 deg/hour; 4 minutes of time per degree.
        correction_min = 4.0 * (longitude - 15.0 * timezone) + equation_of_time(
            day_of_year)
        hour = hour + correction_min / 60.0

    lat = np.radians(latitude)
    decl = np.radians(23.45) * np.sin(np.radians(360.0 * (284 + day_of_year) / 365.0))
    hour_angle = np.radians(15.0 * (hour - 12.0))

    cos_zen = (np.sin(lat) * np.sin(decl)
               + np.cos(lat) * np.cos(decl) * np.cos(hour_angle))
    cos_zen = np.clip(cos_zen, -1.0, 1.0)
    zenith = np.degrees(np.arccos(cos_zen))

    # Azimuth clockwise from north (Sproul 2007 form).
    az = np.arctan2(
        -np.cos(decl) * np.sin(hour_angle),
        np.cos(lat) * np.sin(decl) - np.sin(lat) * np.cos(decl) * np.cos(hour_angle),
    )
    azimuth = np.degrees(az) % 360.0
    return float(zenith), float(azimuth)


def terrain_insolation_factor(ls, zenith_deg, azimuth_deg):
    """Per-cell direct-beam factor on each slope facet (the cosine of incidence).

    ``cos(incidence) = cos(slope) cos(zenith)
                       + sin(slope) sin(zenith) cos(sun_azimuth - aspect)``,
    clamped at 0 (a facet turned away from the sun, or the sun below the horizon,
    receives no direct beam). For flat ground this reduces to ``cos(zenith)``, so
    the factor already carries the diurnal/seasonal sun height. ``ls.aspect`` is
    the downslope azimuth (degrees the slope faces); if absent, cells are treated
    as flat. Returns an array in ``[0, 1]`` shaped like the landscape.
    """
    shape = ls.shape
    if zenith_deg >= 90.0:                      # sun at/below the horizon
        return np.zeros(shape, dtype=float)

    slope_ang = np.arctan(np.asarray(ls.slope_tangent, dtype=float))
    if ls.aspect is None:
        aspect = np.zeros(shape, dtype=float)
        slope_ang = np.zeros(shape, dtype=float)
    else:
        aspect = np.radians(np.asarray(ls.aspect, dtype=float))

    z = np.radians(zenith_deg)
    sun_az = np.radians(azimuth_deg)
    cos_inc = (np.cos(slope_ang) * np.cos(z)
               + np.sin(slope_ang) * np.sin(z) * np.cos(sun_az - aspect))
    return np.clip(cos_inc, 0.0, 1.0)


def canopy_transmission(ls, *, shade=1.0):
    """Per-cell fraction of insolation reaching surface fuels through the canopy.

    ``1 - shade * canopy_cover_fraction``: at ``shade=1`` a 100%-cover cell blocks
    all direct beam, an open cell passes it all. ``shade`` (0-1) tempers that for
    canopies that leak light. Cells with no ``canopy_cover`` band pass everything.
    ``canopy_cover`` is read as a percent (0-100). Returns an array in ``[0, 1]``.
    """
    if ls.canopy_cover is None:
        return np.ones(ls.shape, dtype=float)
    cover = np.clip(np.asarray(ls.canopy_cover, dtype=float) / 100.0, 0.0, 1.0)
    return np.clip(1.0 - float(shade) * cover, 0.0, 1.0)


# --- VPD-based dead fuel moisture (semi-mechanistic alternative to EMC) --------

# Nolan et al. (2016), Eq. 8: the inverse-exponential FMD form of Resco de Dios
# et al. (2015) fitted across SE-Australian fuels. FM(%) = a + b*exp(-c*VPD_kPa),
# r2 = 0.66. These are calibration defaults -- refit with your own
# (VPD, FMC) data for a different fuel/climate.
NOLAN_VPD_COEFFS = (7.86, 140.94, 3.73)


def vapour_pressure_deficit(temp_c, relative_humidity):
    """Vapour pressure deficit (kPa) from air temperature (C) and RH (%).

    ``VPD = es(T) * (1 - RH/100)`` with ``es`` the Magnus saturation vapour
    pressure. Scalars or arrays.
    """
    rh = np.clip(np.asarray(relative_humidity, dtype=float), 0.0, 100.0)
    es_kpa = _es(temp_c) / 10.0                       # hPa -> kPa
    return es_kpa * (1.0 - rh / 100.0)


def dead_fuel_moisture_vpd(temp_c, relative_humidity, *, coeffs=NOLAN_VPD_COEFFS):
    """Dead fine fuel moisture (%) from VPD -- semi-mechanistic VPD model.

    ``FMC = a + b*exp(-c*VPD_kPa)`` (Resco de Dios et al. 2015; coefficients from
    Nolan et al. 2016). The inverse-exponential captures the rapid drying of fine
    dead fuels as the air gets drier; it is calibrated for fire-weather conditions
    and is not meant for near-saturation (VPD -> 0). Scalars or arrays.

    Two calibration caveats with the default SE-Australian coefficients: the model
    **floors at ``a`` (~7.86%)** as VPD grows, so in very dry air it saturates and
    terrain/canopy conditioning gains little leverage (and ``a`` is higher than the
    NFDRS EMC reaches in drought -- refit ``a`` for Mediterranean / US fine fuels).
    Pass your own ``coeffs`` fitted to local ``(VPD, FMC)`` data for other fuels.

    References:
        Resco de Dios, V., et al. 2015. A semi-mechanistic model for predicting the
            moisture content of fine litter. Agricultural and Forest Meteorology
            203: 64-73.
        Nolan, R.H., et al. 2016. Predicting dead fine fuel moisture at regional
            scales using vapour pressure deficit from MODIS and gridded weather
            data. Remote Sensing of Environment 174: 100-108.
    """
    a, b, c = coeffs
    vpd = vapour_pressure_deficit(temp_c, relative_humidity)
    return a + b * np.exp(-c * vpd)


def sun_exposure(ls, *, latitude=None, day_of_year=None, hour=None,
                 longitude=None, timezone=None,
                 zenith_deg=None, azimuth_deg=None, canopy=True, shade=1.0):
    """Per-cell sun-exposure index ``S in [0, 1]`` (terrain beam x canopy openness).

    Give either the sun position directly (``zenith_deg`` / ``azimuth_deg``) or the
    geometry to compute it (``latitude`` / ``day_of_year`` / ``hour``). Set
    ``canopy=False`` to ignore canopy shading. ``S`` drives the near-fuel heating
    in :func:`condition_dead_fuel_moisture`.
    """
    if zenith_deg is None or azimuth_deg is None:
        if latitude is None or day_of_year is None or hour is None:
            raise ValueError(
                "give zenith_deg+azimuth_deg, or latitude+day_of_year+hour")
        zenith_deg, azimuth_deg = solar_position(
            latitude, day_of_year, hour, longitude=longitude, timezone=timezone)
    s = terrain_insolation_factor(ls, zenith_deg, azimuth_deg)
    if canopy:
        s = s * canopy_transmission(ls, shade=shade)
    return s


def condition_dead_fuel_moisture(
    ls,
    *,
    temperature,
    relative_humidity,
    latitude=None,
    day_of_year=None,
    hour=None,
    longitude=None,
    timezone=None,
    zenith_deg=None,
    azimuth_deg=None,
    insolation=None,
    dt_sun: float = DEFAULT_DT_SUN,
    canopy: bool = True,
    shade: float = 1.0,
    reference_elevation=None,
    lapse_rate: float = DEFAULT_LAPSE_RATE,
    model: str = "emc",
    vpd_coeffs=NOLAN_VPD_COEFFS,
    offsets=(0.0, 1.0, 2.0),
):
    """Condition per-cell dead fuel moisture over a landscape (FlamMap-style).

    ``temperature`` (degrees C) and ``relative_humidity`` (%) are the input air
    values -- scalars (uniform weather) or 2D arrays already regridded to ``ls``
    (gridded NWP / reanalysis). The sun exposure that drives the per-cell drying is
    taken from ``insolation`` (a precomputed ``S`` array from :func:`sun_exposure`)
    if given, else computed from the sun position (``zenith_deg`` / ``azimuth_deg``)
    or the geometry (``latitude`` / ``day_of_year`` / ``hour``; pass ``longitude`` +
    ``timezone`` to treat ``hour`` as local clock time).

    ``model`` selects the moisture submodel applied to the per-cell fuel
    microclimate: ``"emc"`` (default) the NFDRS equilibrium moisture content, or
    ``"vpd"`` the semi-mechanistic VPD model (:func:`dead_fuel_moisture_vpd`,
    coefficients ``vpd_coeffs``) -- which is the more accurate point estimator for
    fine dead fuels in recent comparisons but carries fitted constants.

    If ``reference_elevation`` is given and ``ls.elevation`` is present, the input
    temperature is lapse-adjusted to each cell's elevation (``lapse_rate`` degrees C
    per km, vapour pressure held constant) before the solar heating is added -- use
    this when ``temperature`` is a single station value rather than a gridded field.

    Returns ``{"m_1h", "m_10h", "m_100h"}`` of per-cell moisture **fraction**
    arrays, ready to splat into :func:`pyflam.spread_field`. ``offsets`` (% added to
    the base moisture for the slower 10-/100-h classes) mirrors
    :func:`pyflam.atmosphere.dead_fuel_moisture`.
    """
    if model not in ("emc", "vpd"):
        raise ValueError("model must be 'emc' or 'vpd'")
    shape = ls.shape
    t_air = np.broadcast_to(np.asarray(temperature, dtype=float), shape).astype(float)
    rh_in = np.broadcast_to(
        np.asarray(relative_humidity, dtype=float), shape).astype(float)

    # Near-surface water-vapour pressure (held fixed through the adjustments below).
    e_vap = np.clip(rh_in, 0.0, 100.0) / 100.0 * _es(t_air)

    if reference_elevation is not None and ls.elevation is not None:
        dz = np.asarray(ls.elevation, dtype=float) - float(reference_elevation)
        t_air = t_air - lapse_rate * dz / 1000.0

    if insolation is None:
        s = sun_exposure(ls, latitude=latitude, day_of_year=day_of_year, hour=hour,
                         longitude=longitude, timezone=timezone,
                         zenith_deg=zenith_deg, azimuth_deg=azimuth_deg,
                         canopy=canopy, shade=shade)
    else:
        s = np.clip(np.broadcast_to(np.asarray(insolation, dtype=float), shape),
                    0.0, 1.0)

    t_fuel = t_air + float(dt_sun) * s
    rh_fuel = np.clip(100.0 * e_vap / _es(t_fuel), 0.0, 100.0)

    if model == "vpd":
        base = dead_fuel_moisture_vpd(t_fuel, rh_fuel, coeffs=vpd_coeffs)  # percent
    else:
        base = equilibrium_moisture_content(t_fuel, rh_fuel)              # percent
    base = np.maximum(base, 0.0)
    o1, o10, o100 = offsets
    return {
        "m_1h": np.ascontiguousarray((base + o1) / 100.0),
        "m_10h": np.ascontiguousarray((base + o10) / 100.0),
        "m_100h": np.ascontiguousarray((base + o100) / 100.0),
    }


def _as_datetime(t):
    """Coerce a datetime or an ISO-ish ``'YYYY-MM-DD HH:MM'`` string to datetime."""
    if isinstance(t, datetime):
        return t
    s = str(t).strip().replace("T", " ")
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d %H", "%Y-%m-%d"):
        try:
            return datetime.strptime(s, fmt)
        except ValueError:
            continue
    raise ValueError(
        f"could not parse time {t!r}; pass a datetime or 'YYYY-MM-DD HH:MM'")


def condition_from_weather(
    ls,
    *,
    time,
    atmosphere=None,
    temperature=None,
    relative_humidity=None,
    latitude=None,
    longitude=None,
    timezone=None,
    **conditioning_kwargs,
):
    """Condition per-cell dead fuel moisture from **meteorology** (or manual T/RH).

    This is the run-setup front door: instead of typing fixed fuel moistures, derive
    the initial dead fuel moisture for a given **date, time and location** from
    weather, then let terrain insolation + canopy shading
    (:func:`condition_dead_fuel_moisture`) spread it across the landscape. Air
    temperature and RH come from, in order of precedence:

    1. an ``atmosphere`` provider (any :class:`pyflam.atmosphere.AtmosphereProvider`
       -- ``ConstantAtmosphere``, ``GriddedAtmosphere``, or a live
       ``fetch_gfs`` / ``fetch_era5`` / ``open_atmosphere`` source), sampled at
       ``time``. If the landscape is geolocated (has a CRS) the provider is sampled
       **per cell** (``field_on``) so a gridded forecast varies across the domain
       and the terrain conditioning adds sub-grid detail on top; otherwise it is
       sampled at a single point (``state_at`` at ``latitude``/``longitude``);
    2. explicit ``temperature`` (deg C) and ``relative_humidity`` (%) -- the
       **manual fallback** for when no meteo data are available.

    ``time`` is a datetime or ``'YYYY-MM-DD HH:MM'`` string; its day-of-year and hour
    set the sun position. Provider timestamps are UTC, so when a provider is used
    and ``timezone`` is not given it defaults to ``0`` (UTC) with the site longitude,
    which is the correct way to turn a UTC time into true solar time. ``latitude`` /
    ``longitude`` are taken from the geolocated landscape when available, else must
    be supplied. Extra keywords (``model``, ``dt_sun``, ``canopy``, ``shade``,
    ``reference_elevation``, ``offsets`` ...) pass through to
    :func:`condition_dead_fuel_moisture`.

    Returns the ``{"m_1h", "m_10h", "m_100h"}`` per-cell moisture dict.
    """
    when = _as_datetime(time)

    if atmosphere is not None:
        ll = latlon_grid(ls)
        if ll is not None:                       # geolocated: per-cell macro field
            lat2d, lon2d = ll
            field = atmosphere.field_on(ls, when, latlon=ll)
            temperature = field.temperature
            relative_humidity = field.relative_humidity
            if latitude is None:
                latitude = float(np.mean(lat2d))
            if longitude is None:
                longitude = float(np.mean(lon2d))
        else:                                    # no CRS: sample one representative point
            if latitude is None:
                raise ValueError(
                    "landscape has no CRS: pass latitude= (and longitude=) so the "
                    "provider can be sampled and the sun positioned")
            state = atmosphere.state_at(latitude, longitude, when)
            temperature = state.temperature
            relative_humidity = state.relative_humidity
        if timezone is None and longitude is not None:
            timezone = 0.0                       # provider time is UTC

    if temperature is None or relative_humidity is None:
        raise ValueError(
            "provide an `atmosphere` provider, or temperature= and "
            "relative_humidity= (the manual fallback when no meteo data exist)")
    if latitude is None:
        raise ValueError("need latitude= for the sun position")

    hour = when.hour + when.minute / 60.0 + when.second / 3600.0
    return condition_dead_fuel_moisture(
        ls, temperature=temperature, relative_humidity=relative_humidity,
        latitude=latitude, day_of_year=when.timetuple().tm_yday, hour=hour,
        longitude=longitude, timezone=timezone, **conditioning_kwargs)
