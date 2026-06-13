"""Crown fire initiation and spread (roadmap step 4).

Surface fire is only half the story in forested fuels: once the surface
fireline intensity is high enough, fire climbs into the canopy. This module
links pyflam's surface model to the two canonical crown-fire models:

* **Van Wagner (1977)** — crown-fire *initiation*: the critical surface
  fireline intensity ``I_0`` needed to ignite the canopy (a function of canopy
  base height and foliar moisture), and the critical crown spread rate needed to
  *sustain* an active (solid) crown fire (a function of canopy bulk density).
* **Rothermel (1991)** — active crown-fire *spread rate*, the well-known
  ``R_active = 3.34 * R_10`` correlation (fuel model 10, wind reduction 0.4).
* **Scott & Reinhardt (2001)** — ties the two together: classifies the fire as
  *surface / passive (torching) / active*, computes the **crown fraction
  burned** (CFB) and a continuous final spread rate, and defines the
  **torching** and **crowning indices** (the 20-ft wind speeds at which a fire
  transitions to passive and active crowning).

Unit convention
---------------
Van Wagner and Scott & Reinhardt are formulated in SI, and canopy data
(LANDFIRE / LCP) is metric, so **this module works in SI**:

    canopy base height (CBH)   m
    canopy bulk density (CBD)  kg/m^3
    foliar moisture (FMC)      percent
    spread rate                m/min
    fireline intensity         kW/m

The surface model is English-internal, so :func:`crown_fire_behavior` takes the
surface :class:`~pyflam.rothermel.FireBehavior` directly and converts at the
boundary (see :mod:`pyflam.units`). Outputs are SI; convert back with
``pyflam.units`` if you want feet.

References:
    Van Wagner, C.E. 1977. Conditions for the start and spread of crown fire.
        Canadian Journal of Forest Research 7: 23-34.
    Rothermel, R.C. 1991. Predicting behavior and size of crown fires in the
        Northern Rocky Mountains. USDA Forest Service Research Paper INT-438.
    Scott, J.H.; Reinhardt, E.D. 2001. Assessing crown fire potential by linking
        models of surface and crown fire behavior. USDA Forest Service RMRS-RP-29.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from . import fuel_models
from .rothermel import FireBehavior, surface_kernel
from .units import (
    btu_per_ft_s_to_kw_per_m,
    ft_per_min_to_m_per_min,
    kw_per_m_to_btu_per_ft_s,
    m_per_min_to_ft_per_min,
)

# Van Wagner (1977) critical mass flow rate for a solid (active) crown fire,
# kg/m^2/min. R'_active = CRITICAL_MASS_FLOW / CBD.
CRITICAL_MASS_FLOW = 3.0

# Rothermel (1991) active crown-fire spread correlation: R_active = 3.34 * R_10,
# where R_10 is fuel model 10 run with a 0.4 wind-reduction factor.
_ROTHERMEL_CROWN_FACTOR = 3.34
_CROWN_WIND_REDUCTION = 0.4

# Scott & Reinhardt (2001): the crown fraction burned curve is scaled so that
# CFB = 0.9 at the active-crowning spread rate R'_active.
_CFB_AT_ACTIVE = 0.9


@dataclass
class CrownFireBehavior:
    """Crown-fire behavior at a point (SI units).

    ``fire_type`` is one of ``"surface"`` (no crown involvement),
    ``"passive"`` (torching / intermittent crowning) or ``"active"`` (solid
    crown fire). ``crown_fraction_burned`` is 0 for a surface fire and rises
    toward 1 as the fire becomes fully active.
    """

    fire_type: str
    rate_of_spread: float            # m/min, final (Scott & Reinhardt blend)
    surface_rate_of_spread: float    # m/min
    active_rate_of_spread: float     # m/min (Rothermel 1991)
    crown_fraction_burned: float     # 0..1
    fireline_intensity: float        # kW/m, total (surface + crown)
    surface_fireline_intensity: float          # kW/m
    critical_fireline_intensity: float          # kW/m, Van Wagner I_0
    critical_active_rate_of_spread: float       # m/min, Van Wagner R'_active
    initiates: bool                  # surface intensity >= I_0


def critical_fireline_intensity(canopy_base_height: float,
                                foliar_moisture: float) -> float:
    """Van Wagner (1977) critical surface fireline intensity ``I_0`` (kW/m).

    ``canopy_base_height`` in m, ``foliar_moisture`` in percent. Below this
    surface intensity the canopy will not ignite.
    """
    if canopy_base_height <= 0.0:
        return 0.0
    return (0.01 * canopy_base_height * (460.0 + 25.9 * foliar_moisture)) ** 1.5


def critical_active_ros(canopy_bulk_density: float) -> float:
    """Van Wagner (1977) critical spread rate for active crowning (m/min).

    ``R'_active = 3.0 / CBD``. Below this crown spread rate the canopy can't
    carry a solid (active) crown fire, so an initiated fire stays passive.
    ``canopy_bulk_density`` in kg/m^3; returns ``inf`` for zero CBD (no canopy).
    """
    if canopy_bulk_density <= 0.0:
        return math.inf
    return CRITICAL_MASS_FLOW / canopy_bulk_density


def active_crown_ros(
    wind_20ft_ft_per_min: float,
    *,
    m_1h: float,
    m_10h: float,
    m_100h: float,
    m_live_herb: float = 0.0,
    m_live_woody: float = 0.0,
) -> float:
    """Rothermel (1991) active crown-fire spread rate (m/min).

    ``R_active = 3.34 * R_10``, where ``R_10`` is the surface spread rate of
    fuel model 10 driven by the midflame wind ``0.4 * U_20ft``. ``wind_20ft`` is
    the 20-ft wind speed in ft/min (see :mod:`pyflam.units`).
    """
    fm10 = fuel_models.get(10)
    kernel = surface_kernel(
        fm10, m_1h=m_1h, m_10h=m_10h, m_100h=m_100h,
        m_live_herb=m_live_herb, m_live_woody=m_live_woody,
    )
    midflame = _CROWN_WIND_REDUCTION * float(wind_20ft_ft_per_min)
    r10_ft_min = float(kernel.rate_of_spread(midflame, 0.0))
    return _ROTHERMEL_CROWN_FACTOR * ft_per_min_to_m_per_min(r10_ft_min)


def crown_fraction_burned(surface_ros: float, ros_initiation: float,
                          active_ros_crit: float) -> float:
    """Scott & Reinhardt (2001) crown fraction burned (0..1).

    A smooth transition between the surface spread rate at initiation
    (``ros_initiation``) and the active-crowning spread rate
    (``active_ros_crit``); all rates in the same units.
    """
    if surface_ros <= ros_initiation or not math.isfinite(active_ros_crit):
        return 0.0
    span = active_ros_crit - ros_initiation
    if span <= 0.0:
        return 1.0
    a_c = -math.log(1.0 - _CFB_AT_ACTIVE) / span
    return float(min(1.0, 1.0 - math.exp(-a_c * (surface_ros - ros_initiation))))


def crown_fire_behavior(
    surface: FireBehavior,
    *,
    canopy_base_height: float,
    canopy_bulk_density: float,
    foliar_moisture: float,
    wind_20ft_ft_per_min: float,
    canopy_fuel_load: float = 0.0,
    heat_content: float = 18000.0,
    m_1h: float,
    m_10h: float,
    m_100h: float,
    m_live_herb: float = 0.0,
    m_live_woody: float = 0.0,
) -> CrownFireBehavior:
    """Classify and quantify crown fire at one point (Scott & Reinhardt 2001).

    Parameters
    ----------
    surface:
        Surface fire behavior from :func:`pyflam.spread` (English units).
    canopy_base_height, canopy_bulk_density, foliar_moisture:
        Canopy base height (m), bulk density (kg/m^3) and foliar moisture (%).
    wind_20ft_ft_per_min:
        20-ft wind speed (ft/min) driving the Rothermel (1991) active ROS.
    canopy_fuel_load:
        Available canopy fuel load (kg/m^2). Used only for the crown
        contribution to total fireline intensity; 0 omits it.
    heat_content:
        Heat of combustion of canopy fuel (kJ/kg), default 18 000.
    m_1h, m_10h, m_100h, m_live_herb, m_live_woody:
        Dead/live fuel moistures (fractions), for the fuel-model-10 active ROS.
    """
    r_surface = ft_per_min_to_m_per_min(surface.rate_of_spread)
    i_surface = btu_per_ft_s_to_kw_per_m(surface.fireline_intensity)

    i_crit = critical_fireline_intensity(canopy_base_height, foliar_moisture)
    r_active_crit = critical_active_ros(canopy_bulk_density)
    r_active = active_crown_ros(
        wind_20ft_ft_per_min, m_1h=m_1h, m_10h=m_10h, m_100h=m_100h,
        m_live_herb=m_live_herb, m_live_woody=m_live_woody,
    )

    initiates = i_surface >= i_crit and i_crit > 0.0

    # Surface spread rate at which the surface intensity equals I_0. Intensity is
    # linear in ROS (Byram), so R_init scales I_0 / I_surface.
    if i_surface > 0.0:
        r_init = r_surface * i_crit / i_surface
    else:
        r_init = math.inf

    if not initiates:
        cfb = 0.0
        fire_type = "surface"
        r_final = r_surface
    else:
        cfb = crown_fraction_burned(r_surface, r_init, r_active_crit)
        r_final = r_surface + cfb * (r_active - r_surface)
        fire_type = "active" if r_active >= r_active_crit else "passive"

    # Total fireline intensity: surface plus the canopy fuel actually consumed
    # (CFB * canopy load), as Byram intensity I = H * w * R (SI: kW/m).
    crown_load_consumed = cfb * canopy_fuel_load            # kg/m^2
    i_crown = heat_content * crown_load_consumed * (r_final / 60.0)  # kW/m
    i_total = i_surface + i_crown

    return CrownFireBehavior(
        fire_type=fire_type,
        rate_of_spread=r_final,
        surface_rate_of_spread=r_surface,
        active_rate_of_spread=r_active,
        crown_fraction_burned=cfb,
        fireline_intensity=i_total,
        surface_fireline_intensity=i_surface,
        critical_fireline_intensity=i_crit,
        critical_active_rate_of_spread=r_active_crit,
        initiates=initiates,
    )


def torching_index(
    fuel: fuel_models.FuelModel,
    *,
    canopy_base_height: float,
    foliar_moisture: float,
    wind_reduction_factor: float = 0.4,
    slope: float = 0.0,
    m_1h: float,
    m_10h: float,
    m_100h: float,
    m_live_herb: float = 0.0,
    m_live_woody: float = 0.0,
    max_wind_mph: float = 200.0,
) -> float:
    """Torching index: 20-ft wind speed (mph) at which crowning *initiates*.

    The wind at which this fuel's surface fireline intensity first reaches the
    Van Wagner critical intensity ``I_0``. Returns ``inf`` if it can't be
    reached below ``max_wind_mph`` (e.g. nonburnable or very tall canopy).
    """
    i_crit = critical_fireline_intensity(canopy_base_height, foliar_moisture)
    if i_crit <= 0.0 or not fuel.is_burnable:
        return math.inf
    kernel = surface_kernel(
        fuel, m_1h=m_1h, m_10h=m_10h, m_100h=m_100h,
        m_live_herb=m_live_herb, m_live_woody=m_live_woody,
    )

    def surface_intensity_kw(u_mph: float) -> float:
        midflame = wind_reduction_factor * u_mph * 88.0  # ft/min
        fb = kernel.behavior(wind_midflame=midflame, slope=slope)
        return btu_per_ft_s_to_kw_per_m(fb.fireline_intensity)

    return _solve_wind(surface_intensity_kw, i_crit, max_wind_mph)


def crowning_index(
    canopy_bulk_density: float,
    *,
    slope: float = 0.0,
    m_1h: float,
    m_10h: float,
    m_100h: float,
    m_live_herb: float = 0.0,
    m_live_woody: float = 0.0,
    max_wind_mph: float = 200.0,
) -> float:
    """Crowning index: 20-ft wind speed (mph) at which *active* crowning starts.

    The wind at which the Rothermel (1991) active crown ROS first reaches the
    Van Wagner critical active spread rate ``R'_active = 3.0 / CBD``. Returns
    ``inf`` for zero canopy bulk density.
    """
    r_active_crit = critical_active_ros(canopy_bulk_density)
    if not math.isfinite(r_active_crit):
        return math.inf

    fm10 = fuel_models.get(10)
    kernel = surface_kernel(
        fm10, m_1h=m_1h, m_10h=m_10h, m_100h=m_100h,
        m_live_herb=m_live_herb, m_live_woody=m_live_woody,
    )

    def active_ros_m(u_mph: float) -> float:
        midflame = _CROWN_WIND_REDUCTION * u_mph * 88.0  # ft/min
        r10 = float(kernel.rate_of_spread(midflame, slope))
        return _ROTHERMEL_CROWN_FACTOR * ft_per_min_to_m_per_min(r10)

    return _solve_wind(active_ros_m, r_active_crit, max_wind_mph)


def _solve_wind(f, target: float, max_wind_mph: float) -> float:
    """Bisection for the wind (mph) where monotone ``f(wind) == target``."""
    if f(0.0) >= target:
        return 0.0
    if f(max_wind_mph) < target:
        return math.inf
    lo, hi = 0.0, max_wind_mph
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if f(mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def crown_fire_potential(
    ls,
    *,
    foliar_moisture: float,
    wind_20ft_ft_per_min,
    m_1h: float,
    m_10h: float,
    m_100h: float,
    m_live_herb: float = 0.0,
    m_live_woody: float = 0.0,
    wind_midflame=0.0,
    cbh_scale: float = 0.1,
    cbd_scale: float = 0.01,
    heat_content: float = 18000.0,
) -> dict[str, np.ndarray]:
    """Per-cell crown-fire potential over a whole landscape.

    Runs the surface model (uniform moisture, per-cell slope), then applies the
    crown-fire classification of :func:`crown_fire_behavior` cell-by-cell.
    Requires the landscape to carry ``canopy_base_height`` and
    ``canopy_bulk_density`` bands.

    Canopy bands are interpreted as raw LCP/LANDFIRE integers by default and
    rescaled to SI: ``cbh_scale=0.1`` (m*10 -> m) and ``cbd_scale=0.01``
    ((kg/m^3)*100 -> kg/m^3). Pass ``1.0`` if your arrays are already in SI.

    Returns arrays (same shape as the landscape):
      ``fire_type`` (0 surface, 1 passive, 2 active; nonburnable/no-canopy = 0),
      ``rate_of_spread`` (m/min), ``crown_fraction_burned`` (0..1) and
      ``fireline_intensity`` (kW/m, total).
    """
    from .rothermel import surface_kernel as _kernel  # local: avoid cycle noise

    if ls.canopy_base_height is None or ls.canopy_bulk_density is None:
        raise ValueError(
            "crown_fire_potential needs canopy_base_height and "
            "canopy_bulk_density bands on the landscape"
        )

    shape = ls.shape
    fire_type = np.zeros(shape, dtype=np.int8)
    ros = np.zeros(shape, dtype=float)
    cfb = np.zeros(shape, dtype=float)
    fli = np.zeros(shape, dtype=float)

    tan_slope = ls.slope_tangent
    fuel = np.asarray(ls.fuel_model)
    wind_mid = np.broadcast_to(np.asarray(wind_midflame, dtype=float), shape)
    wind20 = np.broadcast_to(np.asarray(wind_20ft_ft_per_min, dtype=float), shape)
    cbh = np.asarray(ls.canopy_base_height, dtype=float) * cbh_scale
    cbd = np.asarray(ls.canopy_bulk_density, dtype=float) * cbd_scale
    # Available canopy fuel load (kg/m^2) ~ CBD * canopy depth, when canopy
    # height is available (same length scale as CBH); else no crown intensity.
    if ls.canopy_height is not None:
        canopy_depth = np.maximum(
            np.asarray(ls.canopy_height, dtype=float) * cbh_scale - cbh, 0.0)
        canopy_load = cbd * canopy_depth
    else:
        canopy_load = np.zeros(shape, dtype=float)
    moist = dict(m_1h=m_1h, m_10h=m_10h, m_100h=m_100h,
                 m_live_herb=m_live_herb, m_live_woody=m_live_woody)

    # Fuel-model-10 active crown ROS depends only on wind + moisture, so compute
    # it once across the grid (vectorized) rather than per cell.
    fm10_kernel = _kernel(fuel_models.get(10), **moist)
    r_active = _ROTHERMEL_CROWN_FACTOR * ft_per_min_to_m_per_min(
        fm10_kernel.rate_of_spread(_CROWN_WIND_REDUCTION * wind20, 0.0)
    )
    r_active_crit = np.where(cbd > 0.0, CRITICAL_MASS_FLOW / np.maximum(cbd, 1e-9),
                             np.inf)
    i_crit = np.where(
        cbh > 0.0,
        (0.01 * cbh * (460.0 + 25.9 * foliar_moisture)) ** 1.5,
        0.0,
    )

    for num in np.unique(fuel):
        num = int(num)
        mask = fuel == num
        try:
            fm = fuel_models.get(num)
        except KeyError:
            continue
        if not fm.is_burnable:
            continue

        kernel = _kernel(fm, **moist)
        r_surf_ft = kernel.rate_of_spread(wind_mid[mask], tan_slope[mask])
        r_surf = ft_per_min_to_m_per_min(np.asarray(r_surf_ft, dtype=float))
        i_surf = btu_per_ft_s_to_kw_per_m(kernel.heat_per_unit_area
                                          * np.asarray(r_surf_ft) / 60.0)

        ic = i_crit[mask]
        rac = r_active[mask]
        racc = r_active_crit[mask]
        initiates = (i_surf >= ic) & (ic > 0.0)

        with np.errstate(divide="ignore", invalid="ignore"):
            r_init = np.where(i_surf > 0.0, r_surf * ic / np.maximum(i_surf, 1e-12),
                              np.inf)
        span = racc - r_init
        a_c = np.where(span > 0.0,
                       -math.log(1.0 - _CFB_AT_ACTIVE) / np.where(span > 0.0, span, 1.0),
                       0.0)
        cfb_cell = np.where(
            initiates & (r_surf > r_init) & np.isfinite(racc),
            np.clip(1.0 - np.exp(-a_c * (r_surf - r_init)), 0.0, 1.0),
            0.0,
        )
        r_final = np.where(initiates, r_surf + cfb_cell * (rac - r_surf), r_surf)
        ftype = np.where(~initiates, 0, np.where(rac >= racc, 2, 1)).astype(np.int8)

        i_crown = heat_content * (cfb_cell * canopy_load[mask]) * (r_final / 60.0)

        ros[mask] = r_final
        cfb[mask] = cfb_cell
        fire_type[mask] = ftype
        fli[mask] = i_surf + i_crown

    return {
        "fire_type": fire_type,
        "rate_of_spread": ros,
        "crown_fraction_burned": cfb,
        "fireline_intensity": fli,
    }
