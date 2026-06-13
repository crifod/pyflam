"""Rothermel (1972) surface fire spread model.

This is the scientific core of FlamMap's "Basic Fire Behavior" outputs. Given a
fuel bed, the dead/live fuel moistures, a midflame wind speed and a slope, it
returns rate of spread, reaction intensity, fireline intensity and flame length.

Everything is computed in Rothermel's native English units (see
:mod:`pyflam.units`). Use :func:`spread` for a high-level call that accepts a
:class:`~pyflam.fuel_models.FuelModel` and friendlier inputs.

References:
    Rothermel, R.C. 1972. A mathematical model for predicting fire spread in
        wildland fuels. USDA Forest Service Research Paper INT-115.
    Albini, F.A. 1976. Estimating wildfire behavior and effects. GTR INT-30.
    Andrews, P.L. 2018. The Rothermel surface fire spread model and associated
        developments: A comprehensive explanation. RMRS-GTR-371.

Equation numbers in comments refer to Andrews (2018), GTR-371.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .fuel_models import FuelModel, SIGMA_10H, SIGMA_100H

# Particle constants (Rothermel 1972), shared by all standard fuels.
PARTICLE_DENSITY = 32.0    # rho_p, oven-dry particle density, lb/ft^3
TOTAL_MINERAL = 0.0555     # S_T
EFFECTIVE_MINERAL = 0.010  # S_e

# Albini (1976) surface-area size class boundaries (1/ft), descending.
_SIZE_CLASS_BOUNDS = (1200.0, 192.0, 96.0, 48.0, 16.0)


def _size_class(sigma: float) -> int:
    for i, bound in enumerate(_SIZE_CLASS_BOUNDS):
        if sigma >= bound:
            return i
    return len(_SIZE_CLASS_BOUNDS)


@dataclass
class FireBehavior:
    """Surface fire behavior at a point (English units)."""

    rate_of_spread: float       # ft/min, in the heading direction
    reaction_intensity: float   # Btu/ft^2/min
    fireline_intensity: float   # Btu/ft/s (Byram)
    flame_length: float         # ft (Byram)
    heat_per_unit_area: float   # Btu/ft^2
    # Intermediate quantities, useful for validation / debugging.
    wind_factor: float          # phi_w
    slope_factor: float         # phi_s
    packing_ratio: float        # beta
    characteristic_sav: float   # sigma-bar, 1/ft


@dataclass
class SurfaceKernel:
    """Wind/slope-independent Rothermel terms for one fuel + moisture state.

    Apply a wind speed (ft/min) and slope (rise/run) — scalar or NumPy array —
    to get rate of spread and the rest of the fire behavior. Arrays let a whole
    landscape, with a per-cell wind field, be evaluated at once.
    """

    r0: float                  # rate of spread with all factors = 1 (ft/min)
    c: float                   # wind factor coefficient C (Eq. 47)
    b: float                   # wind factor exponent B
    e: float                   # wind factor exponent E
    beta_ratio: float          # beta / beta_op
    beta: float                # packing ratio
    reaction_intensity: float  # Btu/ft^2/min
    heat_per_unit_area: float  # Btu/ft^2  (I_R * residence time)
    characteristic_sav: float  # sigma-bar, 1/ft

    def wind_factor(self, wind_midflame):
        u = np.asarray(wind_midflame, dtype=float)
        phi = self.c * np.power(np.where(u > 0.0, u, 0.0), self.b) \
            * self.beta_ratio ** -self.e
        return np.where(u > 0.0, phi, 0.0)

    def slope_factor(self, slope):
        t = np.asarray(slope, dtype=float)
        phi = 5.275 * self.beta ** -0.3 * np.where(t > 0.0, t, 0.0) ** 2
        return np.where(t > 0.0, phi, 0.0)

    def rate_of_spread(self, wind_midflame=0.0, slope=0.0):
        return self.r0 * (1.0 + self.wind_factor(wind_midflame)
                          + self.slope_factor(slope))

    def behavior(self, *, wind_midflame=0.0, slope=0.0) -> FireBehavior:
        """Full scalar fire behavior for one (wind, slope) pair."""
        phi_w = float(self.wind_factor(wind_midflame))
        phi_s = float(self.slope_factor(slope))
        ros = self.r0 * (1.0 + phi_w + phi_s)
        fli = self.heat_per_unit_area * ros / 60.0
        return FireBehavior(
            rate_of_spread=ros,
            reaction_intensity=self.reaction_intensity,
            fireline_intensity=fli,
            flame_length=0.45 * fli ** 0.46 if fli > 0.0 else 0.0,
            heat_per_unit_area=self.heat_per_unit_area,
            wind_factor=phi_w,
            slope_factor=phi_s,
            packing_ratio=self.beta,
            characteristic_sav=self.characteristic_sav,
        )


# Returned for nonburnable fuel models (NB1-NB9) — no fuel, no fire.
_ZERO_BEHAVIOR = FireBehavior(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
_ZERO_KERNEL = SurfaceKernel(0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0)

# Live herbaceous moisture bounds for curing (Scott & Burgan 2005): fully cured
# at/below 30%, fully green at/above 120%, linear in between.
_CURE_WET = 1.20
_CURE_DRY = 0.30


def _cured_fraction(m_live_herb: float) -> float:
    """Fraction of live herbaceous load transferred to dead (0..1)."""
    frac = (_CURE_WET - m_live_herb) / (_CURE_WET - _CURE_DRY)
    return min(max(frac, 0.0), 1.0)


def spread(
    fuel: FuelModel,
    *,
    m_1h: float,
    m_10h: float,
    m_100h: float,
    m_live_herb: float = 0.0,
    m_live_woody: float = 0.0,
    wind_midflame: float = 0.0,
    slope: float = 0.0,
    load_factor: float = 1.0,
) -> FireBehavior:
    """Compute surface fire behavior for a standard fuel model.

    Parameters
    ----------
    fuel:
        The fuel bed (:class:`~pyflam.fuel_models.FuelModel`).
    m_1h, m_10h, m_100h:
        Dead fuel moisture contents, as fractions (e.g. 0.06 for 6%).
    m_live_herb, m_live_woody:
        Live fuel moisture contents, as fractions.
    wind_midflame:
        Midflame wind speed in ft/min (see :mod:`pyflam.units` to convert).
    slope:
        Slope steepness as rise/run (tan of the slope angle), e.g. 0.30 for 30%.
    load_factor:
        Multiplier on all fuel loads (default 1.0). The standard 13/40 fuel
        models are known to *under*-represent real fuel loading (often by 20-30%);
        pass e.g. 1.3 to add 30%. Loading flows through the full Rothermel
        computation, so it raises reaction intensity, heat per unit area and
        fireline intensity (the energy that drives flame length and spotting), and
        shifts the packing ratio.
    """
    if not fuel.is_burnable:
        return _ZERO_BEHAVIOR
    kernel = surface_kernel(
        fuel, m_1h=m_1h, m_10h=m_10h, m_100h=m_100h,
        m_live_herb=m_live_herb, m_live_woody=m_live_woody,
        load_factor=load_factor,
    )
    return kernel.behavior(wind_midflame=wind_midflame, slope=slope)


def surface_kernel(
    fuel: FuelModel,
    *,
    m_1h: float,
    m_10h: float,
    m_100h: float,
    m_live_herb: float = 0.0,
    m_live_woody: float = 0.0,
    load_factor: float = 1.0,
) -> "SurfaceKernel":
    """Precompute the wind/slope-independent part of the Rothermel model.

    The reaction intensity, propagating-flux ratio, heat sink and the wind/slope
    coefficients depend only on the fuel and its moisture, not on wind speed or
    slope. Computing them once and then applying many (wind, slope) pairs is how
    :func:`pyflam.landscape.basic_fire_behavior` stays fast and how a
    spatially-varying wind field (from :mod:`pyflam.windsolver` or
    :mod:`pyflam.cfd`) is applied per cell.

    ``load_factor`` scales all fuel loads (see :func:`spread`) to correct the
    standard models' known low bias; it propagates through the whole kernel.
    """
    if not fuel.is_burnable:
        return _ZERO_KERNEL

    lf = load_factor
    # Dynamic herbaceous load transfer (Scott & Burgan 2005): in dynamic models,
    # a moisture-dependent fraction of the live herb load is "cured" and moves to
    # a dead herbaceous class (same fine SAV, at the 1-h dead moisture).
    cured = _cured_fraction(m_live_herb) if fuel.dynamic else 0.0
    dead_herb = fuel.load_live_herb * cured * lf
    live_herb = fuel.load_live_herb * (1.0 - cured) * lf

    # Assemble particle arrays. dead = [1h, 10h, 100h, cured herb];
    # live = [live herb, woody]. Drop any zero-load class so it contributes
    # nothing and never divides by 0 (a zero-SAV slot would blow up exp(-138/s)).
    dead_load, dead_sav, dead_moist = _pack(
        [fuel.load_1h * lf, fuel.load_10h * lf, fuel.load_100h * lf, dead_herb],
        [fuel.sav_1h, SIGMA_10H, SIGMA_100H, fuel.sav_live_herb],
        [m_1h, m_10h, m_100h, m_1h],
    )
    live_load, live_sav, live_moist = _pack(
        [live_herb, fuel.load_live_woody * lf],
        [fuel.sav_live_herb, fuel.sav_live_woody],
        [m_live_herb, m_live_woody],
    )

    return _kernel_from_arrays(
        dead_load, dead_sav, dead_moist,
        live_load, live_sav, live_moist,
        depth=fuel.depth, mx_dead=fuel.mx_dead,
        heat_dead=fuel.heat_dead, heat_live=fuel.heat_live,
    )


# Quantization steps for grouping per-cell kernel inputs (bins one kernel each).
_KERNEL_QUANTA = {
    "load_factor": 0.05, "m_1h": 0.005, "m_10h": 0.005, "m_100h": 0.005,
    "m_live_herb": 0.02, "m_live_woody": 0.02,
}


def load_factor_groups(load_factor, fuel_num, mask, *, quantum: float = 0.05):
    """Yield ``(submask, lf)`` for one fuel group (scalar / per-fuel dict / array).

    Back-compat helper; :func:`kernel_param_groups` is the general form.
    """
    if isinstance(load_factor, dict):
        yield mask, float(load_factor.get(int(fuel_num), 1.0))
    else:
        for sub, params in kernel_param_groups(mask, {"load_factor": load_factor}):
            yield sub, params["load_factor"]


def kernel_param_groups(mask, params, *, quanta=None):
    """Group cells of one fuel by their per-cell kernel inputs.

    ``params`` maps kernel argument names (``load_factor``, ``m_1h`` ...) to a
    scalar **or a 2D array** on the landscape grid. Cells in ``mask`` are grouped
    so each distinct (quantized) combination of the array-valued inputs is one
    group -- the surface kernel is nonlinear in load and moisture, so these must
    be applied *inside* the kernel, computed once per group. Yields
    ``(submask, resolved)`` where ``resolved`` is all params as scalars (the group
    mean for the array ones). All-scalar params yield a single group (the fast,
    unchanged path).
    """
    quanta = {**_KERNEL_QUANTA, **(quanta or {})}
    array_names = [k for k, v in params.items() if not np.isscalar(v)]
    if not array_names:
        yield mask, {k: float(v) for k, v in params.items()}
        return
    ridx, cidx = np.where(mask)
    keys = np.stack(
        [np.round(np.asarray(params[k], float)[ridx, cidx]
                  / quanta.get(k, 0.01)).astype(np.int64) for k in array_names],
        axis=1)
    uniq, inv = np.unique(keys, axis=0, return_inverse=True)
    inv = inv.ravel()
    for gi in range(len(uniq)):
        sel = inv == gi
        sub = np.zeros_like(mask)
        sub[ridx[sel], cidx[sel]] = True
        resolved = {}
        for k, v in params.items():
            resolved[k] = (float(v) if np.isscalar(v)
                           else float(np.asarray(v, float)[ridx[sel], cidx[sel]].mean()))
        yield sub, resolved


def _pack(load, sav, moist):
    """Keep only the fuel classes with positive load (avoids 0-SAV divisions)."""
    load = np.asarray(load, dtype=float)
    mask = load > 0.0
    return load[mask], np.asarray(sav, float)[mask], np.asarray(moist, float)[mask]


def _category_weights(load, sav):
    """Mean-surface-area weighting within a fuel category (Eq. 53-56)."""
    area = sav * load / PARTICLE_DENSITY  # A_ij
    total = area.sum()
    if total <= 0.0:
        return np.zeros_like(area), 0.0
    return area / total, total  # f_ij, A_j


def _net_load(load, sav, f):
    """Size-class (g-factor) weighted net fuel load for a category (Eq. 59-60)."""
    net = load * (1.0 - TOTAL_MINERAL)
    classes = np.array([_size_class(s) for s in sav])
    out = 0.0
    for k in np.unique(classes):
        mask = classes == k
        g = f[mask].sum()          # g_k = sum of f over particles in class k
        out += g * net[mask].sum()
    return out


def _kernel_from_arrays(
    dead_load, dead_sav, dead_moist,
    live_load, live_sav, live_moist,
    *, depth, mx_dead, heat_dead, heat_live,
) -> SurfaceKernel:
    has_live = live_load.sum() > 0.0

    f_dead, area_dead = _category_weights(dead_load, dead_sav)
    f_live, area_live = _category_weights(live_load, live_sav)
    area_total = area_dead + area_live
    f_cat_dead = area_dead / area_total
    f_cat_live = area_live / area_total

    # Characteristic SAV (Eq. 71-73).
    sigma_dead = float((f_dead * dead_sav).sum())
    sigma_live = float((f_live * live_sav).sum()) if has_live else 0.0
    sigma = f_cat_dead * sigma_dead + f_cat_live * sigma_live

    # Bulk density and packing ratio (Eq. 31, 74).
    total_load = dead_load.sum() + live_load.sum()
    rho_b = total_load / depth
    beta = rho_b / PARTICLE_DENSITY
    beta_op = 3.348 * sigma ** -0.8189          # optimum packing ratio (Eq. 69)

    # --- Moisture damping --------------------------------------------------
    # Category mean moistures (area-weighted, Eq. 66).
    mf_dead = float((f_dead * dead_moist).sum())
    mf_live = float((f_live * live_moist).sum()) if has_live else 0.0

    # Live fuel moisture of extinction (Albini 1976, Eq. 88).
    mx_live = mx_dead
    if has_live:
        fine_dead_w = dead_load * np.exp(-138.0 / dead_sav)
        fine_live_w = live_load * np.exp(-500.0 / live_sav)
        w_ratio = fine_dead_w.sum() / fine_live_w.sum()
        fine_dead_moist = float((dead_moist * fine_dead_w).sum() / fine_dead_w.sum())
        mx_live = max(
            2.9 * w_ratio * (1.0 - fine_dead_moist / mx_dead) - 0.226,
            mx_dead,
        )

    eta_m_dead = _moisture_damping(mf_dead, mx_dead)
    eta_m_live = _moisture_damping(mf_live, mx_live) if has_live else 0.0

    # Mineral damping (Eq. 62), same for both categories here.
    eta_s = min(0.174 * EFFECTIVE_MINERAL ** -0.19, 1.0)

    # --- Reaction intensity (Eq. 58, 75) -----------------------------------
    wn_dead = _net_load(dead_load, dead_sav, f_dead)
    wn_live = _net_load(live_load, live_sav, f_live) if has_live else 0.0

    # Optimum reaction velocity (Eq. 68, 70).
    a = 133.0 * sigma ** -0.7913
    gamma_max = sigma ** 1.5 / (495.0 + 0.0594 * sigma ** 1.5)
    ratio = beta / beta_op
    gamma = gamma_max * ratio ** a * math.exp(a * (1.0 - ratio))

    reaction_intensity = gamma * (
        wn_dead * heat_dead * eta_m_dead * eta_s
        + wn_live * heat_live * eta_m_live * eta_s
    )

    # --- Propagating flux ratio (Eq. 76) -----------------------------------
    xi = math.exp((0.792 + 0.681 * sigma ** 0.5) * (beta + 0.1)) / (
        192.0 + 0.2595 * sigma
    )

    # --- Wind factor coefficients (Eq. 47-49); applied later per (cell) -----
    c = 7.47 * math.exp(-0.133 * sigma ** 0.55)
    b = 0.02526 * sigma ** 0.54
    e = 0.715 * math.exp(-3.59e-4 * sigma)

    # --- Heat sink (Eq. 77, 14, 12) ----------------------------------------
    eps_dead = np.exp(-138.0 / dead_sav)
    eps_live = np.exp(-138.0 / live_sav) if has_live else np.zeros_like(live_sav)
    qig_dead = 250.0 + 1116.0 * dead_moist
    qig_live = 250.0 + 1116.0 * live_moist

    heat_sink_dead = float((f_dead * eps_dead * qig_dead).sum())
    heat_sink_live = float((f_live * eps_live * qig_live).sum()) if has_live else 0.0
    heat_sink = rho_b * (f_cat_dead * heat_sink_dead + f_cat_live * heat_sink_live)

    # --- Base rate of spread with all factors = 1 (Eq. 52) -----------------
    r0 = reaction_intensity * xi / heat_sink if heat_sink > 0.0 else 0.0

    residence_time = 384.0 / sigma if sigma > 0.0 else 0.0  # min (Anderson 1969)
    heat_per_area = reaction_intensity * residence_time     # Btu/ft^2

    return SurfaceKernel(
        r0=r0, c=c, b=b, e=e,
        beta_ratio=ratio, beta=beta,
        reaction_intensity=reaction_intensity,
        heat_per_unit_area=heat_per_area,
        characteristic_sav=sigma,
    )


def _moisture_damping(mf: float, mx: float) -> float:
    """Moisture damping coefficient eta_M (Eq. 64)."""
    if mx <= 0.0:
        return 0.0
    r = min(mf / mx, 1.0)
    return 1.0 - 2.59 * r + 5.11 * r ** 2 - 3.52 * r ** 3
