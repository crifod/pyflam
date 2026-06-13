"""Wind reduction / adjustment factor: 20-ft wind -> midflame wind.

The Rothermel surface model wants the *midflame* wind, but weather and FlamMap's
landscape winds are at 20 ft (6.1 m). The wind adjustment factor (WAF, a.k.a.
wind reduction factor) bridges them. It is not a single number: it depends on how
deep the fuel bed is and whether a tree canopy shelters the surface.

This module implements the standard models FlamMap/BehavePlus use (Albini &
Baughman 1979, as organized by Andrews 2012):

* **Unsheltered** (no/!sparse canopy) — depends only on fuel bed depth, so it
  varies *per fuel model*: shallow litter reduces the wind a lot (WAF ~ 0.25),
  a tall shrub bed much less (WAF ~ 0.55).
* **Sheltered** (under a canopy) — depends on canopy height and cover; used when
  canopy cover clears a threshold and canopy height is known.

Using the right per-cell WAF (instead of a flat 0.4) is roadmap step 3's
remaining piece and the main lever for matching FlamMap's ROS on forested or
mixed landscapes.

Everything is in English units (fuel depth and canopy height in ft); WAF is a
dimensionless multiplier in (0, 1].

References:
    Albini, F.A.; Baughman, R.G. 1979. Estimating windspeeds for predicting
        wildland fire behavior. USDA Forest Service Research Paper INT-221.
    Andrews, P.L. 2012. Modeling wind adjustment factor and midflame wind speed
        for Rothermel's surface fire spread model. USDA FS RMRS-GTR-266.
"""

from __future__ import annotations

import numpy as np

from . import fuel_models

# FlamMap uses the sheltered WAF only where the canopy actually shelters the
# surface: canopy cover at/above this fraction and a known (positive) height.
DEFAULT_COVER_THRESHOLD = 0.05


def unsheltered_waf(fuel_depth_ft):
    """Unsheltered wind adjustment factor from fuel bed depth (ft).

    ``WAF = 1.83 / ln((20 + 0.36 H) / (0.13 H))`` (Albini & Baughman 1979). For
    non-positive depth (e.g. nonburnable) returns 1.0 (no reduction), but such
    cells don't spread anyway.
    """
    h = np.asarray(fuel_depth_ft, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        waf = 1.83 / np.log((20.0 + 0.36 * h) / (0.13 * h))
    return np.clip(np.where(h > 0.0, waf, 1.0), 0.0, 1.0)


def sheltered_waf(canopy_height_ft, canopy_cover_fraction):
    """Sheltered (under-canopy) wind adjustment factor.

    ``WAF = 0.555 / (sqrt(f H) * ln((20 + 0.36 H)/(0.13 H)))`` (Albini &
    Baughman 1979), with ``H`` the canopy height (ft) and ``f`` the canopy cover
    fraction. Returns 1.0 where height is non-positive (formula undefined).
    """
    h = np.asarray(canopy_height_ft, dtype=float)
    f = np.clip(np.asarray(canopy_cover_fraction, dtype=float), 0.0, 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        waf = 0.555 / (np.sqrt(f * h) * np.log((20.0 + 0.36 * h) / (0.13 * h)))
    return np.clip(np.where(h > 0.0, waf, 1.0), 0.0, 1.0)


def wind_adjustment_factor(
    fuel_depth_ft,
    *,
    canopy_height_ft=0.0,
    canopy_cover_fraction=0.0,
    cover_threshold: float = DEFAULT_COVER_THRESHOLD,
):
    """Per-element WAF, choosing sheltered vs unsheltered like FlamMap.

    Sheltered where ``canopy_cover_fraction >= cover_threshold`` and
    ``canopy_height_ft > 0``; unsheltered (fuel-depth based) otherwise. All
    inputs broadcast together; returns an array (or scalar) of WAF in (0, 1].
    """
    depth = np.asarray(fuel_depth_ft, dtype=float)
    ch = np.broadcast_to(np.asarray(canopy_height_ft, dtype=float), depth.shape)
    cc = np.broadcast_to(np.asarray(canopy_cover_fraction, dtype=float), depth.shape)

    unshel = unsheltered_waf(depth)
    use_sheltered = (cc >= cover_threshold) & (ch > 0.0)
    if not np.any(use_sheltered):
        return unshel
    shel = sheltered_waf(ch, cc)
    return np.where(use_sheltered, shel, unshel)


def fuel_depth_field(ls) -> np.ndarray:
    """Per-cell fuel bed depth (ft) from the landscape's fuel-model grid."""
    fuel = np.asarray(ls.fuel_model)
    depth = np.zeros(ls.shape, dtype=float)
    for num in np.unique(fuel):
        try:
            fm = fuel_models.get(int(num))
        except KeyError:
            continue
        depth[fuel == num] = fm.depth
    return depth


def waf_field(
    ls,
    *,
    canopy_height_ft=None,
    cover_threshold: float = DEFAULT_COVER_THRESHOLD,
) -> np.ndarray:
    """Per-cell wind adjustment factor over a landscape.

    Fuel depth comes from the fuel-model grid; canopy cover from
    ``ls.canopy_cover`` (percent -> fraction). Canopy height (ft) may be passed
    explicitly (e.g. already converted from an LCP/LANDFIRE band); if omitted and
    the landscape has no usable height, every cell uses the unsheltered WAF --
    which is the per-fuel-depth variation that matters on canopy-free data.
    """
    depth = fuel_depth_field(ls)
    if canopy_height_ft is None:
        ch = np.zeros(ls.shape, dtype=float)
    else:
        ch = np.broadcast_to(np.asarray(canopy_height_ft, dtype=float), ls.shape)
    if ls.canopy_cover is not None:
        cc = np.clip(np.asarray(ls.canopy_cover, dtype=float) / 100.0, 0.0, 1.0)
    else:
        cc = np.zeros(ls.shape, dtype=float)
    return wind_adjustment_factor(
        depth, canopy_height_ft=ch, canopy_cover_fraction=cc,
        cover_threshold=cover_threshold,
    )


def midflame_field(ls, wind_20ft, *, canopy_height_ft=None,
                   cover_threshold: float = DEFAULT_COVER_THRESHOLD) -> np.ndarray:
    """Midflame wind field = 20-ft wind * per-cell WAF (same units as input)."""
    waf = waf_field(ls, canopy_height_ft=canopy_height_ft,
                    cover_threshold=cover_threshold)
    return np.asarray(wind_20ft, dtype=float) * waf
