"""Unit conversions.

The Rothermel model is formulated in English units (Rothermel 1972), so all
internal computation in :mod:`pyflam.rothermel` uses them. These helpers convert
to/from the SI units that most users (and LANDFIRE data) work in.

Internal (English) units:
    fuel load           lb/ft^2
    SAV ratio (sigma)   1/ft
    fuel bed depth      ft
    heat content        Btu/lb
    rate of spread      ft/min
    reaction intensity  Btu/ft^2/min
    fireline intensity  Btu/ft/s
    wind speed          ft/min  (at midflame height)
"""

from __future__ import annotations

# Load -------------------------------------------------------------------------
TONS_ACRE_TO_LB_FT2 = 2000.0 / 43560.0  # 1 ton/acre -> lb/ft^2


def tons_per_acre_to_lb_per_ft2(load: float) -> float:
    return load * TONS_ACRE_TO_LB_FT2


def kg_per_m2_to_lb_per_ft2(load: float) -> float:
    return load * 0.204816


# Length -----------------------------------------------------------------------
def ft_to_m(x: float) -> float:
    return x * 0.3048


def m_to_ft(x: float) -> float:
    return x / 0.3048


def per_ft_to_per_m(sigma: float) -> float:
    """Surface-area-to-volume ratio: 1/ft -> 1/m."""
    return sigma / 0.3048


def per_m_to_per_ft(sigma: float) -> float:
    return sigma * 0.3048


# Speed ------------------------------------------------------------------------
def ft_per_min_to_m_per_min(r: float) -> float:
    return r * 0.3048


def m_per_min_to_ft_per_min(r: float) -> float:
    return r / 0.3048


def ft_per_min_to_m_per_s(r: float) -> float:
    return r * 0.3048 / 60.0


def mph_to_ft_per_min(u: float) -> float:
    return u * 88.0


def m_per_s_to_ft_per_min(u: float) -> float:
    return u * 60.0 / 0.3048


# Temperature -----------------------------------------------------------------
def kelvin_to_celsius(t: float) -> float:
    return t - 273.15


def celsius_to_fahrenheit(t: float) -> float:
    return t * 9.0 / 5.0 + 32.0


def ft_per_min_to_mph(u: float) -> float:
    return u / 88.0


def chains_per_hour_to_ft_per_min(r: float) -> float:
    """FlamMap's default ROS unit: chains/hour -> ft/min (1 chain = 66 ft)."""
    return r * 66.0 / 60.0


def ft_per_min_to_chains_per_hour(r: float) -> float:
    return r * 60.0 / 66.0


def km_per_h_to_ft_per_min(u: float) -> float:
    return u * 1000.0 / 60.0 / 0.3048


# Intensity --------------------------------------------------------------------
def btu_per_ft_s_to_kw_per_m(i: float) -> float:
    """Byram fireline intensity: Btu/ft/s -> kW/m."""
    return i * 3.46414


def kw_per_m_to_btu_per_ft_s(i: float) -> float:
    """Byram fireline intensity: kW/m -> Btu/ft/s."""
    return i / 3.46414
