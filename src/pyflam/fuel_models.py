"""Standard fire-behavior fuel models.

Implements the original 13 fuel models of Anderson (1982) / Albini (1976),
which is what FlamMap's "Basic Fire Behavior" outputs and BehavePlus use as the
baseline reference set. The newer Scott & Burgan (2005) 40-model set can be added
on top of this same :class:`FuelModel` container later.

References:
    Anderson, H.E. 1982. Aids to determining fuel models for estimating fire
        behavior. USDA Forest Service GTR INT-122.
    Albini, F.A. 1976. Estimating wildfire behavior and effects. USDA Forest
        Service GTR INT-30.

Loads here are stored in lb/ft^2 (converted from the published tons/acre).
SAV ratios are in 1/ft; depth in ft; moisture of extinction as a fraction.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .units import tons_per_acre_to_lb_per_ft2 as _t

# Fixed SAV ratios shared by all standard models (1/ft).
SIGMA_10H = 109.0
SIGMA_100H = 30.0

# Heat content shared by all standard models (Btu/lb).
HEAT_CONTENT = 8000.0


@dataclass(frozen=True)
class FuelModel:
    """A surface fuel bed description.

    Loads are oven-dry, in lb/ft^2. Surface-area-to-volume ratios in 1/ft.
    Depth in ft. ``mx_dead`` (dead fuel moisture of extinction) as a fraction.
    Live SAV ratios apply to the live herbaceous / woody loads when present.
    """

    number: int
    code: str
    name: str
    load_1h: float
    load_10h: float
    load_100h: float
    load_live_herb: float
    load_live_woody: float
    sav_1h: float
    sav_live_herb: float
    sav_live_woody: float
    depth: float
    mx_dead: float
    heat_dead: float = HEAT_CONTENT
    heat_live: float = HEAT_CONTENT
    dynamic: bool = False  # dynamic live herbaceous load transfer (S&B 2005)

    @property
    def is_burnable(self) -> bool:
        return (
            self.load_1h + self.load_10h + self.load_100h
            + self.load_live_herb + self.load_live_woody
        ) > 0.0


def _fm(number, code, name, *, t1, t10, t100, t_herb=0.0, t_woody=0.0,
        sav_1h, sav_herb=0.0, sav_woody=0.0, depth, mx_dead_pct):
    """Build a FuelModel from published tons/acre + percent values."""
    return FuelModel(
        number=number, code=code, name=name,
        load_1h=_t(t1), load_10h=_t(t10), load_100h=_t(t100),
        load_live_herb=_t(t_herb), load_live_woody=_t(t_woody),
        sav_1h=sav_1h, sav_live_herb=sav_herb, sav_live_woody=sav_woody,
        depth=depth, mx_dead=mx_dead_pct / 100.0,
    )


# The original 13 (Anderson 1982). Live loads in models 2,4,5,7,10 are treated
# as live woody here (static, no dead/live moisture transfer).
STANDARD_13: dict[int, FuelModel] = {fm.number: fm for fm in [
    _fm(1, "FM1", "Short grass",
        t1=0.74, t10=0.0, t100=0.0, sav_1h=3500, depth=1.0, mx_dead_pct=12),
    _fm(2, "FM2", "Timber grass and understory",
        t1=2.0, t10=1.0, t100=0.5, t_woody=0.5,
        sav_1h=3000, sav_woody=1500, depth=1.0, mx_dead_pct=15),
    _fm(3, "FM3", "Tall grass",
        t1=3.01, t10=0.0, t100=0.0, sav_1h=1500, depth=2.5, mx_dead_pct=25),
    _fm(4, "FM4", "Chaparral",
        t1=5.01, t10=4.01, t100=2.0, t_woody=5.01,
        sav_1h=2000, sav_woody=1500, depth=6.0, mx_dead_pct=20),
    _fm(5, "FM5", "Brush",
        t1=1.0, t10=0.5, t100=0.0, t_woody=2.0,
        sav_1h=2000, sav_woody=1500, depth=2.0, mx_dead_pct=20),
    _fm(6, "FM6", "Dormant brush, hardwood slash",
        t1=1.5, t10=2.5, t100=2.0, sav_1h=1750, depth=2.5, mx_dead_pct=25),
    _fm(7, "FM7", "Southern rough",
        t1=1.13, t10=1.87, t100=1.5, t_woody=0.37,
        sav_1h=1750, sav_woody=1550, depth=2.5, mx_dead_pct=40),
    _fm(8, "FM8", "Closed timber litter",
        t1=1.5, t10=1.0, t100=2.5, sav_1h=2000, depth=0.2, mx_dead_pct=30),
    _fm(9, "FM9", "Hardwood litter",
        t1=2.92, t10=0.41, t100=0.15, sav_1h=2500, depth=0.2, mx_dead_pct=25),
    _fm(10, "FM10", "Timber litter and understory",
        t1=3.01, t10=2.0, t100=5.01, t_woody=2.0,
        sav_1h=2000, sav_woody=1500, depth=1.0, mx_dead_pct=25),
    _fm(11, "FM11", "Light logging slash",
        t1=1.5, t10=4.51, t100=5.51, sav_1h=1500, depth=1.0, mx_dead_pct=15),
    _fm(12, "FM12", "Medium logging slash",
        t1=4.01, t10=14.03, t100=16.53, sav_1h=1500, depth=2.3, mx_dead_pct=20),
    _fm(13, "FM13", "Heavy logging slash",
        t1=7.01, t10=23.04, t100=28.05, sav_1h=1500, depth=3.0, mx_dead_pct=25),
]}


# --- Scott & Burgan (2005) 40 fuel models -------------------------------------
# Reference: Scott, J.H.; Burgan, R.E. 2005. Standard fire behavior fuel models:
#   a comprehensive set for use with Rothermel's surface fire spread model.
#   USDA Forest Service RMRS-GTR-153.
#
# These differ from the original 13 in two ways pyflam must model:
#   * dynamic models transfer part of the live herbaceous load to a cured (dead)
#     herbaceous class based on live herb moisture (see rothermel._cured_fraction);
#   * the set includes nonburnable (NB) models, which never spread.
# All 40 use 8000 Btu/lb heat content for dead and live.


def _sb(number, code, name, *, t1, t10, t100, t_herb=0.0, t_woody=0.0,
        dynamic, sav_1h, sav_herb=0.0, sav_woody=0.0, depth, mx_dead_pct):
    return FuelModel(
        number=number, code=code, name=name,
        load_1h=_t(t1), load_10h=_t(t10), load_100h=_t(t100),
        load_live_herb=_t(t_herb), load_live_woody=_t(t_woody),
        sav_1h=sav_1h, sav_live_herb=sav_herb, sav_live_woody=sav_woody,
        depth=depth, mx_dead=mx_dead_pct / 100.0, dynamic=dynamic,
    )


def _nb(number, code, name):
    """Nonburnable model: no fuel, never spreads."""
    return FuelModel(
        number=number, code=code, name=name,
        load_1h=0.0, load_10h=0.0, load_100h=0.0,
        load_live_herb=0.0, load_live_woody=0.0,
        sav_1h=1.0, sav_live_herb=0.0, sav_live_woody=0.0,
        depth=0.1, mx_dead=0.10,
    )


STANDARD_40: dict[int, FuelModel] = {fm.number: fm for fm in [
    # Nonburnable
    _nb(91, "NB1", "Urban/developed"),
    _nb(92, "NB2", "Snow/ice"),
    _nb(93, "NB3", "Agricultural"),
    _nb(98, "NB8", "Open water"),
    _nb(99, "NB9", "Bare ground"),
    # Grass (all dynamic)
    _sb(101, "GR1", "Short, sparse dry climate grass",
        t1=0.10, t10=0.0, t100=0.0, t_herb=0.30, dynamic=True,
        sav_1h=2200, sav_herb=2000, depth=0.4, mx_dead_pct=15),
    _sb(102, "GR2", "Low load dry climate grass",
        t1=0.10, t10=0.0, t100=0.0, t_herb=1.0, dynamic=True,
        sav_1h=2000, sav_herb=1800, depth=1.0, mx_dead_pct=15),
    _sb(103, "GR3", "Low load very coarse humid climate grass",
        t1=0.10, t10=0.40, t100=0.0, t_herb=1.50, dynamic=True,
        sav_1h=1500, sav_herb=1300, depth=2.0, mx_dead_pct=30),
    _sb(104, "GR4", "Moderate load dry climate grass",
        t1=0.25, t10=0.0, t100=0.0, t_herb=1.90, dynamic=True,
        sav_1h=2000, sav_herb=1800, depth=2.0, mx_dead_pct=15),
    _sb(105, "GR5", "Low load humid climate grass",
        t1=0.40, t10=0.0, t100=0.0, t_herb=2.50, dynamic=True,
        sav_1h=1800, sav_herb=1600, depth=1.5, mx_dead_pct=40),
    _sb(106, "GR6", "Moderate load humid climate grass",
        t1=0.10, t10=0.0, t100=0.0, t_herb=3.40, dynamic=True,
        sav_1h=2200, sav_herb=2000, depth=1.5, mx_dead_pct=40),
    _sb(107, "GR7", "High load dry climate grass",
        t1=1.0, t10=0.0, t100=0.0, t_herb=5.40, dynamic=True,
        sav_1h=2000, sav_herb=1800, depth=3.0, mx_dead_pct=15),
    _sb(108, "GR8", "High load very coarse humid climate grass",
        t1=0.50, t10=1.0, t100=0.0, t_herb=7.30, dynamic=True,
        sav_1h=1500, sav_herb=1300, depth=4.0, mx_dead_pct=30),
    _sb(109, "GR9", "Very high load humid climate grass",
        t1=1.0, t10=1.0, t100=0.0, t_herb=9.0, dynamic=True,
        sav_1h=1800, sav_herb=1600, depth=5.0, mx_dead_pct=40),
    # Grass-Shrub (all dynamic)
    _sb(121, "GS1", "Low load dry climate grass-shrub",
        t1=0.20, t10=0.0, t100=0.0, t_herb=0.50, t_woody=0.65, dynamic=True,
        sav_1h=2000, sav_herb=1800, sav_woody=1800, depth=0.9, mx_dead_pct=15),
    _sb(122, "GS2", "Moderate load dry climate grass-shrub",
        t1=0.50, t10=0.50, t100=0.0, t_herb=0.60, t_woody=1.0, dynamic=True,
        sav_1h=2000, sav_herb=1800, sav_woody=1800, depth=1.5, mx_dead_pct=15),
    _sb(123, "GS3", "Moderate load humid climate grass-shrub",
        t1=0.30, t10=0.25, t100=0.0, t_herb=1.45, t_woody=1.25, dynamic=True,
        sav_1h=1800, sav_herb=1600, sav_woody=1600, depth=1.8, mx_dead_pct=40),
    _sb(124, "GS4", "High load humid climate grass-shrub",
        t1=1.90, t10=0.30, t100=0.10, t_herb=3.40, t_woody=7.10, dynamic=True,
        sav_1h=1800, sav_herb=1600, sav_woody=1600, depth=2.1, mx_dead_pct=40),
    # Shrub
    _sb(141, "SH1", "Low load dry climate shrub",
        t1=0.25, t10=0.25, t100=0.0, t_herb=0.15, t_woody=1.30, dynamic=True,
        sav_1h=2000, sav_herb=1800, sav_woody=1600, depth=1.0, mx_dead_pct=15),
    _sb(142, "SH2", "Moderate load dry climate shrub",
        t1=1.35, t10=2.40, t100=0.75, t_woody=3.85, dynamic=False,
        sav_1h=2000, sav_woody=1600, depth=1.0, mx_dead_pct=15),
    _sb(143, "SH3", "Moderate load humid climate shrub",
        t1=0.45, t10=3.0, t100=0.0, t_woody=6.20, dynamic=False,
        sav_1h=1600, sav_woody=1400, depth=2.4, mx_dead_pct=40),
    _sb(144, "SH4", "Low load humid climate timber-shrub",
        t1=0.85, t10=1.15, t100=0.20, t_woody=2.55, dynamic=False,
        sav_1h=2000, sav_herb=1800, sav_woody=1600, depth=3.0, mx_dead_pct=30),
    _sb(145, "SH5", "High load dry climate shrub",
        t1=3.60, t10=2.10, t100=0.0, t_woody=2.90, dynamic=False,
        sav_1h=750, sav_woody=1600, depth=6.0, mx_dead_pct=15),
    _sb(146, "SH6", "Low load humid climate shrub",
        t1=2.90, t10=1.45, t100=0.0, t_woody=1.40, dynamic=False,
        sav_1h=750, sav_woody=1600, depth=2.0, mx_dead_pct=30),
    _sb(147, "SH7", "Very high load dry climate shrub",
        t1=3.50, t10=5.30, t100=2.20, t_woody=3.40, dynamic=False,
        sav_1h=750, sav_woody=1600, depth=6.0, mx_dead_pct=15),
    _sb(148, "SH8", "High load humid climate shrub",
        t1=2.05, t10=3.40, t100=0.85, t_woody=4.35, dynamic=False,
        sav_1h=750, sav_woody=1600, depth=3.0, mx_dead_pct=40),
    _sb(149, "SH9", "Very high load humid climate shrub",
        t1=4.50, t10=2.45, t100=0.0, t_herb=1.55, t_woody=7.05, dynamic=True,
        sav_1h=750, sav_herb=1800, sav_woody=1500, depth=4.4, mx_dead_pct=40),
    # Timber-Understory
    _sb(161, "TU1", "Light load dry climate timber-grass-shrub",
        t1=0.20, t10=0.90, t100=1.50, t_herb=0.20, t_woody=0.90, dynamic=True,
        sav_1h=2000, sav_herb=1800, sav_woody=1600, depth=0.6, mx_dead_pct=20),
    _sb(162, "TU2", "Moderate load humid climate timber-shrub",
        t1=0.95, t10=1.80, t100=1.25, t_woody=0.20, dynamic=False,
        sav_1h=2000, sav_woody=1600, depth=1.0, mx_dead_pct=30),
    _sb(163, "TU3", "Moderate load humid climate timber-grass-shrub",
        t1=1.10, t10=0.15, t100=0.25, t_herb=0.65, t_woody=1.10, dynamic=True,
        sav_1h=1800, sav_herb=1600, sav_woody=1400, depth=1.3, mx_dead_pct=30),
    _sb(164, "TU4", "Dwarf conifer with understory",
        t1=4.50, t10=0.0, t100=0.0, t_woody=2.0, dynamic=False,
        sav_1h=2300, sav_woody=2000, depth=0.5, mx_dead_pct=12),
    _sb(165, "TU5", "Very high load dry climate timber-shrub",
        t1=4.0, t10=4.0, t100=3.0, t_woody=3.0, dynamic=False,
        sav_1h=1500, sav_woody=750, depth=1.0, mx_dead_pct=25),
    # Timber Litter
    _sb(181, "TL1", "Low load compact conifer litter",
        t1=1.0, t10=2.20, t100=3.60, dynamic=False,
        sav_1h=2000, depth=0.2, mx_dead_pct=30),
    _sb(182, "TL2", "Low load broadleaf litter",
        t1=1.40, t10=2.30, t100=2.20, dynamic=False,
        sav_1h=2000, depth=0.2, mx_dead_pct=25),
    _sb(183, "TL3", "Moderate load conifer litter",
        t1=0.50, t10=2.20, t100=2.80, dynamic=False,
        sav_1h=2000, depth=0.3, mx_dead_pct=20),
    _sb(184, "TL4", "Small downed logs",
        t1=0.50, t10=1.50, t100=4.20, dynamic=False,
        sav_1h=2000, depth=0.4, mx_dead_pct=25),
    _sb(185, "TL5", "High load conifer litter",
        t1=1.15, t10=2.50, t100=4.40, dynamic=False,
        sav_1h=2000, sav_woody=1600, depth=0.6, mx_dead_pct=25),
    _sb(186, "TL6", "Moderate load broadleaf litter",
        t1=2.40, t10=1.20, t100=1.20, dynamic=False,
        sav_1h=2000, depth=0.3, mx_dead_pct=25),
    _sb(187, "TL7", "Large downed logs",
        t1=0.30, t10=1.40, t100=8.10, dynamic=False,
        sav_1h=2000, depth=0.4, mx_dead_pct=25),
    _sb(188, "TL8", "Long-needle litter",
        t1=5.80, t10=1.40, t100=1.10, dynamic=False,
        sav_1h=1800, depth=0.3, mx_dead_pct=35),
    _sb(189, "TL9", "Very high load broadleaf litter",
        t1=6.65, t10=3.30, t100=4.15, dynamic=False,
        sav_1h=1800, sav_woody=1600, depth=0.6, mx_dead_pct=35),
    # Slash-Blowdown
    _sb(201, "SB1", "Low load activity fuel",
        t1=1.50, t10=3.0, t100=11.0, dynamic=False,
        sav_1h=2000, depth=1.0, mx_dead_pct=25),
    _sb(202, "SB2", "Moderate load activity or low load blowdown",
        t1=4.50, t10=4.25, t100=4.0, dynamic=False,
        sav_1h=2000, depth=1.0, mx_dead_pct=25),
    _sb(203, "SB3", "High load activity or moderate load blowdown",
        t1=5.50, t10=2.75, t100=3.0, dynamic=False,
        sav_1h=2000, depth=1.2, mx_dead_pct=25),
    _sb(204, "SB4", "High load blowdown",
        t1=5.25, t10=3.50, t100=5.25, dynamic=False,
        sav_1h=2000, depth=2.7, mx_dead_pct=25),
]}


# Combined lookup: number -> model (1-13 and 91-204 don't collide).
ALL: dict[int, FuelModel] = {**STANDARD_13, **STANDARD_40}
_BY_CODE: dict[str, FuelModel] = {fm.code.upper(): fm for fm in ALL.values()}


def get(key: int | str) -> FuelModel:
    """Look up any standard fuel model by number or code.

    Accepts a number (1-13 for Anderson 1982, or 91-204 for Scott & Burgan 2005)
    or a code string such as ``"GR1"`` / ``"TL3"``.
    """
    if isinstance(key, str):
        try:
            return _BY_CODE[key.strip().upper()]
        except KeyError:
            raise KeyError(f"Unknown fuel model code {key!r}") from None
    try:
        return ALL[key]
    except KeyError:
        raise KeyError(
            f"Unknown fuel model number {key!r}; expected 1-13 or 91-204"
        ) from None
