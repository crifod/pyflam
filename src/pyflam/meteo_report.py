"""Near-real-time meteorological variation report for a fire run.

Samples an :mod:`pyflam.atmosphere` provider across the run window and tracks how
the fire-weather drivers *change* over time -- the situational picture an analyst
wants before and during a run:

* surface state -- temperature, relative humidity, wind speed and direction;
* dead fuel moisture **per time lag** (1/10/100-h), evolved with the time-lag
  (Nelson-type) model so the larger fuels lag the weather;
* atmospheric conditions for **convection** -- CAPE, CIN, boundary-layer height,
  Monin-Obukhov stability class;
* potential **energy fluxes** -- surface sensible and latent heat flux, and the
  convective plume-enhancement factor those imply.

The report keeps the full time series and a per-variable *variation* summary
(min / max / mean / range / net change / trend), so the headline is what is
shifting and how fast -- e.g. RH dropping, wind backing, fuels drying, the
atmosphere destabilising.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta

import numpy as np

from . import atmosphere as atm

# Variables tracked through the window (numeric ones get a variation summary).
_NUMERIC_VARS = [
    "temperature", "relative_humidity", "wind_speed", "wind_direction",
    "m_1h", "m_10h", "m_100h", "cape", "cin", "boundary_layer_height",
    "obukhov_length", "sensible_heat_flux", "latent_heat_flux",
    "plume_factor",
]
_UNITS = {
    "temperature": "C", "relative_humidity": "%", "wind_speed": "m/s",
    "wind_direction": "deg", "m_1h": "frac", "m_10h": "frac", "m_100h": "frac",
    "cape": "J/kg", "cin": "J/kg", "boundary_layer_height": "m",
    "obukhov_length": "m", "sensible_heat_flux": "W/m2",
    "latent_heat_flux": "W/m2", "plume_factor": "x",
}


@dataclass
class MeteoReport:
    """Time series + variation summary of the fire-weather drivers."""

    times: list
    series: dict          # variable name -> list of values (None where missing)
    stability: list       # stability class per time (categorical)

    def variation(self) -> dict:
        """Per numeric variable: ``{min, max, mean, range, change, trend}``.

        ``change`` is last - first; ``trend`` is "rising"/"falling"/"steady".
        """
        out = {}
        for var, vals in self.series.items():
            v = np.array([x for x in vals if x is not None], dtype=float)
            if v.size == 0:
                continue
            first = next(x for x in vals if x is not None)
            last = next(x for x in reversed(vals) if x is not None)
            rng = float(v.max() - v.min())
            change = float(last - first)
            tol = 1e-9 if rng == 0 else 0.02 * rng
            trend = ("rising" if change > tol else
                     "falling" if change < -tol else "steady")
            out[var] = {"min": float(v.min()), "max": float(v.max()),
                        "mean": float(v.mean()), "range": rng,
                        "change": change, "trend": trend}
        return out

    def to_records(self) -> list:
        """One dict per timestep (time + all variables + stability)."""
        recs = []
        for i, t in enumerate(self.times):
            r = {"time": t, "stability": self.stability[i]}
            for var, vals in self.series.items():
                r[var] = vals[i]
            recs.append(r)
        return recs

    def summary(self) -> str:
        var = self.variation()
        n = len(self.times)
        span = ""
        if self.times and self.times[0] is not None:
            span = f" {self.times[0]} -> {self.times[-1]}"
        lines = [f"Meteo variation over {n} steps{span}:"]
        order = ["temperature", "relative_humidity", "wind_speed",
                 "wind_direction", "m_1h", "m_10h", "m_100h", "cape",
                 "boundary_layer_height", "sensible_heat_flux", "plume_factor"]
        for v in order:
            if v not in var:
                continue
            s = var[v]
            u = _UNITS.get(v, "")
            lines.append(
                f"  {v:20s} {s['min']:8.2f} .. {s['max']:8.2f} {u:5s} "
                f"(mean {s['mean']:.2f}, change {s['change']:+.2f}, {s['trend']})")
        if self.stability:
            changes = [self.stability[0]] + [
                b for a, b in zip(self.stability, self.stability[1:]) if a != b]
            lines.append(f"  stability            {' -> '.join(changes)}")
        return "\n".join(lines)


def meteo_report(
    atmosphere,
    *,
    location,
    start_time: datetime,
    end_time: datetime | None = None,
    total_minutes: float | None = None,
    step_minutes: float = 60.0,
    wind_reduction_factor: float = 0.4,
    evolve_fuel_moisture: bool = True,
    z0: float = 0.1,
) -> MeteoReport:
    """Build a :class:`MeteoReport` by sampling the atmosphere across the window.

    ``atmosphere`` is an :class:`pyflam.atmosphere.AtmosphereProvider`,
    ``location`` a ``(lat, lon)``. The window runs from ``start_time`` for
    ``total_minutes`` (or to ``end_time``) at ``step_minutes`` intervals. Dead
    fuel moisture is evolved with the time-lag model (``evolve_fuel_moisture``);
    set ``False`` to report the instantaneous equilibrium each step instead.
    """
    if total_minutes is None:
        if end_time is None:
            raise ValueError("give end_time or total_minutes")
        total_minutes = (end_time - start_time).total_seconds() / 60.0
    lat, lon = location

    series = {v: [] for v in _NUMERIC_VARS}
    times, stability = [], []
    fm_model = None
    t = 0.0
    while t <= total_minutes + 1e-9:
        clock = start_time + timedelta(minutes=t)
        st = atmosphere.state_at(lat, lon, clock)
        times.append(clock)
        stability.append(atm.stability_class(st, z0=z0))

        if fm_model is None:
            fm_model = atm.DeadFuelMoistureModel.equilibrium(st)
            m = {"m_1h": fm_model.m_1h, "m_10h": fm_model.m_10h,
                 "m_100h": fm_model.m_100h}
        elif evolve_fuel_moisture:
            m = fm_model.update(st, step_minutes)
        else:
            m = atm.dead_fuel_moisture(st)

        ob = atm.obukhov_length(st, z0=z0)
        vals = {
            "temperature": st.temperature,
            "relative_humidity": st.relative_humidity,
            "wind_speed": st.wind_speed,
            "wind_direction": st.wind_direction,
            "m_1h": m["m_1h"], "m_10h": m["m_10h"], "m_100h": m["m_100h"],
            "cape": st.cape, "cin": st.cin,
            "boundary_layer_height": st.boundary_layer_height,
            "obukhov_length": (ob if np.isfinite(ob) else None),
            "sensible_heat_flux": st.sensible_heat_flux,
            "latent_heat_flux": st.latent_heat_flux,
            "plume_factor": atm.convective_plume_factor(st),
        }
        for v in _NUMERIC_VARS:
            series[v].append(vals.get(v))
        t += step_minutes

    return MeteoReport(times=times, series=series, stability=stability)
