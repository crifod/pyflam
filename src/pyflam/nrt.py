"""Near-real-time run product: weather -> spread -> perimeter -> reports, in one call.

:func:`run_realtime` is the operational entry point. Given a landscape, an
ignition and an atmosphere provider (a live forecast from
:func:`pyflam.atmosphere.fetch_gfs`, an ERA5 reanalysis slice, or any provider),
it:

1. builds the **meteo variation report** over the run window
   (:mod:`pyflam.meteo_report`);
2. **spreads the fire** through the evolving weather with the time-marched
   coupling (:func:`pyflam.pyroconvection.fire_atmosphere_march`) -- optionally
   with the fire-plume CFD feedback;
3. extracts the **perimeter** and runs the **operative driving-force analysis**
   by sector (:mod:`pyflam.operative`);

and returns a :class:`RunProduct` bundling the arrival-time raster, both reports,
and a one-call GeoJSON export for a mapping front-end.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import timedelta

from .atmosphere import spread_inputs_from_state
from .meteo_report import meteo_report
from .operative import analyze_perimeter
from .pyroconvection import _uniform_wind_field, fire_atmosphere_march


@dataclass
class RunProduct:
    """The bundled output of a near-real-time run."""

    arrival_time: object              # 2D array, minutes
    meteo: object                     # MeteoReport
    operative: object                 # OperativeReport (or None if nothing burned)
    location: tuple
    start_time: object
    total_time: float

    def operative_geojson(self, ls, **kwargs) -> dict:
        if self.operative is None:
            raise ValueError("no perimeter to export (nothing burned)")
        return self.operative.to_geojson(ls, **kwargs)

    def write_geojson(self, ls, path, **kwargs) -> None:
        if self.operative is None:
            raise ValueError("no perimeter to export (nothing burned)")
        self.operative.write_geojson(ls, path, **kwargs)

    def summary(self) -> str:
        import numpy as np
        burned = int(np.isfinite(self.arrival_time).sum())
        parts = [f"Near-real-time run from {self.start_time} "
                 f"for {self.total_time:g} min at {self.location}:",
                 f"  burned cells: {burned}", "", self.meteo.summary()]
        if self.operative is not None:
            parts += ["", self.operative.summary()]
        else:
            parts += ["", "  (no perimeter: fire did not establish)"]
        return "\n".join(parts)


def run_realtime(
    ls,
    ignitions,
    *,
    atmosphere,
    location,
    start_time,
    total_time: float,
    dt: float = 30.0,
    m_live_herb: float = 0.0,
    m_live_woody: float = 0.0,
    load_factor: float = 1.0,
    wind_reduction_factor: float = 0.4,
    meteo_step_minutes: float = 60.0,
    use_plume: bool = False,
    spatial: bool = False,
    subsectors: int = 1,
    **march_kwargs,
) -> RunProduct:
    """Run a complete near-real-time simulation and return a :class:`RunProduct`.

    ``atmosphere`` is a :class:`pyflam.atmosphere.AtmosphereProvider`, ``location``
    a ``(lat, lon)``, ``start_time`` a datetime. The fire spreads for
    ``total_time`` minutes in ``dt``-minute increments through the evolving
    weather. ``use_plume`` enables the buoyant-CFD fire-plume feedback (needs
    OpenFOAM and a DEM); otherwise the atmospheric wind drives spread directly.
    ``spatial`` samples the atmosphere per cell (gridded weather). The operative
    analysis uses the conditions at the analysis (final) time.
    """
    lat, lon = location

    # 1. meteo variation report over the window.
    meteo = meteo_report(
        atmosphere, location=location, start_time=start_time,
        total_minutes=total_time, step_minutes=meteo_step_minutes,
        wind_reduction_factor=wind_reduction_factor)

    # 2. spread through the evolving weather. Wind-model routing:
    #   spatial + use_plume -> per-cell atmosphere with the fire plume superposed;
    #   spatial             -> per-cell atmosphere (no plume);
    #   use_plume           -> single-column atmosphere + CFD plume coupling;
    #   else                -> atmospheric wind only (no CFD).
    if spatial:
        wind_provider = None
    elif use_plume:
        wind_provider = None          # default -> couple_fire_wind (OpenFOAM)
    else:
        def wind_provider(ls_, intensity, active, spd, dirn):
            return _uniform_wind_field(ls_, spd, dirn)

    march = fire_atmosphere_march(
        ls, ignitions, total_time=total_time, dt=dt, atmosphere=atmosphere,
        location=location, start_time=start_time, spatial=spatial,
        plume=(use_plume and spatial),
        m_live_herb=m_live_herb, m_live_woody=m_live_woody, load_factor=load_factor,
        wind_reduction_factor=wind_reduction_factor, wind_provider=wind_provider,
        **march_kwargs)
    arrival = march["arrival_time"]

    # 3. operative driving-force analysis at the final perimeter.
    final = atmosphere.state_at(lat, lon, start_time + timedelta(minutes=total_time))
    si = spread_inputs_from_state(final, wind_reduction_factor=wind_reduction_factor)
    try:
        op = analyze_perimeter(
            ls, arrival, total_time, subsectors=subsectors,
            wind_midflame=si["wind_midflame"], wind_direction=si["wind_direction"],
            m_1h=si["m_1h"], m_10h=si["m_10h"], m_100h=si["m_100h"],
            m_live_herb=m_live_herb, m_live_woody=m_live_woody,
            load_factor=load_factor)
    except ValueError:
        op = None                     # nothing burned -> no perimeter

    return RunProduct(arrival_time=arrival, meteo=meteo, operative=op,
                      location=location, start_time=start_time,
                      total_time=float(total_time))
