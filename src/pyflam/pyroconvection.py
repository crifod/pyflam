"""Fire -> atmosphere coupling: the fire's energy drives a buoyant plume that
re-shapes the local wind (pyroconvection), which in turn drives the fire.

The earlier wind models are one-way: terrain + ambient weather set a wind, and the
fire spreads in it. But an intense fire is not a passive tracer -- its convective
heat release builds a buoyant plume that pulls air inward (indrafts) and lifts it
over the front, locally strengthening and bending the wind. This module closes
that loop through the existing RANS solver:

    fireline intensity  ->  ground heat-flux field  ->  buoyant RANS (OpenFOAM)
        ->  modified near-surface wind  ->  (re-run spread)

The heat-flux step is the physical bridge. Byram's fireline intensity ``I`` is an
energy release **per unit length of fire front**; spread over a grid cell of side
``L`` the cell-averaged convective heat flux into the atmosphere is
``q = chi_c * I / L`` (W/m^2), with ``chi_c`` the convective fraction. That field
becomes a spatially-varying ground sensible-heat flux for
``buoyantBoussinesqSimpleFoam`` (see :mod:`pyflam.cfd`), so the solved wind
includes the plume the fire itself creates.

This is a **quasi-steady** coupling (a steady RANS plume for a frozen fire state),
matching the rest of :mod:`pyflam.cfd`; it captures indraft/updraft structure, not
transient pyro-convective bursts. Needs OpenFOAM (via the ``openfoam`` wrapper).

References:
    Byram, G.M. 1959. Combustion of forest fuels. In: Forest fire: control and use.
    Clark, T.L.; Coen, J.; Latham, D. 2004. Description of a coupled
        atmosphere-fire model. Int. J. Wildland Fire 13. (Fire-atmosphere coupling.)
"""

from __future__ import annotations

import math

import numpy as np

from . import cfd

# Byram fireline intensity: Btu/ft/s -> W/m (3.46414 kW/m per Btu/ft/s).
_BTU_FT_S_TO_W_M = 3464.14
DEFAULT_CONVECTIVE_FRACTION = 0.7


def fire_heat_flux(
    fireline_intensity,
    cellsize_m: float,
    *,
    convective_fraction: float = DEFAULT_CONVECTIVE_FRACTION,
    active_mask=None,
) -> np.ndarray:
    """Cell-averaged convective heat flux (W/m^2) from fireline intensity.

    ``fireline_intensity`` is Byram intensity (Btu/ft/s, per unit length of
    front), e.g. ``SpreadField.fireline_intensity``. The cell-averaged flux is
    ``q = convective_fraction * I[W/m] / cellsize``. ``active_mask`` (bool, same
    shape) restricts the flux to the actively burning cells -- pass a fire
    perimeter / front so unburned ground stays at ambient.
    """
    i_w_m = np.asarray(fireline_intensity, dtype=float) * _BTU_FT_S_TO_W_M
    q = convective_fraction * i_w_m / float(cellsize_m)
    if active_mask is not None:
        q = np.where(np.asarray(active_mask, dtype=bool), q, 0.0)
    return q


def couple_fire_wind(
    ls,
    fireline_intensity,
    *,
    speed: float,
    direction: float,
    active_mask=None,
    convective_fraction: float = DEFAULT_CONVECTIVE_FRACTION,
    ambient_heat_flux: float = 0.0,
    z0: float | None = None,
    output_height: float = 6.1,
    **cfd_kwargs,
):
    """Solve the pyroconvective wind for a fire on a landscape.

    Turns ``fireline_intensity`` (Btu/ft/s, e.g. ``spread_field(...).
    fireline_intensity``) into a ground heat-flux field and runs a buoyant RANS
    simulation, returning the modified :class:`pyflam.wind.WindField` at
    ``output_height`` AGL -- the wind including the fire's own plume. Restrict the
    heat to the active fire with ``active_mask`` (e.g. ``perimeter_mask(arrival,
    t)``). ``speed``/``direction`` are the ambient inflow; ``ambient_heat_flux``
    (W/m^2) is the atmosphere's own background surface heating (from
    :func:`pyflam.atmosphere.ambient_surface_heat_flux`) which the fire's plume
    develops into -- it is added everywhere. Extra ``cfd_kwargs`` pass through to
    :func:`pyflam.cfd.solve_rans`.

    Requires OpenFOAM. Feed the result back into :func:`pyflam.spread_field` as a
    ``wind_midflame`` field (via :meth:`WindField.to_midflame`) to let the plume
    influence the next increment of spread.
    """
    if ls.elevation is None:
        raise ValueError("landscape has no elevation band for the RANS solver")
    q = fire_heat_flux(fireline_intensity, ls.cellsize_x,
                       convective_fraction=convective_fraction,
                       active_mask=active_mask)
    if ambient_heat_flux:
        q = q + float(ambient_heat_flux)        # atmosphere's background heating
    if z0 is None:
        from .windsolver import roughness_from_fuel
        z0 = np.asarray(roughness_from_fuel(ls.fuel_model), dtype=float)
    return cfd.solve_rans(
        np.asarray(ls.elevation, dtype=float), ls.cellsize_x,
        speed=speed, direction=direction, z0=z0,
        buoyant=True, surface_heat_flux=q, output_height=output_height,
        west=ls.west, north=ls.north, crs=ls.crs, **cfd_kwargs)


def _uniform_wind_field(ls, speed, direction):
    """A uniform :class:`~pyflam.wind.WindField` (m/s) on the landscape grid."""
    from .wind import WindField
    shape = ls.shape
    return WindField(
        speed=np.full(shape, float(speed)),
        direction=np.full(shape, float(direction)),
        cellsize=ls.cellsize_x, west=ls.west, north=ls.north,
        speed_units="m/s", crs=ls.crs)


_MS_PER_FT_MIN = 0.3048 / 60.0


def _windfield_to_uv(wf):
    """Wind field -> (u east, v north) components in m/s."""
    spd = wf.speed_ft_per_min() * _MS_PER_FT_MIN
    toward = np.radians((np.asarray(wf.direction, float) + 180.0) % 360.0)
    return spd * np.sin(toward), spd * np.cos(toward)


def _uv_to_windfield(u, v, like):
    """(u, v) m/s components -> a WindField (m/s) on ``like``'s grid."""
    from .wind import WindField
    spd = np.hypot(u, v)
    toward = np.degrees(np.arctan2(u, v)) % 360.0
    return WindField(speed=spd, direction=(toward + 180.0) % 360.0,
                     cellsize=like.cellsize, west=like.west, north=like.north,
                     speed_units="m/s", crs=like.crs)


def _stabilize_wind(wf_new, wf_prev, ls, *, relax, max_factor, ambient_mean_ms):
    """Damp and cap the plume-fed wind to keep the coupling feedback bounded.

    Under-relaxation blends this increment's wind with the previous one
    (vectorially, in u/v) -- ``relax=1`` is no damping; smaller damps harder. The
    cap limits each cell's speed to ``max_factor * ambient mean speed``, so the
    crowning -> plume -> wind -> crown loop cannot run away. A no-op at the defaults
    (``relax=1``, ``max_factor=inf``).
    """
    wf = wf_new if wf_new.shape == ls.shape else wf_new.to_landscape(ls)
    if relax < 1.0 and wf_prev is not None:
        wp = wf_prev if wf_prev.shape == ls.shape else wf_prev.to_landscape(ls)
        up, vp = _windfield_to_uv(wp)
        un, vn = _windfield_to_uv(wf)
        wf = _uv_to_windfield((1.0 - relax) * up + relax * un,
                              (1.0 - relax) * vp + relax * vn, like=wf)
    if math.isfinite(max_factor) and ambient_mean_ms > 0.0:
        cap = max_factor * ambient_mean_ms
        u, v = _windfield_to_uv(wf)
        spd = np.hypot(u, v)
        scale = np.where(spd > cap, cap / np.maximum(spd, 1e-9), 1.0)
        wf = _uv_to_windfield(u * scale, v * scale, like=wf)
    return wf


def superpose_plume(ambient_wf, base_wf, fire_wf):
    """Add a plume perturbation (fire - base) onto an ambient wind field.

    All three are :class:`~pyflam.wind.WindField` on the same grid; returns the
    merged field ``ambient + (fire - base)`` in (u, v) space. This is the linear
    superposition that lets a *spatially varying* atmospheric wind and the fire's
    own plume coexist (the plume perturbation is computed at a single inflow but
    added to the per-cell ambient wind).
    """
    au, av = _windfield_to_uv(ambient_wf)
    bu, bv = _windfield_to_uv(base_wf)
    fu, fv = _windfield_to_uv(fire_wf)
    return _uv_to_windfield(au + (fu - bu), av + (fv - bv), ambient_wf)


def merge_plume_wind(ls, ambient_windfield, fireline_intensity, *,
                     active_mask=None, convective_fraction=DEFAULT_CONVECTIVE_FRACTION,
                     z0=None, output_height=6.1, **cfd_kwargs):
    """Spatial ambient wind + the fire's plume, via CFD perturbation superposition.

    Runs the buoyant RANS twice at the ambient *mean* inflow -- once without the
    fire and once with it -- and adds the difference (the plume's indraft/updraft
    perturbation) to the per-cell ``ambient_windfield``. So the result carries
    both the spatially-varying forecast wind and the fire-induced plume flow.
    Needs OpenFOAM.
    """
    import math
    au, av = _windfield_to_uv(ambient_windfield)
    mu, mv = float(np.nanmean(au)), float(np.nanmean(av))
    mean_speed = math.hypot(mu, mv)
    mean_from = (math.degrees(math.atan2(mu, mv)) + 180.0) % 360.0
    common = dict(speed=mean_speed, direction=mean_from, z0=z0,
                  output_height=output_height, ambient_heat_flux=0.0, **cfd_kwargs)
    base = couple_fire_wind(ls, np.zeros_like(np.asarray(fireline_intensity, float)),
                            **common)
    fire = couple_fire_wind(ls, fireline_intensity, active_mask=active_mask,
                            convective_fraction=convective_fraction, **common)
    return superpose_plume(ambient_windfield, base, fire)


def _field_from_wind(ls, wf, *, wind_reduction_factor, moist, load_factor):
    """Build a spread field from a (possibly per-cell) WindField."""
    from .mtt import spread_field
    wfl = wf if wf.shape == ls.shape else wf.to_landscape(ls)
    midflame = wfl.speed_ft_per_min() * wind_reduction_factor
    return spread_field(ls, wind_midflame=midflame, wind_direction=wfl.direction,
                        load_factor=load_factor, **moist)


def fire_atmosphere_march(
    ls,
    ignitions,
    *,
    total_time: float,
    dt: float,
    speed: float | None = None,
    direction: float | None = None,
    m_1h: float | None = None,
    m_10h: float | None = None,
    m_100h: float | None = None,
    m_live_herb: float = 0.0,
    m_live_woody: float = 0.0,
    load_factor: float = 1.0,
    wind_reduction_factor: float = 0.4,
    flame_residence: float | None = None,
    ring: int = 2,
    wind_provider=None,
    ambient_wind_provider=None,
    atmosphere=None,
    location=None,
    start_time=None,
    spatial: bool = False,
    plume: bool = False,
    crown: bool = False,
    crown_spread: str = "cruz2005",
    foliar_moisture: float | None = None,
    cbh_scale: float = 0.1,
    cbd_scale: float = 0.01,
    crown_heat_content: float = 18000.0,
    wind_relax: float = 1.0,
    max_wind_factor: float = math.inf,
    pyroconvection: bool = False,
    profile=None,
    state=None,
    return_history: bool = False,
    **cfd_kwargs,
):
    """Time-march fire growth coupled to the fire's own plume wind.

    Grows the fire in ``dt``-minute increments; before each increment the wind is
    re-solved from the *current* fire state (the active front's fireline intensity
    drives a buoyant plume), so spread responds to the wind the fire is making and
    vice-versa. The first increment uses the ambient wind (no fire yet).

    ``ignitions`` is an iterable of ``(row, col)``. ``wind_provider(ls, intensity,
    active_mask, speed, direction)`` returns the plume-modified
    :class:`~pyflam.wind.WindField`; it defaults to :func:`couple_fire_wind` (needs
    OpenFOAM) but can be injected (e.g. for testing or a different solver).
    ``ambient_wind_provider(ls, speed, direction)`` gives the no-fire wind
    (default: uniform). ``flame_residence`` (minutes; default ``dt``) sets how long
    a cell keeps emitting heat after it burns (the active flaming band).

    **Atmosphere-driven runs.** Pass an ``atmosphere`` provider (see
    :mod:`pyflam.atmosphere`) with ``location=(lat, lon)`` and a ``start_time``
    (datetime); each increment then re-reads the atmospheric state at the
    advancing clock time and derives the wind, dead fuel moisture and background
    convective heat flux from it -- so the fire responds to evolving weather, for
    a near-real-time forecast or a reanalysis. Without ``atmosphere`` the fixed
    ``speed``/``direction``/``m_1h``/``m_10h``/``m_100h`` are used throughout.
    With ``spatial=True`` the atmosphere is sampled **per cell** (gridded weather:
    wind and fuel moisture vary across the domain). With ``spatial=True`` *and*
    ``plume=True`` the fire's plume perturbation (a buoyant-CFD fire-vs-no-fire
    difference) is **superposed** onto the per-cell atmospheric wind
    (:func:`merge_plume_wind`) -- spatial weather and the fire's own plume
    together. (``spatial=False`` with the default CFD ``wind_provider`` is the
    single-column plume-coupling path.)

    **Crown fire (plan step 3).** With ``crown=True`` each increment builds a
    *crown-aware* spread field (:func:`pyflam.crown_spread_field`) from the current
    plume-modified wind: cells that crown spread at the crown rate (``crown_spread``,
    default Cruz 2005) and carry the much-higher crown fireline intensity -- which
    feeds straight into the plume and the ember spotting, closing the
    crowning -> stronger plume -> faster crown feedback. Needs ``foliar_moisture``
    and canopy base-height / bulk-density bands on ``ls``; the output then also
    carries the final ``fire_type`` raster.

    **Feedback stability.** The plume->wind->spread loop is positive feedback, so two
    bounds keep it from running away: ``wind_relax`` (0-1) under-relaxes the wind by
    blending each increment with the previous (``1`` = none, default; ``~0.5`` damps),
    and ``max_wind_factor`` caps each cell's speed at that multiple of the ambient
    mean wind (default ``inf`` = no cap). ``return_history`` adds ``mean_wind`` per
    increment so a run's boundedness can be checked.

    **Pyroconvection.** With ``pyroconvection=True`` the fireline intensity fed to
    the plume each step is scaled by :func:`pyflam.convective_plume_factor` -- so an
    unstable, dry, pyroCb-prone atmosphere drives a *stronger* plume (and hence
    spotting/indraft), a stable one a weaker one. Pass a vertical ``profile``
    (:class:`pyflam.AtmosphericProfile`) to use the literature's pyroCb predictors
    (LCL / inverted-V / Continuous Haines) rather than surface CAPE alone. The
    surface state is the ``atmosphere`` state (scalar/point runs) or a manually
    supplied ``state`` (:class:`pyflam.AtmosphericState`); the output then also
    carries the ``pyroconvection`` potential and (with a profile) the
    ``pyrocb_firepower_threshold``. Only the plume coupling is scaled -- the spread
    field's own intensity (flame length, etc.) is unchanged.

    Returns a dict with ``arrival_time`` (minutes); with ``return_history`` also
    ``winds``/``fields``/``times``/``mean_wind``/``plume_factor`` lists, one per
    increment.
    """
    if atmosphere is None and None in (speed, direction, m_1h, m_10h, m_100h):
        raise ValueError(
            "provide speed/direction/m_1h/m_10h/m_100h, or an `atmosphere` "
            "provider with `location` (and `start_time` for time-varying weather)")
    if crown:
        if foliar_moisture is None:
            raise ValueError("crown=True needs foliar_moisture")
        if ls.canopy_base_height is None or ls.canopy_bulk_density is None:
            raise ValueError("crown=True needs canopy_base_height and "
                             "canopy_bulk_density bands on the landscape")

    from datetime import timedelta

    from .mtt import (
        build_traveltime_graph, minimum_travel_time, _dijkstra_with_start_times,
    )

    ambient_q = [0.0]                      # current background heat flux (W/m^2)
    if wind_provider is None:
        def wind_provider(ls_, intensity, active_mask, spd, dirn):
            return couple_fire_wind(ls_, intensity, speed=spd, direction=dirn,
                                    active_mask=active_mask,
                                    ambient_heat_flux=ambient_q[0], **cfd_kwargs)
    if ambient_wind_provider is None:
        ambient_wind_provider = _uniform_wind_field
    if flame_residence is None:
        flame_residence = dt

    fixed_moist = dict(m_1h=m_1h, m_10h=m_10h, m_100h=m_100h,
                       m_live_herb=m_live_herb, m_live_woody=m_live_woody)

    crown_state = {"fire_type": None, "crown_fraction_burned": None}

    def build_field(wf, moist):
        """Surface spread field, or a crown-aware one when ``crown=True``."""
        if not crown:
            return _field_from_wind(
                ls, wf, wind_reduction_factor=wind_reduction_factor,
                moist=moist, load_factor=load_factor)
        from .crownfire import crown_spread_field
        wfl = wf if wf.shape == ls.shape else wf.to_landscape(ls)
        w20 = wfl.speed_ft_per_min()                  # the 20-ft (6.1 m) wind
        caf = crown_spread_field(
            ls, wind_midflame=w20 * wind_reduction_factor,
            wind_direction=wfl.direction, wind_20ft_ft_per_min=w20,
            foliar_moisture=foliar_moisture, crown_spread=crown_spread,
            load_factor=load_factor, cbh_scale=cbh_scale, cbd_scale=cbd_scale,
            heat_content=crown_heat_content, **moist)
        crown_state["fire_type"] = caf.fire_type
        crown_state["crown_fraction_burned"] = caf.crown_fraction_burned
        return caf.field

    pyro_state = [state]                  # representative AtmosphericState (updated below)

    def resolve(sim_min):
        """(moist, speed, direction, ambient_heat_flux, windfield) for this time.

        ``windfield`` is a per-cell WindField in ``spatial`` mode (the atmosphere
        varies across the domain), else ``None`` and scalar speed/direction apply.
        Also stashes the scalar/point atmospheric state for the pyroconvection
        plume factor (left as the manual ``state`` in spatial mode).
        """
        if atmosphere is None:
            return fixed_moist, speed, direction, 0.0, None
        from .atmosphere import (
            ambient_surface_heat_flux, dead_fuel_moisture, wind_field_from_state,
        )
        clock = (start_time + timedelta(minutes=sim_min)
                 if start_time is not None else None)
        m_live = {"m_live_herb": m_live_herb, "m_live_woody": m_live_woody}
        if spatial:
            fld = atmosphere.field_on(ls, clock)
            m = {**dead_fuel_moisture(fld), **m_live}
            aq = float(np.mean(ambient_surface_heat_flux(fld)))
            return m, None, None, aq, wind_field_from_state(fld, ls)
        lat, lon = (location or (None, None))
        st = atmosphere.state_at(lat, lon, clock)
        pyro_state[0] = st
        m = {**dead_fuel_moisture(st), **m_live}
        return m, st.wind_speed, st.wind_direction, \
            ambient_surface_heat_flux(st), None

    def plume_factor():
        """Convective plume-intensity multiplier for this step (>=1 pyroconvective)."""
        if not pyroconvection or pyro_state[0] is None:
            return 1.0
        from .atmosphere import convective_plume_factor
        return float(convective_plume_factor(pyro_state[0], profile=profile))

    nrows, ncols = ls.shape
    n = nrows * ncols

    def _mean_ms(w):
        return float(np.nanmean(w.speed_ft_per_min())) * _MS_PER_FT_MIN

    def _ambient_mean_ms(spd_scalar, wf_a):
        if wf_a is not None:
            return _mean_ms(wf_a)
        return float(spd_scalar) if spd_scalar is not None else 0.0

    # Increment 0: ambient wind, no plume.
    moist, spd, dirn, ambient_q[0], wf_atm = resolve(0.0)
    wf = wf_atm if wf_atm is not None else ambient_wind_provider(ls, spd, dirn)
    field = build_field(wf, moist)
    arrival = minimum_travel_time(field, ignitions, max_time=dt, ring=ring)

    wf_prev = wf
    history = {"winds": [wf], "fields": [field], "times": [dt],
               "mean_wind": [_mean_ms(wf)], "plume_factor": [plume_factor()]}
    t = dt
    while t < total_time - 1e-9:
        t_next = min(t + dt, total_time)
        burned = np.isfinite(arrival) & (arrival <= t)
        if not burned.any():
            break
        # Refresh the atmospheric forcing for this increment (evolving weather).
        moist, spd, dirn, ambient_q[0], wf_atm = resolve(t)
        # Active flaming band drives the plume; re-solve the wind from it. A
        # pyroconvective atmosphere strengthens the plume (scale only the coupling
        # input, not the spread field's own intensity).
        active = burned & (arrival > t - flame_residence)
        pf = plume_factor()
        plume_intensity = pf * np.asarray(field.fireline_intensity, dtype=float)
        if wf_atm is not None:                # spatial atmosphere
            wf = (merge_plume_wind(ls, wf_atm, plume_intensity,
                                   active_mask=active, **cfd_kwargs)
                  if plume else wf_atm)       # merge the plume onto the ambient field
        else:
            wf = wind_provider(ls, plume_intensity, active, spd, dirn)
        # Bound the feedback: under-relax against the previous wind and cap the speed.
        wf = _stabilize_wind(wf, wf_prev, ls, relax=wind_relax,
                             max_factor=max_wind_factor,
                             ambient_mean_ms=_ambient_mean_ms(spd, wf_atm))
        wf_prev = wf
        field = build_field(wf, moist)
        graph = build_traveltime_graph(field, ring=ring)
        rb, cb = np.where(burned)
        sources = (rb * ncols + cb).tolist()
        starts = arrival[rb, cb].tolist()
        extended = _dijkstra_with_start_times(
            graph, n, sources, starts, t_next).reshape(ls.shape)
        arrival = np.minimum(arrival, extended)
        if return_history:
            history["winds"].append(wf)
            history["fields"].append(field)
            history["times"].append(t_next)
            history["mean_wind"].append(_mean_ms(wf))
            history["plume_factor"].append(pf)
        t = t_next

    out = {"arrival_time": arrival}
    if pyroconvection and pyro_state[0] is not None:
        from .atmosphere import pyroconvection_potential, pyrocb_firepower_threshold
        out["pyroconvection"] = pyroconvection_potential(pyro_state[0], profile=profile)
        if profile is not None:
            out["pyrocb_firepower_threshold"] = pyrocb_firepower_threshold(
                pyro_state[0], profile)
    if crown:
        out["fire_type"] = crown_state["fire_type"]
        out["crown_fraction_burned"] = crown_state["crown_fraction_burned"]
    if return_history:
        out.update(history)
    return out
