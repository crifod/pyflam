"""Pyroconvection-type classification + fuel gate (shared compute core).

The column classification, the .lcp fuel sampling and the fuel-gated fireline-
intensity grid used by the daily Tuscany product (``tests/pyroconv_daily.py``) and
by the GUI's preview page, factored out as parameter-driven functions (no module
globals, no plotting) so both callers share exactly the same science.

The numbers follow Castellnou et al. (2022): standard-pressure-level heights,
the 10 MW/m (=1e4 kW/m) fire-power gate, and the LCL/ABL + mixed-layer + cap
thresholds inside :func:`pyflam.atmosphere.pyroconvection_type`.
"""

from __future__ import annotations

import numpy as np

# Standard geopotential heights (m) of the pressure levels used (ICON/GFS levels).
STD_LEVEL_HEIGHT_M = {1000: 110.0, 850: 1457.0, 700: 3012.0, 500: 5574.0}
DEFAULT_LEVELS = (850, 700, 500)
FLI_GATE_KW = 1.0e4              # 10 MW/m minimum fire power for any pyroCu
ABL_MIN_M = 600.0               # below this ABL depth, held at surface plume


def _open_grib_subset(path, bbox):
    """Open one GRIB var, subset to ``bbox`` (north, west, south, east). -> (ds, var)."""
    import xarray as xr
    n, w, s, e = bbox
    ds = xr.open_dataset(path, engine="cfgrib", backend_kwargs={"indexpath": ""})
    la = ds["latitude"].values
    ds = (ds.sel(latitude=slice(n, s), longitude=slice(w, e)) if la[0] > la[-1]
          else ds.sel(latitude=slice(s, n), longitude=slice(w, e)))
    return ds, list(ds.data_vars)[0]


def read_icon2i(files, bbox, run_dt, valid_dt, hours, levels=DEFAULT_LEVELS):
    """Read ICON-2I GRIB fields for an AOI and selected valid hours.

    ``files`` is the dict returned by
    :func:`pyflam.atmosphere.fetch_icon2i_mistral` (keys ``T850``/``T700``/``T500``,
    ``T2M``, ``TD2M``, ``U10``, ``V10``); ``bbox`` is ``(north, west, south, east)``;
    ``run_dt``/``valid_dt`` are the run and valid datetimes; ``hours`` the valid
    hours (UTC) to extract. Returns a dict with 1D ``lat``/``lon``, the time-index
    list ``idx``, and stacked fields ``T`` (per level), ``T2m``, ``Td``, ``U``, ``V``.
    """
    import numpy as np

    T = {}
    for p in levels:
        ds, v = _open_grib_subset(files[f"T{p}"], bbox)
        T[p] = ds[v].values - 273.15
    t2, v = _open_grib_subset(files["T2M"], bbox); T2m = t2[v].values - 273.15
    td, v = _open_grib_subset(files["TD2M"], bbox); Td = td[v].values - 273.15
    u, vu = _open_grib_subset(files["U10"], bbox); U = u[vu].values
    vv, vn = _open_grib_subset(files["V10"], bbox); V = vv[vn].values
    lat = t2["latitude"].values; lon = t2["longitude"].values
    sh = (t2["step"].values / np.timedelta64(1, "h")).astype(int)
    idx = []
    for h in hours:
        target = int((valid_dt.replace(hour=h) - run_dt).total_seconds() // 3600)
        m = np.where(sh == target)[0]
        if not m.size:
            raise ValueError(f"valid {h:02d}Z needs forecast step +{target} h, "
                             f"beyond this run (max +{int(sh.max())} h)")
        idx.append(int(m[0]))
    return dict(lat=lat, lon=lon, idx=idx, T=T, T2m=T2m, Td=Td, U=U, V=V)


def classify(T2m, RH, Tl, *, levels=DEFAULT_LEVELS, Z=None, fli=None,
             fli_gate_kw=FLI_GATE_KW):
    """Per-cell pyroconvection class (0..4) for one time slice.

    ``T2m`` (deg C) and ``RH`` (%) are 2D surface fields; ``Tl`` maps each pressure
    level (hPa) to its 2D temperature field (deg C). ``fli`` (optional, kW/m) gates
    each cell to a surface plume below ``fli_gate_kw``. Returns an ``int16`` array
    of class *levels* (via ``PYROCONVECTION_TYPE_LEVEL``).
    """
    from pyflam.atmosphere import (
        lcl_height_m, theta_kelvin, pyroconvection_type,
        PYROCONVECTION_TYPE_LEVEL)

    Z = Z or STD_LEVEL_HEIGHT_M
    levels = tuple(levels)
    nlat, nlon = np.asarray(T2m).shape
    out = np.zeros((nlat, nlon), np.int16)

    th_sfc = theta_kelvin(T2m, 1000.0)
    th = {p: theta_kelvin(Tl[p], float(p)) for p in levels}
    ml_grad = (th[850] - th_sfc) / (Z[850] - Z[1000])
    gamma = (th[500] - th[700]) / (Z[500] - Z[700])
    lcl = lcl_height_m(T2m, RH)

    abl = np.full((nlat, nlon), Z[700])
    for (zlo, tlo), (zhi, thi) in [((Z[1000], th_sfc), (Z[850], th[850])),
                                   ((Z[850], th[850]), (Z[700], th[700]))]:
        crossed = (thi > th_sfc + 0.5) & (abl == Z[700])
        denom = np.where(thi != tlo, thi - tlo, 1e9)
        frac = np.clip((th_sfc + 0.5 - tlo) / denom, 0.0, 1.0)
        abl = np.where(crossed, zlo + frac * (zhi - zlo), abl)
    abl = np.maximum(abl, 100.0)
    ratio = lcl / abl

    for a in range(nlat):
        for b in range(nlon):
            if abl[a, b] < ABL_MIN_M:
                continue
            fkw = None if fli is None else float(fli[a, b])
            out[a, b] = PYROCONVECTION_TYPE_LEVEL[pyroconvection_type(
                lcl_abl_ratio=float(ratio[a, b]),
                ml_theta_gradient=float(ml_grad[a, b]),
                gamma_theta=float(gamma[a, b]),
                fireline_intensity_kw=fkw, fli_threshold_kw=fli_gate_kw)]
    return out


def _is_known_fuel(n) -> bool:
    from pyflam import fuel_models
    try:
        fuel_models.get(int(n))
        return True
    except KeyError:
        return False


def lcp_fields(ls, lat, lon, transformer):
    """Sample a landscape's fuel/canopy bands onto a ``(lat, lon)`` weather grid.

    ``ls`` is a :class:`pyflam.Landscape`; ``lat``/``lon`` are 1D coordinate arrays
    of the weather grid; ``transformer`` is a ``pyproj``-style object with a
    ``.transform(lon, lat) -> (x, y)`` method into the landscape's CRS. Returns
    ``(fields_dict, burnable_mask)`` where ``fields_dict`` has ``fuel``, ``slope``,
    ``cbh``, ``cbd`` and ``ch`` 2D arrays on the weather grid (canopy bands rescaled
    from .lcp integer units to SI), and ``burnable_mask`` is indexed by fuel number.
    """
    from pyflam import fuel_models

    n = int(np.asarray(ls.fuel_model).max()) + 1
    burn = np.array([fuel_models.get(int(k)).is_burnable if _is_known_fuel(k) else False
                     for k in range(n)])
    lon2d, lat2d = np.meshgrid(lon, lat)
    x, y = transformer.transform(lon2d.ravel(), lat2d.ravel())
    c = ((x - ls.west) / ls.cellsize_x).astype(int)
    r = ((ls.north - y) / ls.cellsize_y).astype(int)
    ins = (r >= 0) & (r < ls.shape[0]) & (c >= 0) & (c < ls.shape[1])
    rr = np.clip(r, 0, ls.shape[0] - 1)
    cc = np.clip(c, 0, ls.shape[1] - 1)
    fn = np.full(r.shape, -1)
    fn[ins] = np.asarray(ls.fuel_model)[rr[ins], cc[ins]]
    sh = lat2d.shape

    def band(arr, scale=1.0):
        return (np.asarray(arr)[rr, cc] / scale).reshape(sh)

    fields = dict(
        fuel=fn.reshape(sh),
        slope=ls.slope_tangent[rr, cc].reshape(sh),
        cbh=band(ls.canopy_base_height, 10.0) if ls.canopy_base_height is not None else np.zeros(sh),
        cbd=band(ls.canopy_bulk_density, 100.0) if ls.canopy_bulk_density is not None else np.zeros(sh),
        ch=band(ls.canopy_height, 10.0) if ls.canopy_height is not None else np.zeros(sh),
    )
    return fields, burn


def fli_grid(T2m, RH, wsp, lf, burn, *, m_live_herb=0.70, m_live_woody=0.90,
             wind_reduction_factor=0.4):
    """Fuel-gated fireline-intensity grid (kW/m) from forecast moisture/wind.

    Runs Rothermel surface spread (and Cruz-2005 crown where a canopy stack
    exists) on the sampled .lcp fuels (``lf`` from :func:`lcp_fields`) using the
    forecast 2 m T/RH (equilibrium moisture) and 10 m wind speed ``wsp`` (m/s).
    """
    import pyflam
    from pyflam import units
    from pyflam.atmosphere import equilibrium_moisture_content

    nlat, nlon = np.asarray(T2m).shape
    out = np.zeros((nlat, nlon))
    for a in range(nlat):
        for b in range(nlon):
            fnum = int(lf["fuel"][a, b])
            if fnum < 0 or fnum >= burn.size or not burn[fnum]:
                continue
            fuel = pyflam.get_fuel_model(fnum)
            emc = equilibrium_moisture_content(T2m[a, b], RH[a, b]) / 100.0
            u = max(float(wsp[a, b]), 0.0)
            surf = pyflam.spread(
                fuel, m_1h=emc, m_10h=emc + 0.01, m_100h=emc + 0.02,
                m_live_herb=m_live_herb, m_live_woody=m_live_woody,
                wind_midflame=units.m_per_s_to_ft_per_min(u) * wind_reduction_factor,
                slope=float(lf["slope"][a, b]))
            fli = units.btu_per_ft_s_to_kw_per_m(surf.fireline_intensity)
            if lf["cbh"][a, b] > 0 and lf["cbd"][a, b] > 0:
                cr = pyflam.crown_fire_behavior(
                    surf, canopy_base_height=lf["cbh"][a, b],
                    canopy_bulk_density=lf["cbd"][a, b], foliar_moisture=100.0,
                    wind_20ft_ft_per_min=units.m_per_s_to_ft_per_min(u * 0.87),
                    canopy_fuel_load=max(lf["cbd"][a, b]
                                         * max(lf["ch"][a, b] - lf["cbh"][a, b], 0), 0),
                    m_1h=emc, m_10h=emc + 0.01, m_100h=emc + 0.02,
                    m_live_herb=m_live_herb, m_live_woody=m_live_woody,
                    crown_spread="cruz2005")
                fli = max(fli, cr.fireline_intensity)
            out[a, b] = fli
    return out
