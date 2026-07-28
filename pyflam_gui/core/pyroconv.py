# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (c) 2026 Cristiano Foderi <cristiano.foderi@gmail.com>

"""Pyroconvection-type classification + fuel gate (shared compute core).

The column classification, the .lcp fuel sampling and the fuel-gated fireline-
intensity grid used by the daily Tuscany product (``tests/pyroconv_daily.py``) and
by the GUI's preview page, factored out as parameter-driven functions (no module
globals, no plotting) so both callers share exactly the same science.

Two classification paths live here:

* :func:`classify` -- the original three-diagnostic path. Standard-atmosphere level
  heights, a theta-excess ABL, and the published Castellnou ladder. Needs only
  temperature on 3 pressure levels, so it works on any T-only source (GFS, the old
  ICON-2I field set) and is what the reference cases are validated against.
* :func:`classify_profile` -- the profile path. Real per-cell heights from
  geopotential, a bulk-Richardson ABL, an exact (Bolton) LCL, moisture at the ABL
  top and, where the levels resolve it, the shear-maximum height. Needs the full
  :data:`pyflam.atmosphere.ICON2I_PROFILE_FIELDS` download. This is what the daily
  product now runs.

Both share the 10 MW/m (=1e4 kW/m) fire-power gate on the .lcp fuels, which is what
separates the *expected* map from the *potential* upper bound.
"""

from __future__ import annotations

import numpy as np

# Standard geopotential heights (m) of the pressure levels used (ICON/GFS levels).
# Only the legacy :func:`classify` path uses these; :func:`classify_profile` takes
# real heights from the model's geopotential instead.
STD_LEVEL_HEIGHT_M = {1000: 110.0, 850: 1457.0, 700: 3012.0, 500: 5574.0}
DEFAULT_LEVELS = (850, 700, 500)
FLI_GATE_KW = 1.0e4              # 10 MW/m minimum fire power for any pyroCu
# Plausibility floor on the diagnosed ABL, matching the Rib solver's own clip
# (``atmosphere._ABL_MIN_M``). It guards against a near-zero depth turning LCL/ABL into a
# meaningless quotient; it is deliberately **not** a discriminating gate.
#
# It was 600 m until 2026-07-26, held there as "mixing too shallow to classify". That had no
# basis in the source method -- Castellnou et al. (2022) sec.2.4.2 specifies the Rib procedure
# (Ri_b > 0.33, search from ~400 m AGL) and imposes no minimum depth -- and it is contradicted
# by the campaign's own data: of the eight ambient sondes released beside real pyroconvective
# wildfires (Castellnou Ribau et al. 2025, AMT 18, 7805-7831; Zenodo 15264835), five sit at an
# ambient ABL of 202-391 m and would have been discarded unclassified, though their fires grew
# observed fire-induced boundary layers of 451-2850 m. A gate that rejects the regime the
# instrument campaign was built to sample cannot be defended, so it is gone.
#
# Removing it does not resurrect the collapsed evening hours on its own: the mixed-layer
# gradient stays undefined where there is no mixed layer to fit, and those columns remain
# PYROCONV_NODATA. That is the honest outcome -- see docs/graf_vs_pyflam_2026-07-26.md.
ABL_MIN_M = 150.0
# Levels a least-squares mixed-layer dtheta/dz needs inside the layer. A hard floor, not a
# preference: below it :func:`_masked_lin_slope` returns nan, ``valid`` goes False and the cell
# drops out of the class map entirely -- which is why the *share of columns clearing it*, and
# not the typical level count, is what says whether the gradient is measured or unsupported.
ML_FIT_MIN_PTS = 3
# Sentinel for "this column could not be classified" -- no usable profile, a non-finite
# diagnostic, or an ABL below the mixing floor. Distinct from class 0 (surface plume), which
# is a *diagnosis*: the ladder ran and found no significant cloud development. Folding the two
# together paints unclassifiable ground as quiet, which is the opposite of an honest blank --
# most visibly around sunset, when the mixed layer collapses and the whole domain drops out.
PYROCONV_NODATA = -1
# Reference plume potential-temperature excess (K) for the dry-pyrocloud fireABL diagnostic:
# a fixed intense-fire forcing applied to every cell (the "potential"/upper-bound framing).
#
# Prescribed in kelvin, not derived from a heat flux. Until 2026-07-26 this was a 200 W/m^2
# reference flux converted by ``fire_parcel_theta_excess``, which (a) omitted the specific heat
# capacity and so returned J/kg rather than kelvin, overstating the scale ~88x, and (b) is a
# mixed-layer *turbulence* similarity scale even when dimensionally correct -- a few tenths of
# a kelvin, an order of magnitude below what a plume core actually is. No cell-averaged flux
# reconciles the two: it spreads the fire's heat over ground that is mostly not burning.
#
# 10 K is the intense end of the only direct measurements available: the GRAF in-plume campaign
# recorded 0.1-13.1 K (median 4.0) as the in-plume-minus-environment excess over the lowest
# 200 m, across 10 soundings at 7 wildfires (Castellnou Ribau et al. 2025, AMT 18, 7805-7831).
# The upper-bound framing takes the intense end, as the old flux comment intended.
#
# Treat the resulting decoupling ratio qualitatively. Even when forced with each fire's OWN
# measured excess, the encroachment model reproduces observed fireABL at ~151 % relative error
# on the 4 usable campaign observations (mean |error| 1156 m against observations of
# 371-2850 m). The forcing is now on measured ground; the model it feeds is not yet skilful.
_REFERENCE_THETA_EXCESS_K = 10.0

# Profile path: the levels the ICON-2I open-data archive publishes below 500 hPa.
PROFILE_LEVELS = (1000, 925, 850, 700, 500)

# Minimum usable levels for a classification. The reference implementation demands 4,
# but it is written for ERA5's 19 levels, where losing the underground ones is free.
# Here the archive gives 5: above ~750 m of terrain the 1000 *and* 925 hPa levels are
# below ground, leaving 850/700/500 -- so a threshold of 4 would silently blank the
# whole Apennine ridge, i.e. the highest-relief ground in the domain. Three levels plus
# the 2 m / 10 m surface state still bracket the ABL and the cap layer; the profile is
# coarser there, and the exported ABL raster is the place to see it.
MIN_PROFILE_LEVELS = 3


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


# --- profile path -------------------------------------------------------------

def _as_height_m(da):
    """Geometric height (m) from a DataArray that may hold a geopotential (m2/s2).

    Decided from the GRIB's declared units / standard name -- ``m**2 s**-2`` or
    ``geopotential`` means divide by g. Guessing from the magnitude does not work: at
    850 hPa the geopotential is ~15362 m2/s2, and a threshold set high enough to catch
    that would also "convert" a genuine 500 hPa height, while one set low enough to
    leave heights alone silently passes 15362 through as 15 km and throws the level
    away. The magnitude check survives only as a fallback for files with no units.
    """
    from pyflam.atmosphere import _G

    a = np.asarray(da.values, float)
    units = str(da.attrs.get("units", "")).replace(" ", "").lower()
    name = str(da.attrs.get("standard_name", "")).lower()
    if units in ("m**2s**-2", "m2s-2", "m^2s^-2", "j/kg", "jkg-1") or name == "geopotential":
        return a / _G
    if units in ("m", "metres", "meters", "gpm") or "height" in name:
        return a
    return a / _G if np.nanmedian(np.abs(a)) > 1.0e4 else a     # last resort


def read_icon2i_profile(files, bbox, run_dt, valid_dt, hours, levels=PROFILE_LEVELS):
    """Read the full ICON-2I profile (geopotential/T/RH/wind + surface) for an AOI.

    ``files`` is the dict returned by :func:`pyflam.atmosphere.fetch_icon2i_mistral`
    with ``fields=ICON2I_PROFILE_FIELDS``. Returns a dict whose per-level fields are
    stacked level-major with levels **ascending in height** (i.e. descending in
    pressure): ``z_asl``/``T``/``RH``/``U``/``V`` are ``(nlev, ntime, ny, nx)``, and
    ``hsurf``/``frland`` are static 2-D fields.

    ICON's FI is a *geopotential* (m2/s2); some archives serve a geopotential *height*
    (m) under the same short name. They differ by a factor g, and both are plausible
    magnitudes -- 850 hPa is 15362 m2/s2, which is also a believable height in metres
    if you are not careful -- so the units are read from the GRIB metadata rather than
    guessed from the magnitude. See :func:`_as_height_m`.
    """
    asc = sorted(levels, reverse=True)          # 1000 -> 500 = ascending in height
    z, T, RH, U, V = [], [], [], [], []
    ref = None
    for p in asc:
        dz, v = _open_grib_subset(files[f"FI{p}"], bbox)
        z.append(_as_height_m(dz[v]))
        dt_, v = _open_grib_subset(files[f"T{p}"], bbox); T.append(dt_[v].values)
        dr, v = _open_grib_subset(files[f"RELHUM{p}"], bbox); RH.append(dr[v].values)
        du, v = _open_grib_subset(files[f"U{p}"], bbox); U.append(du[v].values)
        dv, v = _open_grib_subset(files[f"V{p}"], bbox); V.append(dv[v].values)
        ref = dt_

    t2, v = _open_grib_subset(files["T2M"], bbox); T2m = t2[v].values
    td, v = _open_grib_subset(files["TD2M"], bbox); Td2m = td[v].values
    u1, v = _open_grib_subset(files["U10"], bbox); U10 = u1[v].values
    v1, vn = _open_grib_subset(files["V10"], bbox); V10 = v1[vn].values
    ps, v = _open_grib_subset(files["PS"], bbox); PS = ps[v].values
    hs, v = _open_grib_subset(files["HSURF"], bbox); HS = _as_height_m(hs[v])
    fl, v = _open_grib_subset(files["FRLAND"], bbox); FL = np.asarray(fl[v].values, float)

    lat = ref["latitude"].values
    lon = ref["longitude"].values
    sh = (ref["step"].values / np.timedelta64(1, "h")).astype(int)
    idx = []
    for h in hours:
        target = int((valid_dt.replace(hour=h) - run_dt).total_seconds() // 3600)
        m = np.where(sh == target)[0]
        if not m.size:
            raise ValueError(f"valid {h:02d}Z needs forecast step +{target} h, "
                             f"beyond this run (max +{int(sh.max())} h)")
        idx.append(int(m[0]))

    return dict(lat=lat, lon=lon, idx=idx, levels=tuple(asc),
                z_asl=np.stack(z), T=np.stack(T), RH=np.stack(RH),
                U=np.stack(U), V=np.stack(V),
                T2m=T2m, Td2m=Td2m, U10=U10, V10=V10, PS=PS,
                hsurf=HS, frland=FL)


def _edge_value(z, x, *, top: bool):
    """Value of ``x`` at the lowest (or highest) *usable* level of each column.

    Unusable levels are ``nan`` (underground, or missing), and which ones they are
    varies per cell -- a 1000 hPa level is below ground over the Apennines but not
    over the coast -- so the first/last valid index has to be found per column rather
    than assumed to be 0 / -1.
    """
    ok = np.isfinite(z) & np.isfinite(x)
    any_ok = ok.any(axis=0)
    order = ok[::-1] if top else ok
    idx = np.argmax(order, axis=0)
    if top:
        idx = (ok.shape[0] - 1) - idx
    val = np.take_along_axis(x, idx[None, ...], axis=0)[0]
    return np.where(any_ok, val, np.nan)


def _interp_at(z, x, ztarget):
    """Per-column linear interpolation of ``x(z)`` to a 2-D target height.

    ``z``/``x`` are ``(nlev, ny, nx)`` with ``z`` ascending along axis 0; ``ztarget``
    is ``(ny, nx)``. Levels that are unusable in a given column carry ``nan`` and are
    skipped, so a column whose bottom level is underground still interpolates from its
    remaining levels. Targets outside the usable span clamp to the nearest usable
    level. Vectorised over the grid (the level loop is short -- 5 levels -- so this
    stays cheap).
    """
    out = np.full(ztarget.shape, np.nan)
    for k in range(z.shape[0] - 1):
        z0, z1 = z[k], z[k + 1]
        x0, x1 = x[k], x[k + 1]
        span = z1 - z0
        good = (np.isfinite(z0) & np.isfinite(z1) & np.isfinite(x0) & np.isfinite(x1)
                & (np.abs(span) > 1e-6))
        f = np.divide(ztarget - z0, span, out=np.zeros_like(ztarget), where=good)
        seg = np.isnan(out) & good & (ztarget >= z0) & (ztarget <= z1)
        out = np.where(seg, x0 + np.clip(f, 0.0, 1.0) * (x1 - x0), out)

    z_lo = _edge_value(z, z, top=False); x_lo = _edge_value(z, x, top=False)
    z_hi = _edge_value(z, z, top=True); x_hi = _edge_value(z, x, top=True)
    out = np.where(np.isnan(out) & (ztarget <= z_lo), x_lo, out)
    out = np.where(np.isnan(out) & (ztarget >= z_hi), x_hi, out)
    return out


def _layer_gradient(z, x, zlo, zhi):
    """d<x>/dz (per metre) across a per-cell layer ``[zlo, zhi]`` (2-D bounds)."""
    dz = zhi - zlo
    lo = _interp_at(z, x, zlo)
    hi = _interp_at(z, x, zhi)
    return np.where(dz > 1.0, (hi - lo) / np.where(dz > 1.0, dz, 1.0), np.nan)


def _layer_mean(z, x, zlo, zhi, samples: int = 5):
    """Mean of ``x`` across a per-cell layer, from evenly spaced sample heights."""
    acc = np.zeros(zlo.shape)
    for i in range(samples):
        zi = zlo + (zhi - zlo) * (i / (samples - 1.0))
        acc += _interp_at(z, x, zi)
    return acc / samples


def _masked_lin_slope(z, x, mask, min_pts: int = ML_FIT_MIN_PTS):
    """Per-cell least-squares slope dx/dz over the masked levels (vectorised).

    ``z``/``x``/``mask`` are ``(nlev, ny, nx)``. Returns a 2-D slope, ``nan`` where a
    column has fewer than ``min_pts`` masked levels or zero spread in z. This is the
    honest mixed-layer gradient -- it needs enough *real* levels inside the layer to
    fit a line, which only the model-level (ICON-EU) path provides; on 5 pressure
    levels it falls back to nan and the caller uses a proxy instead.
    """
    w = (mask & np.isfinite(z) & np.isfinite(x)).astype(float)
    n = w.sum(0)
    safe = np.maximum(n, 1.0)
    zbar = (w * z).sum(0) / safe
    xbar = (w * x).sum(0) / safe
    cov = (w * (z - zbar) * (x - xbar)).sum(0)
    var = (w * (z - zbar) ** 2).sum(0)
    ok = (n >= min_pts) & (var > 1e-6)
    return np.where(ok, cov / np.where(var > 1e-6, var, 1.0), np.nan)


ML_METHODS = ("fit_in_ml", "surface_to_parcel", "surface_to_abl", "mid_layer")


def profile_diagnostics(d, si, *, thresholds=None, shear=True,
                        ml_method: str = "surface_to_parcel"):
    """The five column diagnostics for one time slice of :func:`read_icon2i_profile`.

    Returns a dict of 2-D fields: ``abl``, ``lcl``, ``lcl_ratio``, ``ml_grad``
    (mixed-layer dtheta/dz), ``gamma`` (free-troposphere cap), ``rh_top`` (RH at the
    ABL top), ``shear_dist`` (distance from the shear maximum to the ABL/LCL, in ABL
    units) and ``valid`` (enough usable levels to classify).

    Also returns ``parcel_ml``, the well-mixed-layer depth.

    ``ml_method`` selects the mixed-layer stability diagnostic.

    **Read this before trusting the ML-stability gate.** The mixed-layer dtheta/dz is
    *not measurable* from a 5-pressure-level archive. Validated against IGRA radiosondes
    (JJA 12Z, period of record, n=2835 soundings): a 5-level column holds only **0-2 model
    levels inside the mixed layer** -- one or none in 92% of coastal columns, three or more
    in just 8.5% of inland ones -- and the lowest of those sits in the superadiabatic
    surface layer. Every *direct* estimate consequently fails. Fitting a gradient across the
    in-ML levels scores a Youden J of ~0.00 (no skill at all); taking theta at the Rib ABL
    top instead folds in the entrainment jump Delta-theta, which Castellnou et al. (2022,
    sec.2.1.1) treat as a state variable *separate* from the gradient the ladder conditions
    on (mixed-layer slab model; Vila-Guerau de Arellano et al. 2015).

    The default is therefore an honest **proxy**, not a measurement:

    * ``"surface_to_parcel"`` (**default**) -- bulk theta difference from 2 m to the parcel
      mixing depth. Since the parcel top is *defined* as theta_sfc + 0.5 K, this is
      identically ``0.5 / depth``: a **mixing-depth proxy** for ML stability (at the 1.1e-3
      threshold it says exactly "parcel mixing depth >= 455 m"), *not* an estimate of
      dtheta/dz. It is the default because it is the only candidate with any skill --
      against the radiosonde truth it scores J = 0.29 inland / 0.55 coastal (r = +0.50),
      versus J = 0.14 / 0.42 for ``surface_to_abl``. Deep mixing does genuinely imply an
      unstable ML; that correlation, not a gradient measurement, is what this rests on.
    * ``"surface_to_abl"`` -- bulk gradient to theta at the Rib ABL top. Superseded: the Rib
      height sits above the entrainment zone, so theta there is post-jump.
    * ``"mid_layer"`` -- the reference implementation's 0.2-0.8*ABL window. Superseded: at
      ~900 m level spacing it straddles the entrainment zone.

    Known bias of the default: it **over-flags**. Inland its specificity is only 0.29 --
    of the columns that are truly ML-stable it still calls ~71% pyroCu-capable (sensitivity
    0.99, so it rarely *misses* a capable column). Treat the ML-stability gate as a weak
    filter and do not read the resulting class *counts* quantitatively. Fixing this needs
    more vertical levels (ERA5's 37 pressure levels, or ICON native model levels), which
    this open-data archive does not publish.

    The Rib ABL remains the right height for the LCL/ABL and shear/ABL *ratios* -- that is
    what the paper uses it for; it is only the wrong upper bound for this diagnostic.

    The other layers follow the reference: the cap over ABL+200 m to ABL+1200 m, and
    the ABL-top RH as the mean over ABL +/- 150 m.

    ``shear`` is honoured only if the profile has enough levels to locate a shear
    maximum; with ICON-2I's 5-6 levels it does not, so ``shear_dist`` comes back all
    ``nan`` and the classifier drops to the four-diagnostic ladder. See
    :func:`pyflam.atmosphere.shear_height_adaptive` for why this is deliberate.
    """
    if ml_method not in ML_METHODS:
        raise ValueError(f"unknown ml_method {ml_method!r}; expected one of {ML_METHODS}")
    from pyflam.atmosphere import (
        _SHEAR_MIN_LEVELS, bulk_richardson_abl_grid, dewpoint_from_rh,
        lcl_height_bolton_m, relative_humidity_from_dewpoint,
        specific_humidity_from_rh, theta_kelvin, virtual_potential_temperature,
        DEFAULT_PYROCONV_THRESHOLDS, shear_height_window, parcel_mixing_depth_grid,
        fire_induced_abl_grid,
    )
    th = thresholds or DEFAULT_PYROCONV_THRESHOLDS

    levels = np.asarray(d["levels"], float)             # hPa, ascending in height
    z_agl = d["z_asl"][:, si] - d["hsurf"][None, ...]   # (nlev, ny, nx)
    T = d["T"][:, si]                                   # K
    RH = np.clip(d["RH"][:, si], 0.0, 100.0)            # %
    U, V = d["U"][:, si], d["V"][:, si]
    p_pa = levels[:, None, None] * 100.0

    T2m, Td2m = d["T2m"][si], d["Td2m"][si]             # K
    U10, V10 = d["U10"][si], d["V10"][si]
    ps_pa = d["PS"][si]
    if np.nanmedian(ps_pa) < 2000.0:                    # served in hPa
        ps_pa = ps_pa * 100.0

    # Drop levels below ground / above the useful column, exactly as the reference:
    # a 1000 hPa level is underground over the Apennines.
    ok = (np.isfinite(z_agl) & np.isfinite(T) & np.isfinite(RH)
          & (z_agl >= 20.0) & (z_agl <= th.max_abl_m + 3500.0)
          & (p_pa <= ps_pa[None, ...] + 100.0))
    valid = ok.sum(axis=0) >= MIN_PROFILE_LEVELS
    # Unusable levels become nan (not a sentinel height): every downstream routine is
    # nan-aware per column, so a cell whose 1000 hPa level is underground still uses
    # the levels it does have.
    z_use = np.where(ok, z_agl, np.nan)
    T_use = np.where(ok, T, np.nan)
    RH_use = np.where(ok, RH, np.nan)

    theta = theta_kelvin(T_use - 273.15, levels[:, None, None])
    q = specific_humidity_from_rh(RH_use, T_use, p_pa)
    theta_v = virtual_potential_temperature(theta, q)

    ps_hpa = ps_pa / 100.0
    rh_sfc = relative_humidity_from_dewpoint(T2m - 273.15, Td2m - 273.15)
    q_sfc = specific_humidity_from_rh(rh_sfc, T2m, ps_pa)
    theta_v_sfc = virtual_potential_temperature(theta_kelvin(T2m - 273.15, ps_hpa), q_sfc)

    abl = bulk_richardson_abl_grid(
        z_use, theta_v, U, V, theta_v_surface=theta_v_sfc,
        wind_u_surface=U10, wind_v_surface=V10)
    lcl = lcl_height_bolton_m(T2m, Td2m, ps_pa)

    # theta (dry, not virtual): the ladder's thresholds are defined on theta.
    theta_sfc = theta_kelvin(T2m - 273.15, ps_hpa)
    parcel_ml = parcel_mixing_depth_grid(z_use, theta, theta_sfc)

    def _bulk_to(top):
        return np.where(np.isfinite(top),
                        (_interp_at(z_use, theta, top) - theta_sfc)
                        / np.maximum(top - 2.0, 1.0), np.nan)

    if ml_method == "fit_in_ml":
        # A genuine least-squares dtheta/dz across the levels inside the mixed layer.
        # Needs many levels to have skill; on 5 pressure levels it returns nan almost
        # everywhere and should not be the default -- it is here for the model-level
        # (ICON-EU) path, which shares this helper via iconeu_diagnostics.
        in_ml = (z_use >= 80.0) & (z_use <= parcel_ml[None, ...])
        ml_grad = _masked_lin_slope(z_use, theta, in_ml)
    elif ml_method == "surface_to_parcel":
        ml_grad = _bulk_to(parcel_ml)          # below the entrainment jump
    elif ml_method == "surface_to_abl":
        ml_grad = _bulk_to(abl)                # includes the jump -- superseded
    else:
        abl_lo = np.maximum(100.0, 0.20 * abl)
        abl_hi = np.maximum(abl_lo + 150.0, 0.80 * abl)
        abl_hi = np.minimum(abl_hi, np.where(abl > 250.0, abl - 50.0, abl))
        ml_grad = _layer_gradient(z_use, theta, abl_lo, abl_hi)

    gamma = _layer_gradient(z_use, theta, abl + 200.0, abl + 1200.0)
    rh_top = _layer_mean(z_use, RH_use, np.maximum(50.0, abl - 150.0), abl + 150.0)

    # Dry-pyrocloud diagnostic (DIAGNOSTIC ONLY, no class label): the fire-induced
    # boundary layer a reference intense fire would grow by sensible heat alone, and
    # how far it decouples above the ambient ABL. This is the DRY counterpart to the
    # moist LCL/cap/shear ladder -- it fires on deep, hot, dry columns the ladder
    # scores low (Castellnou et al. 2022; Castellnou Ribau et al. 2024). Like the
    # "potential" panel it assumes a fire everywhere (here a fixed reference flux),
    # so it is an upper bound, not an expectation. The dry/moist split and any
    # in-plume-LCL offset are deliberately NOT applied (offset not supported by the
    # GRAF labels; see atmosphere.fire_induced_abl_grid).
    theta_mean_below = _layer_mean(z_use, theta,
                                   np.full(abl.shape, 50.0), np.maximum(abl, 200.0))
    fireabl = fire_induced_abl_grid(
        z_use, theta, theta_mean_below=theta_mean_below,
        theta_excess=np.full(abl.shape, _REFERENCE_THETA_EXCESS_K), blh=abl)
    with np.errstate(invalid="ignore", divide="ignore"):
        decoupling = fireabl / np.where(abl > 0, abl, np.nan)

    shear_dist = np.full(abl.shape, np.nan)
    n_levels = int(ok.sum(axis=0).max()) if ok.size else 0
    if shear and n_levels >= _SHEAR_MIN_LEVELS:
        ny, nx = abl.shape
        for a in range(ny):
            for b in range(nx):
                if not valid[a, b]:
                    continue
                zs = shear_height_window(z_use[:, a, b], U[:, a, b], V[:, a, b])
                if np.isfinite(zs) and abl[a, b] > 0:
                    shear_dist[a, b] = min(abs(zs - abl[a, b]),
                                           abs(zs - lcl[a, b])) / abl[a, b]

    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = lcl / np.where(abl > 0, abl, np.nan)
    return dict(abl=abl, lcl=lcl, lcl_ratio=ratio, ml_grad=ml_grad, gamma=gamma,
                rh_top=rh_top, shear_dist=shear_dist, valid=valid,
                parcel_ml=parcel_ml, fireabl=fireabl, decoupling=decoupling,
                n_levels=n_levels,
                # No least-squares fit on this path -- the ML gradient is a mixing-depth
                # proxy, so there is no fit support to report.
                ml_fit_support=None)


def classify_profile(diag, *, fli=None, fli_gate_kw=FLI_GATE_KW, ladder="adaptive",
                     thresholds=None, abl_min_m=ABL_MIN_M):
    """Per-cell pyroconvection class (0..4) from :func:`profile_diagnostics`.

    ``ladder`` is passed through to :func:`pyflam.atmosphere.pyroconvection_type`;
    the default ``"adaptive"`` runs the richest ladder the diagnostics support (five
    where the shear height resolved, else four), which on ICON-2I open data means the
    four-diagnostic ladder. ``fli`` (kW/m, optional) gates each cell to a surface
    plume below ``fli_gate_kw`` -- pass it for the *expected* map, omit it for the
    *potential* upper bound.

    ``abl_min_m`` is a plausibility floor, not a depth gate -- see :data:`ABL_MIN_M` for why
    the 600 m version was removed. Callers wanting a stricter floor pass their own.

    Cells that cannot be classified -- no usable profile, a non-finite diagnostic, or an ABL
    below ``abl_min_m`` -- come back as :data:`PYROCONV_NODATA`, **not** as class 0. The two
    say different things: class 0 is a diagnosis (the ladder ran and found a surface plume),
    nodata is its absence. They shared the value 0 until now, which drew unclassifiable ground
    in the "no significant convection" colour; around sunset the mixed layer collapses across
    the whole domain, so the map read as a forecast of a quiet evening when it was really an
    empty one. A fuel-gated cell below ``fli_gate_kw`` stays class 0 -- that *is* a diagnosis.

    Also returns the ladder actually used, so the product can state it.
    """
    from pyflam.atmosphere import (
        pyroconvection_type, PYROCONVECTION_TYPE_LEVEL, DEFAULT_PYROCONV_THRESHOLDS)
    th = thresholds or DEFAULT_PYROCONV_THRESHOLDS

    abl, ratio = diag["abl"], diag["lcl_ratio"]
    ml, gamma, rh_top = diag["ml_grad"], diag["gamma"], diag["rh_top"]
    sd, valid = diag["shear_dist"], diag["valid"]
    ny, nx = abl.shape
    out = np.full((ny, nx), PYROCONV_NODATA, np.int16)
    used = set()

    for a in range(ny):
        for b in range(nx):
            if not valid[a, b] or not np.isfinite(abl[a, b]) or abl[a, b] < abl_min_m:
                continue
            if not np.isfinite([ratio[a, b], ml[a, b], gamma[a, b]]).all():
                continue
            rh = float(rh_top[a, b]) if np.isfinite(rh_top[a, b]) else None
            sdv = float(sd[a, b]) if np.isfinite(sd[a, b]) else None
            fkw = None if fli is None else float(fli[a, b])
            eff = ladder
            if ladder == "adaptive":
                eff = "shear" if sdv is not None else "noshear" if rh is not None \
                    else "castellnou"
            used.add(eff)
            out[a, b] = PYROCONVECTION_TYPE_LEVEL[pyroconvection_type(
                lcl_abl_ratio=float(ratio[a, b]), ml_theta_gradient=float(ml[a, b]),
                gamma_theta=float(gamma[a, b]), ladder=eff, rh_top_abl=rh,
                shear_distance=sdv, fireline_intensity_kw=fkw,
                fli_threshold_kw=fli_gate_kw, thresholds=th)]
    return out, sorted(used)


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


# --- ICON-EU model-level path (the hybrid product's atmosphere) ----------------
#
# ICON-2I open data publishes 5 pressure levels, on which the mixed-layer dtheta/dz is
# not measurable (see profile_diagnostics and scripts/validation/). ICON-EU is coarser
# horizontally (6.5 km) but publishes the native model levels -- ~10 inside the mixed
# layer over Tuscany -- so the gradient becomes a real least-squares fit and the Rib ABL
# bias drops from ~-700 m to ~-70 m against radiosondes. The hybrid product takes the
# profile diagnostics from here and keeps ICON-2I's 2.2 km surface fields for the fuel
# gate, where fine terrain actually matters.

def _open_iconeu(path, bbox):
    """Open one decompressed ICON-EU GRIB, subset to ``bbox`` (n, w, s, e). -> (da)."""
    import xarray as xr
    n, w, s, e = bbox
    ds = xr.open_dataset(path, engine="cfgrib", backend_kwargs={"indexpath": ""})
    v = list(ds.data_vars)[0]
    la = ds["latitude"].values
    d = (ds.sel(latitude=slice(n, s), longitude=slice(w, e)) if la[0] > la[-1]
         else ds.sel(latitude=slice(s, n), longitude=slice(w, e)))
    return d[v]


def read_icon_eu(files, bbox, levels):
    """Read an ICON-EU model-level column stack for an AOI.

    ``files`` is the dict from :func:`pyflam.atmosphere.fetch_icon_eu`; ``bbox`` is
    ``(north, west, south, east)``; ``levels`` the full model levels (ascending in
    height, e.g. 74..51). Full-level heights AGL are reconstructed from the HHL half
    levels (full level k = mean of half levels k and k+1, minus the surface height).
    Returns a dict with 1-D ``lat``/``lon``, level-major ``(nlev, ny, nx)`` stacks
    ``z``/``T``/``QV``/``P``/``U``/``V``, and the 2-D surface fields.
    """
    n, w, s, e = bbox
    half = {k: np.asarray(_open_iconeu(files[f"HHL{k}"], bbox).values, float)
            for k in tuple(levels) + (max(levels) + 1,)}
    orog = half[max(levels) + 1]                        # HHL at the surface half level
    z = np.stack([0.5 * (half[k] + half[k + 1]) - orog for k in levels])

    def stack(var):
        return np.stack([np.asarray(_open_iconeu(files[f"{var}{k}"], bbox).values, float)
                         for k in levels])

    ref = _open_iconeu(files["T_2M"], bbox)
    return dict(
        lat=ref["latitude"].values, lon=ref["longitude"].values, levels=tuple(levels),
        z=z, T=stack("T"), QV=stack("QV"), P=stack("P"), U=stack("U"), V=stack("V"),
        T2m=np.asarray(_open_iconeu(files["T_2M"], bbox).values, float),
        Td2m=np.asarray(_open_iconeu(files["TD_2M"], bbox).values, float),
        PS=np.asarray(_open_iconeu(files["PS"], bbox).values, float),
        U10=np.asarray(_open_iconeu(files["U_10M"], bbox).values, float),
        V10=np.asarray(_open_iconeu(files["V_10M"], bbox).values, float),
        frland=np.asarray(_open_iconeu(files["FR_LAND"], bbox).values, float),
        orog=orog)


def iconeu_diagnostics(d, *, thresholds=None, ml_method="fit_in_ml", shear=True,
                       abl_method="rib"):
    """Profile diagnostics from an ICON-EU model-level stack (:func:`read_icon_eu`).

    Same output dict as :func:`profile_diagnostics` (``abl``, ``parcel_ml``, ``lcl``,
    ``lcl_ratio``, ``ml_grad``, ``gamma``, ``rh_top``, ``shear_dist``, ``valid``,
    ``n_levels``), so :func:`classify_profile` consumes it unchanged. ``n_levels`` is the
    **land-only** median count of model levels inside the mixed layer, and
    ``ml_fit_support`` the share of land columns holding at least ``ML_FIT_MIN_PTS`` of
    them -- the number that says whether the gradient is measured, since columns below the
    floor return a nan slope and leave the class map. The default
    ``ml_method="fit_in_ml"`` is the genuine least-squares mixed-layer gradient, which is
    only trustworthy *because* this is the model-level path (many levels inside the ML);
    the pressure-level path cannot use it (see :func:`profile_diagnostics`).
    """
    from pyflam.atmosphere import (
        theta_kelvin, specific_humidity_from_rh, virtual_potential_temperature,
        saturation_vapour_pressure_pa, bulk_richardson_abl_grid, parcel_mixing_depth_grid,
        lcl_height_bolton_m, relative_humidity_from_dewpoint, DEFAULT_PYROCONV_THRESHOLDS,
        shear_height_grid, shear_distance_grid, _EPSILON, max_rh_abl_grid,
        entrainment_jump_grid, fire_cape_grid, residual_layer_grid,
        pyrocb_firepower_threshold_grid)
    th = thresholds or DEFAULT_PYROCONV_THRESHOLDS

    z, T, QV, P, U, V = d["z"], d["T"], d["QV"], d["P"], d["U"], d["V"]
    T2m, Td2m, ps = d["T2m"], d["Td2m"], d["PS"]
    U10, V10 = d["U10"], d["V10"]

    theta = theta_kelvin(T - 273.15, P / 100.0)
    thv = virtual_potential_temperature(theta, QV)
    e = QV * P / (_EPSILON + (1.0 - _EPSILON) * QV)
    RH = np.clip(100.0 * e / saturation_vapour_pressure_pa(T), 1.0, 100.0)

    rh_s = relative_humidity_from_dewpoint(T2m - 273.15, Td2m - 273.15)
    q_s = specific_humidity_from_rh(rh_s, T2m, ps)
    theta_sfc = theta_kelvin(T2m - 273.15, ps / 100.0)
    thv_s = virtual_potential_temperature(theta_sfc, q_s)

    abl_rib = bulk_richardson_abl_grid(z, thv, U, V, theta_v_surface=thv_s,
                                       wind_u_surface=U10, wind_v_surface=V10)
    # ABL depth. Default "rib": the bulk-Richardson depth.
    #
    # "maxrh_floored" implements the criterion the source method uses -- Castellnou Ribau et al.
    # (2025) sec. 2.6 take the boundary-layer top as the height of maximum relative humidity,
    # "supplemented with numerical calculations using the bulk Richardson number" -- with that
    # supplement applied as a floor, since a moisture maximum below the dynamically diagnosed
    # mixing top is not a capping inversion. It was briefly the default and is **not**, because
    # measured on the ladder this product actually runs it is worse.
    #
    # Against 26 ambient campaign sondes at GRAF-labelled fires, on the **adaptive** ladder
    # (which resolved to the 5-diagnostic "shear" path in 51 of 52 cases):
    #
    #     rib             11/26 exact, 20/26 within one class, mean bias +0.12
    #     maxrh_floored    7/26 exact, 19/26 within one class, mean bias +0.58
    #
    # The ordering reverses on the 3-diagnostic "castellnou" ladder (+1.38 for maxrh_floored
    # against +1.81 for rib), which is what an earlier analysis measured and got wrong. The
    # reason is that a deeper ABL shifts both the cap layer (abl+200 to abl+1200) and the
    # RH-top window upward, and on the adaptive ladder those two gates carry the
    # discrimination -- improving the LCL/ABL ratio while degrading the gates is a net loss.
    # The castellnou ladder has no RH-top gate, so it cannot show this.
    #
    # Keep "maxrh_floored" available: it is the source method's own criterion, and on a ladder
    # that does not lean on the RH-top gate it is the better depth. See
    # docs/graf_vs_pyflam_2026-07-26.md sec. 17.
    if abl_method == "maxrh_floored":
        z_rh = max_rh_abl_grid(z, RH)
        abl = np.where(np.isfinite(z_rh) & (z_rh >= abl_rib), z_rh, abl_rib)
    elif abl_method == "rib":
        abl = abl_rib
    else:
        raise ValueError(f"unknown abl_method {abl_method!r}; expected 'maxrh_floored' or 'rib'")
    parcel_ml = parcel_mixing_depth_grid(z, theta, theta_sfc)
    lcl = lcl_height_bolton_m(T2m, Td2m, ps)

    # Item 5 -- the residual layer. By day this equals parcel_ml (the theta minimum is at the
    # surface); after the CBL decays it is the depth of the layer that retains the day's
    # near-neutral profile, while parcel_ml collapses to the shallow nocturnal stable layer.
    # Fitting the mixed-layer gradient inside a collapsed parcel depth has no levels to work
    # with and returns nan, which is what left the evening hours unclassifiable. Fitting it
    # over the residual layer is well-posed around the clock.
    resid_ml = residual_layer_grid(z, theta)
    fit_depth = np.where(np.isfinite(resid_ml), np.maximum(parcel_ml, resid_ml), parcel_ml)

    in_ml = (z >= 80.0) & (z <= fit_depth[None, ...])
    n_in_ml = in_ml.sum(0)
    # Summarise the fit's support over LAND only. The summer Tyrrhenian carries a shallow,
    # stably stratified marine layer -- a different regime, and one the fuel gate never
    # classifies -- so a whole-grid median understates the resolution on the ground the
    # product is actually about (today: land 14 levels vs sea 6 at 12Z). ``ml_fit_support``
    # is the share of land columns clearing ML_FIT_MIN_PTS; it, not the median, is what
    # licenses calling the gradient measured, since the columns below the floor are dropped.
    frl = d.get("frland")
    land = (np.asarray(frl, float) >= 0.5) if frl is not None else np.ones(n_in_ml.shape, bool)
    if not land.any():
        land = np.ones(n_in_ml.shape, bool)
    n_land = n_in_ml[land]
    n_levels_land = int(np.median(n_land))
    ml_fit_support = float(np.mean(n_land >= ML_FIT_MIN_PTS))
    if ml_method == "fit_in_ml":
        ml_grad = _masked_lin_slope(z, theta, in_ml)
    elif ml_method == "surface_to_parcel":
        ml_grad = np.where(np.isfinite(parcel_ml),
                           (_interp_at(z, theta, parcel_ml) - theta_sfc)
                           / np.maximum(parcel_ml - 2.0, 1.0), np.nan)
    else:
        raise ValueError(f"iconeu_diagnostics: unsupported ml_method {ml_method!r}")

    gamma = _layer_gradient(z, theta, abl + 200.0, abl + 1200.0)
    rh_top = _layer_mean(z, RH, np.maximum(50.0, abl - 150.0), abl + 150.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = lcl / np.where(abl > 0, abl, np.nan)

    # The fifth diagnostic: the model levels resolve a shear-maximum height (the coarse
    # pressure-level path cannot), so classify_profile(ladder="adaptive") runs the full
    # 5-diagnostic ladder where shear_dist is finite.
    if shear:
        shear_dist = shear_distance_grid(shear_height_grid(z, U, V), abl, lcl)
    else:
        shear_dist = np.full(abl.shape, np.nan)

    valid = np.isfinite(abl) & np.isfinite(ratio) & np.isfinite(ml_grad)

    # Dry-pyrocloud fireABL diagnostic (see profile_diagnostics for the rationale):
    # the sensible-heat-only decoupling of a reference intense fire. Model levels give
    # the real theta(z) stack, so this is the higher-fidelity path for the diagnostic.
    from pyflam.atmosphere import fire_induced_abl_grid
    theta_mean_below = _layer_mean(z, theta,
                                   np.full(abl.shape, 50.0), np.maximum(abl, 200.0))
    theta_excess = np.full(abl.shape, _REFERENCE_THETA_EXCESS_K)
    fireabl = fire_induced_abl_grid(
        z, theta, theta_mean_below=theta_mean_below,
        theta_excess=theta_excess, blh=abl)
    with np.errstate(invalid="ignore", divide="ignore"):
        decoupling = fireabl / np.where(abl > 0, abl, np.nan)

    # Items 3-4: the two variables the source method computes and pyflam did not. Both are
    # DIAGNOSTIC for now -- exposed and exported, but not yet wired into the ladder, so this
    # commit changes no class. See docs/graf_vs_pyflam_2026-07-26.md for why they are the
    # candidates for the over-extent, and why gating on them is a separate decision.
    delta_theta = entrainment_jump_grid(z, theta, abl)
    firecape = fire_cape_grid(z, thv, theta_excess=theta_excess)
    # Can the reference fire's parcel clear the capping jump at all? This is the penetration
    # test of Castellnou et al. (2022) sec.2.1.2, expressed as a ratio: >= 1 means the fire's
    # theta excess exceeds the entrainment-zone jump it has to cross.
    with np.errstate(invalid="ignore", divide="ignore"):
        penetration = theta_excess / np.where(delta_theta > 0, delta_theta, np.nan)

    # PyroCb Firepower Threshold (Tory & Kepert 2021 eq 31) -- the minimum TOTAL firepower
    # this column needs for pyroCb, in GW. DIAGNOSTIC: exported, never gated on. The authors
    # state the absolute values are unverified and recommend the field for *relative* threat,
    # and comparing it to a fire needs a total power (Byram intensity x head-fire length, or
    # their appendix D area-burned form), which this product does not compute per cell.
    #
    # Needs the profile up to the -20 C electrification level, which is why
    # ICON_EU_MODEL_LEVELS reaches ~8.9 km; on a shallower stack every column returns nan.
    pft = pyrocb_firepower_threshold_grid(z, T, d["P"], d["QV"], U, V,
                                          surface_pressure_pa=ps)

    return dict(abl=abl, abl_rib=abl_rib, lcl=lcl, lcl_ratio=ratio, ml_grad=ml_grad,
                gamma=gamma,
                rh_top=rh_top, shear_dist=shear_dist, valid=valid,
                parcel_ml=parcel_ml, residual_ml=resid_ml, fireabl=fireabl,
                decoupling=decoupling, delta_theta=delta_theta, firecape=firecape,
                penetration=penetration,
                pft_gw=pft["pft_gw"], z_fc=pft["z_fc_m"],
                delta_theta_fc=pft["delta_theta_fc_k"], u_ml=pft["u_ml_ms"],
                n_levels=n_levels_land, ml_fit_support=ml_fit_support)


def regrid_to(field, lat_src, lon_src, lat_dst, lon_dst):
    """Bilinear-interpolate a 2-D field from one lat/lon grid onto another.

    Used to move ICON-EU diagnostics (6.5 km) onto the ICON-2I 2.2 km grid the hybrid
    product renders on. Points outside the source grid extrapolate from the edge rather
    than go blank (the ICON-EU domain amply covers Tuscany, so this only bites at the
    very border).
    """
    from scipy.interpolate import RegularGridInterpolator

    la, fld = lat_src, np.asarray(field, float)
    if la[0] > la[-1]:
        la = la[::-1]
        fld = fld[::-1]
    f = RegularGridInterpolator((la, lon_src), fld, bounds_error=False, fill_value=None)
    LO, LA = np.meshgrid(lon_dst, lat_dst)
    return f(np.stack([LA.ravel(), LO.ravel()], -1)).reshape(LA.shape)


def regrid_diagnostics(diag, lat_src, lon_src, lat_dst, lon_dst, *, abl_min_m=ABL_MIN_M):
    """Move a whole diagnostics dict onto a target grid and rebuild ``valid``.

    Interpolates each 2-D diagnostic field (``shear_dist`` included, so the regridded product
    keeps the 5-diagnostic ladder where the model levels resolved a shear maximum), then
    recomputes ``valid`` on the target grid (finite ABL/ratio/gradient, ABL above the floor)
    rather than interpolating a boolean. ``n_levels`` is carried through unchanged.
    """
    fields = ("abl", "lcl", "lcl_ratio", "ml_grad", "gamma", "rh_top", "parcel_ml",
              "shear_dist", "fireabl", "decoupling", "residual_ml", "delta_theta",
              "firecape", "penetration", "pft_gw", "z_fc", "delta_theta_fc", "u_ml",
              "abl_rib")
    out = {k: regrid_to(diag[k], lat_src, lon_src, lat_dst, lon_dst)
           for k in fields if k in diag}
    out["valid"] = (np.isfinite(out["abl"]) & np.isfinite(out["lcl_ratio"])
                    & np.isfinite(out["ml_grad"]) & (out["abl"] >= abl_min_m))
    out["n_levels"] = diag["n_levels"]
    out["ml_fit_support"] = diag.get("ml_fit_support")
    return out
