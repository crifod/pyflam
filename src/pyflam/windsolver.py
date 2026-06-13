"""Native mass-consistent (diagnostic) wind model.

A pure-Python reimplementation of the physics behind WindNinja's *default*
solver: the Sasaki variational "mass-consistent" method. It takes a single input
wind plus a DEM and returns a terrain-aware gridded wind field that conserves
mass (incompressible continuity), with no external binary.

Method (Sasaki 1970; Sherman 1978; Ross et al. 1988; Forthofer et al. 2014):

1. Build a first-guess 3D wind field ``u0`` from the input wind with a log-law
   vertical profile, draped over the terrain. Over slopes this guess is *not*
   divergence-free.
2. Find the field ``u`` closest to ``u0`` that satisfies ``div(u) = 0``. With a
   Lagrange multiplier ``lambda`` this reduces to an anisotropic Poisson problem

       lambda_xx + lambda_yy + TR * lambda_zz = -2 div(u0)

   solved here by finite volumes on a terrain-masked Cartesian grid, then

       u = u0 + 1/2 grad_h(lambda),   w = w0 + TR/2 lambda_z.

   ``TR = (alpha1/alpha2)^2`` is the stability ratio (1.0 = neutral).

Boundary conditions: no-flux (Neumann) at the terrain surface and bottom; open
(Dirichlet ``lambda = 0``) on the lateral faces and the domain top.

This is the most physically grounded — and most code — module in pyflam. It is a
*structured, staircase-terrain* discretization; WindNinja uses a smoother
terrain-following FEM mesh, so expect qualitative agreement (ridge speed-up,
valley channeling) rather than cell-exact equality. Validate against WindNinja's
``*_vel``/``*_ang`` grids; the solver self-checks that its output is discretely
divergence-free (see :func:`solve_mass_consistent` ``return_diagnostics``).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import cg, spsolve

from .wind import WindField

VON_KARMAN = 0.41


# Rough surface-roughness length z0 (m) by fuel-model family. Approximate values
# in the spirit of WindNinja's grass/brush/trees roughness classes.
def roughness_from_fuel(fuel_model: np.ndarray) -> np.ndarray:
    """Map fuel-model numbers to a roughness length z0 (m), per cell."""
    fm = np.asarray(fuel_model)
    z0 = np.full(fm.shape, 0.1, dtype=float)  # default / mixed
    def _in(lo, hi):
        return (fm >= lo) & (fm <= hi)
    z0[_in(1, 3)] = 0.03        # grass (FM1-3)
    z0[_in(101, 109)] = 0.03    # GR
    z0[_in(4, 7)] = 0.30        # brush/chaparral (FM4-7)
    z0[(_in(121, 124)) | (_in(141, 149))] = 0.30  # GS, SH
    z0[_in(8, 13)] = 1.0        # timber/slash litter (FM8-13)
    z0[(_in(161, 165)) | (_in(181, 189)) | (_in(201, 204))] = 1.0  # TU, TL, SB
    return z0


@dataclass
class _Grid:
    """A terrain-masked Cartesian grid and its cell bookkeeping."""

    fluid: np.ndarray   # (nz, ny, nx) bool
    idmap: np.ndarray   # (nz, ny, nx) int, cell index or -1
    z_centers: np.ndarray  # (nz,) absolute heights of cell centres
    elev: np.ndarray    # (ny, nx)
    h: float            # horizontal cell size (m)
    dz: float
    n: int              # number of fluid cells

    @property
    def coords(self):
        return np.where(self.fluid)  # (kk, jj, ii) in index order


def _build_grid(elev, h, *, dz, top_margin, max_layers):
    z_floor = float(elev.min())
    z_top = float(elev.max()) + top_margin
    nz = int(math.ceil((z_top - z_floor) / dz))
    if nz > max_layers:
        dz = (z_top - z_floor) / max_layers
        nz = max_layers
    z_centers = z_floor + (np.arange(nz) + 0.5) * dz
    fluid = z_centers[:, None, None] >= elev[None, :, :]
    idmap = np.full(fluid.shape, -1, dtype=np.int64)
    idmap[fluid] = np.arange(int(fluid.sum()))
    return _Grid(fluid=fluid, idmap=idmap, z_centers=z_centers, elev=elev,
                 h=h, dz=dz, n=int(fluid.sum()))


def _log_profile(z_agl, z0, speed_ref, z_ref):
    """Log-law wind speed at height ``z_agl`` (AGL), matched to ``speed_ref``."""
    num = np.log((np.maximum(z_agl, 0.0) + z0) / z0)
    den = math.log((z_ref + float(np.mean(z0))) / float(np.mean(z0)))
    return speed_ref * num / den


def _first_guess(grid: _Grid, speed, direction_from, z0, z_ref):
    """Horizontal first-guess wind (east, north components) on the 3D grid."""
    phi = math.radians(direction_from)
    comp_e = -speed * math.sin(phi)   # wind blows *toward* (from + 180)
    comp_n = -speed * math.cos(phi)
    z_agl = grid.z_centers[:, None, None] - grid.elev[None, :, :]
    z0b = np.broadcast_to(z0, grid.elev.shape)[None, :, :]
    speed_ref_unit = _log_profile(z_agl, z0b, 1.0, z_ref)  # profile shape only
    u0 = comp_e * speed_ref_unit
    v0 = comp_n * speed_ref_unit
    return u0 * grid.fluid, v0 * grid.fluid


def _assemble(grid: _Grid, u0, v0, tr):
    """Build the sparse Poisson system A lambda = b for the mass correction."""
    kk, jj, ii = grid.coords
    nz, ny, nx = grid.fluid.shape
    h, dz = grid.h, grid.dz
    a_h = dz                      # horizontal face coefficient
    a_z = tr * h * h / dz         # vertical face coefficient
    area_h = h * dz
    n = grid.n

    diag = np.zeros(n)
    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    data: list[np.ndarray] = []
    b = np.zeros(n)

    # axis: 0=z, 1=y(north), 2=x(east); sign +/-1; coefficient; velocity comp.
    dirs = [
        (0, +1, a_z, None), (0, -1, a_z, None),     # vertical (w0 = 0)
        (1, +1, a_h, v0), (1, -1, a_h, v0),         # north-south
        (2, +1, a_h, u0), (2, -1, a_h, u0),         # east-west
    ]
    comp_self_cache = {id(v0): v0[kk, jj, ii], id(u0): u0[kk, jj, ii]}

    for axis, s, a, comp in dirs:
        nk, nj, ni = kk.copy(), jj.copy(), ii.copy()
        if axis == 0:
            nk = kk + s
        elif axis == 1:
            nj = jj + s
        else:
            ni = ii + s

        in_dom = ((nk >= 0) & (nk < nz) & (nj >= 0) & (nj < ny)
                  & (ni >= 0) & (ni < nx))
        nid = np.full(n, -2, dtype=np.int64)  # -2 = out of domain
        nid[in_dom] = grid.idmap[nk[in_dom], nj[in_dom], ni[in_dom]]
        fluid_nb = nid >= 0
        # Out-of-domain classification: bottom (axis z, going down) is ground
        # (Neumann); all other out-of-domain faces are open (Dirichlet 0).
        out = ~in_dom
        ground = out & (axis == 0) & (s < 0)
        open_nb = out & ~ground

        # Matrix: fluid neighbour -> off-diagonal + diagonal; open -> diagonal.
        diag[fluid_nb] -= a
        rows.append(np.where(fluid_nb)[0])
        cols.append(nid[fluid_nb])
        data.append(np.full(int(fluid_nb.sum()), a))
        diag[open_nb] -= a

        # RHS: outward first-guess flux through this face (walls -> 0).
        if comp is not None:
            comp_self = comp_self_cache[id(comp)]
            face_vel = np.zeros(n)
            face_vel[fluid_nb] = 0.5 * (comp_self[fluid_nb]
                                        + comp[nk[fluid_nb], nj[fluid_nb], ni[fluid_nb]])
            face_vel[open_nb] = comp_self[open_nb]
            b += -2.0 * s * face_vel * area_h

    rows.append(np.arange(n))
    cols.append(np.arange(n))
    data.append(diag)
    a_mat = sp.coo_matrix((np.concatenate(data),
                           (np.concatenate(rows), np.concatenate(cols))),
                          shape=(n, n)).tocsr()
    return a_mat, b


def _neighbor_lambda(grid: _Grid, lam_field, axis, s):
    """Neighbour lambda values with BCs: solid -> self (Neumann), open -> 0."""
    kk, jj, ii = grid.coords
    nz, ny, nx = grid.fluid.shape
    nk, nj, ni = kk.copy(), jj.copy(), ii.copy()
    if axis == 0:
        nk = kk + s
    elif axis == 1:
        nj = jj + s
    else:
        ni = ii + s
    self_lam = lam_field[kk, jj, ii]
    out = ((nk < 0) | (nk >= nz) | (nj < 0) | (nj >= ny) | (ni < 0) | (ni >= nx))
    res = np.empty(grid.n)
    inb = ~out
    vals = np.where(grid.fluid[nk[inb], nj[inb], ni[inb]],
                    lam_field[nk[inb], nj[inb], ni[inb]],
                    self_lam[inb])  # solid -> self (Neumann)
    res[inb] = vals
    ground = out & (axis == 0) & (s < 0)
    res[out & ~ground] = 0.0          # open -> Dirichlet 0
    res[ground] = self_lam[ground]    # ground -> Neumann
    return res


def solve_mass_consistent(
    elevation: np.ndarray,
    cellsize: float,
    *,
    speed: float,
    direction: float,
    z0: float | np.ndarray = 0.1,
    reference_height: float = 6.1,
    output_height: float = 6.1,
    stability_ratio: float = 1.0,
    dz: float | None = None,
    top_margin: float | None = None,
    max_layers: int = 40,
    west: float = 0.0,
    north: float = 0.0,
    speed_units: str = "m/s",
    crs: object | None = None,
    return_diagnostics: bool = False,
):
    """Solve the mass-consistent wind field over a DEM.

    Parameters
    ----------
    elevation:
        2D DEM in metres (row 0 = north).
    cellsize:
        Horizontal cell size in metres.
    speed, direction:
        Input wind speed (in m/s) and direction it blows *from* (deg from north).
    z0:
        Surface roughness length (m), scalar or per-cell (see
        :func:`roughness_from_fuel`).
    reference_height, output_height:
        Heights AGL (m) of the input wind and of the returned wind.
    stability_ratio:
        ``TR`` in the Poisson equation; 1.0 neutral, <1 more stable (vertical
        adjustment suppressed), >1 less stable.

    Returns a :class:`~pyflam.wind.WindField` (and a diagnostics dict if
    ``return_diagnostics``).
    """
    elev = np.asarray(elevation, dtype=float)
    relief = float(elev.max() - elev.min())
    if top_margin is None:
        top_margin = max(3.0 * relief, 200.0)
    if dz is None:
        dz = cellsize
    z0_arr = np.broadcast_to(np.asarray(z0, dtype=float), elev.shape)

    grid = _build_grid(elev, cellsize, dz=dz, top_margin=top_margin,
                       max_layers=max_layers)
    u0, v0 = _first_guess(grid, speed, direction, z0_arr, reference_height)

    a_mat, b = _assemble(grid, u0, v0, stability_ratio)
    # A is symmetric negative-definite (open boundaries make it nonsingular);
    # solve the SPD system (-A) lambda = (-b).
    neg_a = (-a_mat).tocsr()
    if grid.n <= 4000:
        lam = spsolve(neg_a, -b)
    else:
        lam, info = cg(neg_a, -b, rtol=1e-8, maxiter=5000)
        if info != 0:  # pragma: no cover - fallback for hard cases
            lam = spsolve(neg_a, -b)

    lam_field = np.zeros(grid.fluid.shape)
    kk, jj, ii = grid.coords
    lam_field[kk, jj, ii] = lam

    # Recover corrected velocities at cell centres.
    h = grid.h
    dlam_dx = (_neighbor_lambda(grid, lam_field, 2, +1)
               - _neighbor_lambda(grid, lam_field, 2, -1)) / (2 * h)
    dlam_dy = (_neighbor_lambda(grid, lam_field, 1, +1)
               - _neighbor_lambda(grid, lam_field, 1, -1)) / (2 * h)
    u = u0[kk, jj, ii] + 0.5 * dlam_dx
    v = v0[kk, jj, ii] + 0.5 * dlam_dy

    u_field = np.zeros(grid.fluid.shape)
    v_field = np.zeros(grid.fluid.shape)
    u_field[kk, jj, ii] = u
    v_field[kk, jj, ii] = v

    speed2d, dir2d = _extract_at_height(
        grid, u_field, v_field, output_height, z0_arr)

    wf = WindField(
        speed=speed2d, direction=dir2d, cellsize=cellsize,
        west=west, north=north, speed_units=speed_units,
        height=output_height, height_units="m", crs=crs,
    )
    if not return_diagnostics:
        return wf

    diag = {
        "max_divergence": _max_flux_divergence(
            grid, u0, v0, lam_field, stability_ratio),
        "n_cells": grid.n,
        "n_layers": grid.fluid.shape[0],
        "dz": grid.dz,
        "lambda_field": lam_field,
    }
    return wf, diag


def _extract_at_height(grid: _Grid, u_field, v_field, output_height, z0):
    """Sample horizontal wind at ``output_height`` AGL for every column.

    Above the lowest cell centre we interpolate linearly in height; below it we
    extrapolate with the log law (the surface layer), which is how a near-surface
    target like 6.1 m (20 ft) is reached from a coarser vertical grid.
    """
    ny, nx = grid.elev.shape
    speed2d = np.zeros((ny, nx))
    dir2d = np.zeros((ny, nx))
    for j in range(ny):
        for i in range(nx):
            ks = np.where(grid.fluid[:, j, i])[0]
            if ks.size == 0:
                continue
            z_agl = grid.z_centers[ks] - grid.elev[j, i]
            ue_col, vn_col = u_field[ks, j, i], v_field[ks, j, i]
            if output_height < z_agl[0]:
                z0c = float(z0[j, i])
                factor = (math.log((output_height + z0c) / z0c)
                          / math.log((z_agl[0] + z0c) / z0c))
                ue, vn = ue_col[0] * factor, vn_col[0] * factor
            else:
                ue = np.interp(output_height, z_agl, ue_col)
                vn = np.interp(output_height, z_agl, vn_col)
            speed2d[j, i] = math.hypot(ue, vn)
            dir2d[j, i] = math.degrees(math.atan2(-ue, -vn)) % 360.0
    return speed2d, dir2d


def _max_flux_divergence(grid: _Grid, u0, v0, lam_field, tr):
    """Largest |3D face-flux divergence| of the corrected field (must be ~0).

    Sums, over all six faces of every fluid cell, the outward normal flux of the
    corrected wind ``u0 + 1/2 grad(lambda)`` (vertical correction scaled by the
    stability ratio ``tr``). This is the true mass-consistency check: it should
    be machine-zero where ``A lambda = b`` was solved exactly.
    """
    kk, jj, ii = grid.coords
    nz, ny, nx = grid.fluid.shape
    h, dz = grid.h, grid.dz
    area_h, area_z = h * dz, h * h
    lam_c = lam_field[kk, jj, ii]
    total = np.zeros(grid.n)

    # (axis, sign, advected component, face area, correction coefficient).
    coef_h = 0.5 * area_h / h      # outward correction flux per unit (lam_n-lam_c)
    coef_z = 0.5 * tr * area_z / dz
    z0 = np.zeros(grid.n)
    specs = [
        (0, +1, z0, area_z, coef_z), (0, -1, z0, area_z, coef_z),
        (1, +1, v0[kk, jj, ii], area_h, coef_h), (1, -1, v0[kk, jj, ii], area_h, coef_h),
        (2, +1, u0[kk, jj, ii], area_h, coef_h), (2, -1, u0[kk, jj, ii], area_h, coef_h),
    ]
    raw = {0: None, 1: v0, 2: u0}
    for axis, s, comp_self, area, coef in specs:
        nk, nj, ni = kk.copy(), jj.copy(), ii.copy()
        if axis == 0:
            nk = kk + s
        elif axis == 1:
            nj = jj + s
        else:
            ni = ii + s
        out = ((nk < 0) | (nk >= nz) | (nj < 0) | (nj >= ny)
               | (ni < 0) | (ni >= nx))
        inb = ~out
        is_fluid = np.zeros(grid.n, bool)
        is_fluid[inb] = grid.fluid[nk[inb], nj[inb], ni[inb]]
        ff = inb & is_fluid
        ground = out & (axis == 0) & (s < 0)
        open_nb = out & ~ground

        # Advective (first-guess) outward flux: needs the outward sign s.
        adv = np.zeros(grid.n)
        if raw[axis] is not None:
            comp = raw[axis]
            adv[ff] = 0.5 * (comp_self[ff] + comp[nk[ff], nj[ff], ni[ff]])
            adv[open_nb] = comp_self[open_nb]
        # Correction outward flux: (lam_n - lam_c) already encodes direction.
        corr = np.zeros(grid.n)
        corr[ff] = lam_field[nk[ff], nj[ff], ni[ff]] - lam_c[ff]
        corr[open_nb] = 0.0 - lam_c[open_nb]

        total += s * adv * area + coef * corr
    return float(np.max(np.abs(total)))


def wind_field_from_landscape(
    ls,
    *,
    speed: float,
    direction: float,
    z0: float | np.ndarray | None = None,
    **kwargs,
) -> WindField:
    """Mass-consistent wind for a :class:`~pyflam.landscape.Landscape`.

    Uses the landscape's elevation band and georeferencing. If ``z0`` is omitted
    it is derived from the fuel-model grid via :func:`roughness_from_fuel`.
    """
    if ls.elevation is None:
        raise ValueError("landscape has no elevation band for the wind solver")
    if z0 is None:
        z0 = roughness_from_fuel(ls.fuel_model)
    return solve_mass_consistent(
        np.asarray(ls.elevation, dtype=float), ls.cellsize_x,
        speed=speed, direction=direction, z0=z0,
        west=ls.west, north=ls.north, crs=ls.crs, **kwargs,
    )
