"""Fire growth: directional elliptical spread + Minimum Travel Time (steps 3 & 5).

This is the spread/perimeter engine — the piece that turns per-cell spread
*rates* into a spreading *fire*: a fire-arrival-time surface across the
landscape, and the perimeters that fall out of it.

Two ideas, both from Finney:

* **Directional elliptical spread** (Finney 1998, FARSITE). Under a given wind
  and slope a point fire grows as an ellipse, fastest in the heading direction.
  pyflam combines the Rothermel wind and slope factors as *vectors* (wind
  blowing downwind, slope pushing upslope) to get the heading direction and the
  maximum spread rate, and derives the ellipse's eccentricity from the effective
  wind speed (Anderson 1983 length-to-breadth ratio). The spread rate in any
  direction ``psi`` off the heading is then
  ``R(psi) = R_max (1 - e) / (1 - e cos psi)``. This is also roadmap step 3
  (directional spread from wind bearing + aspect).

* **Minimum Travel Time** (Finney 2002). Fire arrival time is the shortest
  *time* path from the ignition to every cell, where the time to cross a
  straight segment is its length divided by the elliptical spread rate in that
  segment's direction. That is a shortest-path problem, solved here with
  Dijkstra over a lattice of travel directions. The result is identical in
  spirit to FlamMap's MTT output: an arrival-time raster from which fire size,
  perimeters and (later) flow paths are read.

Everything is in pyflam's native English units: spread rates ft/min, distances
ft, arrival times minutes, directions degrees clockwise from north.

References:
    Finney, M.A. 1998. FARSITE: Fire Area Simulator -- model development and
        evaluation. USDA Forest Service Research Paper RMRS-RP-4.
    Finney, M.A. 2002. Fire growth using minimum travel time methods. Canadian
        Journal of Forest Research 32: 1420-1424.
    Anderson, H.E. 1983. Predicting wind-driven wildland fire size and shape.
        USDA Forest Service Research Paper INT-305.
"""

from __future__ import annotations

import heapq
import math
from dataclasses import dataclass
from math import gcd

import numpy as np

from . import fuel_models
from .rothermel import kernel_param_groups, surface_kernel

try:                                    # optional: JIT the Eikonal sweep (pyflam[accel])
    import numba
    _HAVE_NUMBA = True
except Exception:                       # pragma: no cover - numba is optional
    _HAVE_NUMBA = False

# Anderson (1983) length-to-breadth ratio is capped here, as in FARSITE, so the
# ellipse can't become an unphysical sliver at extreme winds.
_MAX_LENGTH_BREADTH = 8.0


@dataclass
class SpreadField:
    """Per-cell elliptical spread template for a landscape.

    ``ros_max`` is the heading (maximum) spread rate (ft/min); ``heading`` is the
    azimuth it points *toward* (degrees clockwise from north); ``eccentricity``
    sets the ellipse shape (0 = circle, → 1 = long and thin). Nonburnable cells
    have ``ros_max == 0`` and act as barriers to fire growth.
    """

    ros_max: np.ndarray
    eccentricity: np.ndarray
    heading: np.ndarray          # degrees, direction of maximum spread (toward)
    cellsize_x: float
    cellsize_y: float
    fireline_intensity: np.ndarray | None = None  # Btu/ft/s at head, for spotting

    @property
    def shape(self) -> tuple[int, int]:
        return self.ros_max.shape

    def directional_ros(self, azimuth) -> np.ndarray:
        """Spread rate (ft/min) in compass direction ``azimuth`` per cell."""
        psi = np.radians(np.asarray(azimuth, dtype=float) - self.heading)
        e = self.eccentricity
        return self.ros_max * (1.0 - e) / (1.0 - e * np.cos(psi))


def _length_to_breadth(eff_wind_mph):
    """Anderson (1983) ellipse length-to-breadth ratio from effective wind."""
    u = np.asarray(eff_wind_mph, dtype=float)
    lb = 0.936 * np.exp(0.2566 * u) + 0.461 * np.exp(-0.1548 * u) - 0.397
    return np.clip(lb, 1.0, _MAX_LENGTH_BREADTH)


def spread_field(
    ls,
    *,
    m_1h: float,
    m_10h: float,
    m_100h: float,
    m_live_herb: float = 0.0,
    m_live_woody: float = 0.0,
    wind_midflame=0.0,
    wind_direction=0.0,
    load_factor: float = 1.0,
) -> SpreadField:
    """Build the per-cell elliptical spread template for a landscape.

    Combines the Rothermel wind factor (toward the downwind direction) and slope
    factor (toward upslope, from the landscape ``aspect``) as vectors to get each
    cell's maximum spread rate and heading. ``wind_midflame`` (ft/min) and
    ``wind_direction`` (degrees the wind blows *from*, met convention) may each
    be a scalar or a 2D field. If the landscape has no ``aspect`` band, the slope
    contribution is aligned with the wind (the maximum-spread assumption).

    The heavy Rothermel terms are computed once per unique fuel model and applied
    vectorized, exactly like :func:`pyflam.basic_fire_behavior`.
    """
    shape = ls.shape
    ros_max = np.zeros(shape, dtype=float)
    ecc = np.zeros(shape, dtype=float)
    heading = np.zeros(shape, dtype=float)
    fli = np.zeros(shape, dtype=float)        # head fireline intensity, Btu/ft/s

    tan_slope = ls.slope_tangent
    fuel = np.asarray(ls.fuel_model)
    wind_mid = np.broadcast_to(np.asarray(wind_midflame, dtype=float), shape)
    wind_dir = np.broadcast_to(np.asarray(wind_direction, dtype=float), shape)
    wind_toward = np.radians((wind_dir + 180.0) % 360.0)
    if ls.aspect is not None:
        upslope = np.radians((np.asarray(ls.aspect, dtype=float) + 180.0) % 360.0)
    else:
        upslope = wind_toward  # no aspect: push slope along the wind
    moist = dict(m_1h=m_1h, m_10h=m_10h, m_100h=m_100h,
                 m_live_herb=m_live_herb, m_live_woody=m_live_woody)

    for num in np.unique(fuel):
        num = int(num)
        mask = fuel == num
        try:
            fm = fuel_models.get(num)
        except KeyError:
            continue
        if not fm.is_burnable:
            continue

        lf = load_factor if not isinstance(load_factor, dict) \
            else load_factor.get(num, 1.0)
        for sub, p in kernel_param_groups(mask, {"load_factor": lf, **moist}):
            kernel = surface_kernel(fm, **p)
            phi_w = np.asarray(kernel.wind_factor(wind_mid[sub]), dtype=float)
            phi_s = np.asarray(kernel.slope_factor(tan_slope[sub]), dtype=float)

            # Vector sum of the wind push (downwind) and slope push (upslope).
            wt = wind_toward[sub]
            up = upslope[sub]
            vx = phi_w * np.sin(wt) + phi_s * np.sin(up)
            vy = phi_w * np.cos(wt) + phi_s * np.cos(up)
            phi_eff = np.hypot(vx, vy)
            head = np.degrees(np.arctan2(vx, vy)) % 360.0

            r = kernel.r0 * (1.0 + phi_eff)
            ros_max[sub] = r
            heading[sub] = head
            fli[sub] = kernel.heat_per_unit_area * r / 60.0   # Byram, Btu/ft/s

            # Effective wind speed = the wind that alone would give this combined
            # factor (invert the Rothermel wind factor), then Anderson L/B -> ecc.
            if kernel.c > 0.0 and kernel.b > 0.0:
                u_ftmin = np.where(
                    phi_eff > 0.0,
                    (phi_eff * kernel.beta_ratio ** kernel.e / kernel.c)
                    ** (1.0 / kernel.b),
                    0.0,
                )
            else:
                u_ftmin = np.zeros_like(phi_eff)
            lb = _length_to_breadth(u_ftmin / 88.0)        # ft/min -> mph
            ecc[sub] = np.sqrt(1.0 - 1.0 / lb ** 2)

    return SpreadField(
        ros_max=ros_max, eccentricity=ecc, heading=heading,
        cellsize_x=ls.cellsize_x, cellsize_y=ls.cellsize_y,
        fireline_intensity=fli,
    )


def _template(ring: int) -> list[tuple[int, int]]:
    """Primitive lattice travel directions out to ``ring`` cells.

    Only primitive (gcd==1) offsets are kept, so colinear duplicates like (2,0)
    are dropped while a spread of directions and segment lengths is retained.
    ``ring=2`` gives the 16-direction template; larger rings reduce the angular
    discretization error of the shortest-path search at higher cost.
    """
    offs = []
    for dr in range(-ring, ring + 1):
        for dc in range(-ring, ring + 1):
            if dr == 0 and dc == 0:
                continue
            if gcd(abs(dr), abs(dc)) == 1:
                offs.append((dr, dc))
    return offs


def _axis_slices(d: int, n: int) -> tuple[slice, slice]:
    """Aligned (source, target) slices along one axis for a shift of ``d``."""
    if d >= 0:
        return slice(0, n - d), slice(d, n)
    return slice(-d, n), slice(0, n + d)


def _offset_blocks(dr: int, dc: int, nrows: int, ncols: int,
                   r0: int, r1: int):
    """Aligned (source, target) index windows for offset ``(dr, dc)``.

    Restricts the source rows to the band ``[r0, r1)``. Returns
    ``(src_rows, src_cols, tgt_rows, tgt_cols)`` slices, or ``None`` if the band
    contributes no cells for this offset.
    """
    sr, tr = _axis_slices(dr, nrows)
    sc, tc = _axis_slices(dc, ncols)
    a = max(sr.start, r0)
    b = min(sr.stop, r1)
    if b <= a:
        return None
    off = tr.start - sr.start                      # target row = source row + off
    return (slice(a, b), sc, slice(a + off, b + off), tc)


def build_traveltime_graph(field: SpreadField, *, ring: int = 2,
                           chunk_rows: int | None = None):
    """Sparse directed graph of segment travel times (minutes) between cells.

    Node ``id = row * ncols + col``. For every lattice direction in the template,
    the directed edge ``A -> B`` is weighted by ``dist / R`` where ``R`` is the
    harmonic mean of the elliptical spread rate at ``A`` and ``B`` evaluated in
    the ``A -> B`` direction. Nonburnable cells contribute no edges, so they are
    barriers. The graph is built with vectorized NumPy slicing (no Python
    per-cell loop), which is what lets MTT scale to large grids.

    For huge / high-resolution landscapes the build is **chunked**: edges are
    counted first (using only the burnable mask -- an edge exists iff both
    endpoints burn), the exact COO arrays are preallocated once, then filled
    band-by-band. That avoids the ~2x peak of concatenating per-direction arrays,
    so memory stays close to the size of the final graph. ``chunk_rows`` sets the
    band height; ``None`` auto-chunks above a few million cells, ``0`` forces the
    single-shot path.

    Returns a :class:`scipy.sparse.csr_matrix` of shape ``(N, N)``.
    """
    from scipy import sparse

    nrows, ncols = field.shape
    n = nrows * ncols
    ros_max = np.asarray(field.ros_max, dtype=float)
    ecc = np.asarray(field.eccentricity, dtype=float)
    heading = np.asarray(field.heading, dtype=float)
    burnable = ros_max > 0.0
    # int32 node ids keep the edge arrays compact for huge grids; weights stay
    # float64 so arrival times are full precision.
    node_id = np.arange(n, dtype=np.int32).reshape(nrows, ncols)
    csx, csy = field.cellsize_x, field.cellsize_y
    offsets = _template(ring)

    # Per-offset distance and travel azimuth (deg from north).
    geom = [(dr, dc, math.hypot(dc * csx, -dr * csy),
             math.degrees(math.atan2(dc * csx, -dr * csy)) % 360.0)
            for dr, dc in offsets]

    if chunk_rows is None:
        # Auto: chunk above ~3M cells, ~1M cells per band.
        chunk_rows = max(1, 1_000_000 // ncols) if n > 3_000_000 else 0
    if chunk_rows <= 0 or chunk_rows >= nrows:
        bands = [(0, nrows)]
    else:
        bands = [(r0, min(r0 + chunk_rows, nrows))
                 for r0 in range(0, nrows, chunk_rows)]

    def edge_weights(dr, dc, dist, az, blocks):
        sr, sc, trr, tc = blocks
        e_s, h_s = ecc[sr, sc], heading[sr, sc]
        e_t, h_t = ecc[trr, tc], heading[trr, tc]
        with np.errstate(divide="ignore", invalid="ignore"):
            r_s = ros_max[sr, sc] * (1.0 - e_s) / (
                1.0 - e_s * np.cos(math.radians(az) - np.radians(h_s)))
            r_t = ros_max[trr, tc] * (1.0 - e_t) / (
                1.0 - e_t * np.cos(math.radians(az) - np.radians(h_t)))
            r_seg = 2.0 * r_s * r_t / (r_s + r_t)
            w = dist / r_seg
        return w

    # Pass 1: per-node out-degree, using only the burnable mask (an edge exists
    # iff both endpoints burn). This sizes the CSR exactly and lets us assemble
    # it directly -- no COO intermediate, so peak memory is ~the graph itself.
    deg = np.zeros(n, dtype=np.int64)
    for r0, r1 in bands:
        for dr, dc, _dist, _az in geom:
            blocks = _offset_blocks(dr, dc, nrows, ncols, r0, r1)
            if blocks is None:
                continue
            sr, sc, trr, tc = blocks
            valid = burnable[sr, sc] & burnable[trr, tc]
            if valid.any():
                ids = node_id[sr, sc][valid]
                deg += np.bincount(ids, minlength=n)

    indptr = np.empty(n + 1, dtype=np.int64)
    indptr[0] = 0
    np.cumsum(deg, out=indptr[1:])
    nnz = int(indptr[-1])
    if nnz == 0:
        return sparse.csr_matrix((n, n))

    indices = np.empty(nnz, dtype=np.int32)
    data = np.empty(nnz, dtype=np.float64)
    cursor = indptr[:-1].copy()            # running write position per source row

    # Pass 2: place each edge directly at its CSR slot. Within one offset every
    # source node appears at most once, so the scattered writes never collide.
    for r0, r1 in bands:
        for dr, dc, dist, az in geom:
            blocks = _offset_blocks(dr, dc, nrows, ncols, r0, r1)
            if blocks is None:
                continue
            sr, sc, trr, tc = blocks
            valid = burnable[sr, sc] & burnable[trr, tc]
            if not valid.any():
                continue
            ids = node_id[sr, sc][valid]
            w = edge_weights(dr, dc, dist, az, blocks)
            pos = cursor[ids]
            indices[pos] = node_id[trr, tc][valid]
            data[pos] = w[valid]
            cursor[ids] = pos + 1

    return sparse.csr_matrix((data, indices, indptr), shape=(n, n))


def _shift_inf(arr: np.ndarray, dr: int, dc: int) -> np.ndarray:
    """``out[r, c] = arr[r + dr, c + dc]`` with out-of-range entries set to inf."""
    nr, nc = arr.shape
    out = np.full_like(arr, np.inf)
    rd0, rd1 = max(0, -dr), min(nr, nr - dr)
    cd0, cd1 = max(0, -dc), min(nc, nc - dc)
    if rd0 < rd1 and cd0 < cd1:
        out[rd0:rd1, cd0:cd1] = arr[rd0 + dr:rd1 + dr, cd0 + dc:cd1 + dc]
    return out


# Eight compass neighbours ordered around the cell; consecutive pairs span the
# eight triangular simplices used by the semi-Lagrangian update.
_EIKONAL_NEIGHBORS = [(0, 1), (-1, 1), (-1, 0), (-1, -1),
                      (0, -1), (1, -1), (1, 0), (1, 1)]


def _eikonal_samples(csx, csy, alpha_samples):
    """Tabulate the (simplex, alpha) sources: neighbour offsets, distance, bearing."""
    dar, dac, dbr, dbc, als, dists, bears = [], [], [], [], [], [], []
    alphas = np.linspace(0.0, 1.0, alpha_samples)
    for i in range(8):
        a_off, b_off = _EIKONAL_NEIGHBORS[i], _EIKONAL_NEIGHBORS[(i + 1) % 8]
        ax, ay = a_off[1] * csx, -a_off[0] * csy
        bx, by = b_off[1] * csx, -b_off[0] * csy
        for al in alphas:
            yx, yy = al * ax + (1 - al) * bx, al * ay + (1 - al) * by
            dist = math.hypot(yx, yy)
            if dist == 0.0:
                continue
            bearing = math.degrees(math.atan2(-yx, -yy)) % 360.0
            dar.append(a_off[0]); dac.append(a_off[1])
            dbr.append(b_off[0]); dbc.append(b_off[1])
            als.append(float(al)); dists.append(dist); bears.append(bearing)
    return (np.array(dar, np.int64), np.array(dac, np.int64),
            np.array(dbr, np.int64), np.array(dbc, np.int64),
            np.array(als, np.float64), np.array(dists, np.float64),
            np.array(bears, np.float64))


_NB8 = np.array(_EIKONAL_NEIGHBORS, dtype=np.int64)   # 8 compass offsets, for the heap


if _HAVE_NUMBA:
    @numba.njit(cache=True)
    def _eikonal_heap_numba(arrival, ros_max, ecc, heading, dar, dac, dbr, dbc,
                            alphas, dists, bearings, nb8, src_r, src_c, max_time):
        """Heap-based narrow-band anisotropic Fast Marching of the semi-Lagrangian
        update. Each cell is *accepted* (finalised) the first time it is popped, and
        a cell's tentative value is built only from already-accepted neighbours --
        the Fast-Marching causality that makes acceptance final, so the non-decreasing
        pop order lets a finite ``max_time`` prune exactly like Dijkstra. Only the
        advancing front is ever touched, so it is far cheaper than sweeping the whole
        grid for a small bounded fire. (Accept-on-pop is exact for the isotropic part
        and an Ordered-Upwind approximation for the anisotropic part.)"""
        nr, nc = arrival.shape
        ns = alphas.shape[0]
        deg2rad = math.pi / 180.0
        accepted = np.zeros((nr, nc), np.uint8)
        cap = 1024
        ht = np.empty(cap, np.float64)
        hr = np.empty(cap, np.int64)
        hc = np.empty(cap, np.int64)
        size = 0
        for i in range(src_r.shape[0]):                 # push ignition cells (T=0)
            ht[size] = arrival[src_r[i], src_c[i]]
            hr[size] = src_r[i]; hc[size] = src_c[i]
            j = size; size += 1
            while j > 0:
                p = (j - 1) // 2
                if ht[p] <= ht[j]:
                    break
                ht[p], ht[j] = ht[j], ht[p]
                hr[p], hr[j] = hr[j], hr[p]
                hc[p], hc[j] = hc[j], hc[p]
                j = p

        while size > 0:
            t = ht[0]; r = hr[0]; c = hc[0]             # pop min
            size -= 1
            ht[0] = ht[size]; hr[0] = hr[size]; hc[0] = hc[size]
            i = 0
            while True:
                l = 2 * i + 1; rr = 2 * i + 2; sm = i
                if l < size and ht[l] < ht[sm]:
                    sm = l
                if rr < size and ht[rr] < ht[sm]:
                    sm = rr
                if sm == i:
                    break
                ht[sm], ht[i] = ht[i], ht[sm]
                hr[sm], hr[i] = hr[i], hr[sm]
                hc[sm], hc[i] = hc[i], hc[sm]
                i = sm
            if accepted[r, c]:
                continue                                # already finalised (stale)
            if t > max_time:
                break                                   # prune: nothing left is closer
            accepted[r, c] = 1                          # accept-on-pop: value is final

            for k in range(8):                          # relax the 8 neighbours of (r,c)
                vr = r + nb8[k, 0]; vc = c + nb8[k, 1]
                if vr < 0 or vr >= nr or vc < 0 or vc >= nc:
                    continue
                if accepted[vr, vc]:
                    continue
                rm = ros_max[vr, vc]
                if rm <= 0.0:
                    continue
                e = ecc[vr, vc]; hd = heading[vr, vc]
                best = arrival[vr, vc]
                improved = False
                for s in range(ns):
                    ra = vr + dar[s]; ca = vc + dac[s]
                    rb = vr + dbr[s]; cb = vc + dbc[s]
                    # only accepted neighbours carry usable (final) times
                    ta = (arrival[ra, ca] if 0 <= ra < nr and 0 <= ca < nc
                          and accepted[ra, ca] else 1e300)
                    tb = (arrival[rb, cb] if 0 <= rb < nr and 0 <= cb < nc
                          and accepted[rb, cb] else 1e300)
                    al = alphas[s]
                    if al == 0.0:
                        interp = tb
                    elif al == 1.0:
                        interp = ta
                    elif ta >= 1e300 or tb >= 1e300:
                        continue
                    else:
                        interp = al * ta + (1.0 - al) * tb
                    if interp >= 1e300:
                        continue
                    denom = 1.0 - e * math.cos(deg2rad * (bearings[s] - hd))
                    speed = rm * (1.0 - e) / denom
                    if speed <= 0.0:
                        continue
                    cand = interp + dists[s] / speed
                    if cand < best:
                        best = cand
                        improved = True
                if improved:
                    arrival[vr, vc] = best
                    if size >= cap:                     # grow heap arrays
                        ncap = cap * 2
                        nt = np.empty(ncap, np.float64); nt[:size] = ht[:size]; ht = nt
                        n2 = np.empty(ncap, np.int64); n2[:size] = hr[:size]; hr = n2
                        n3 = np.empty(ncap, np.int64); n3[:size] = hc[:size]; hc = n3
                        cap = ncap
                    ht[size] = best; hr[size] = vr; hc[size] = vc
                    j = size; size += 1
                    while j > 0:
                        p = (j - 1) // 2
                        if ht[p] <= ht[j]:
                            break
                        ht[p], ht[j] = ht[j], ht[p]
                        hr[p], hr[j] = hr[j], hr[p]
                        hc[p], hc[j] = hc[j], hc[p]
                        j = p
        return arrival

    @numba.njit(cache=True)
    def _eikonal_gauss_seidel(arrival, ros_max, ecc, heading, dar, dac, dbr, dbc,
                              alphas, dists, bearings, r0, r1, c0, c1,
                              max_passes, tol):
        """Gauss-Seidel fast sweeping of the semi-Lagrangian anisotropic-Eikonal
        update (in place), restricted to the ``[r0:r1, c0:c1]`` box. Four alternating
        sweep directions per pass propagate information across the box, so it
        converges in far fewer passes than a Jacobi iteration; the box lets the heap
        backend correct only the (small) burned region."""
        nr, nc = arrival.shape
        ns = alphas.shape[0]
        deg2rad = math.pi / 180.0
        for _ in range(max_passes):
            max_change = 0.0
            for sr in range(2):
                for sc in range(2):
                    rlo = r0 if sr == 0 else r1 - 1
                    rhi = r1 if sr == 0 else r0 - 1
                    rstep = 1 if sr == 0 else -1
                    clo = c0 if sc == 0 else c1 - 1
                    chi = c1 if sc == 0 else c0 - 1
                    cstep = 1 if sc == 0 else -1
                    for r in range(rlo, rhi, rstep):
                        for c in range(clo, chi, cstep):
                            rm = ros_max[r, c]
                            if rm <= 0.0:
                                continue
                            cur = arrival[r, c]
                            best = cur
                            e = ecc[r, c]
                            hd = heading[r, c]
                            for s in range(ns):
                                ra = r + dar[s]; ca = c + dac[s]
                                rb = r + dbr[s]; cb = c + dbc[s]
                                ta = (arrival[ra, ca]
                                      if 0 <= ra < nr and 0 <= ca < nc else 1e300)
                                tb = (arrival[rb, cb]
                                      if 0 <= rb < nr and 0 <= cb < nc else 1e300)
                                al = alphas[s]
                                if al == 0.0:
                                    interp = tb
                                elif al == 1.0:
                                    interp = ta
                                elif ta >= 1e300 or tb >= 1e300:
                                    interp = 1e300
                                else:
                                    interp = al * ta + (1.0 - al) * tb
                                if interp >= 1e300:
                                    continue
                                denom = 1.0 - e * math.cos(
                                    deg2rad * (bearings[s] - hd))
                                speed = rm * (1.0 - e) / denom
                                if speed <= 0.0:
                                    continue
                                cand = interp + dists[s] / speed
                                if cand < best:
                                    best = cand
                            if best < cur:
                                if cur - best > max_change:
                                    max_change = cur - best
                                arrival[r, c] = best
            if max_change < tol:
                break
        return arrival
else:                                   # pragma: no cover - exercised without numba
    _eikonal_gauss_seidel = None
    _eikonal_heap_numba = None


def anisotropic_eikonal(
    field: SpreadField,
    ignitions,
    *,
    max_time: float = math.inf,
    alpha_samples: int = 9,
    max_iter: int | None = None,
    tol: float = 1e-7,
    backend: str | None = None,
) -> np.ndarray:
    """Fire arrival time by an anisotropic-Eikonal (Finsler) front solver.

    An alternative to the Dijkstra-on-a-lattice :func:`minimum_travel_time` that
    solves the same problem as a front rather than a graph: the arrival time is the
    geodesic distance in the Finsler metric whose unit ball is the cell's elliptical
    fire shape. It uses a **semi-Lagrangian** update -- the characteristic may arrive
    from any point on the segment between two adjacent neighbours (a 1-D minimisation
    over ``alpha_samples`` interpolation points per simplex), not only from the
    discrete lattice directions -- iterated with fast sweeping. That removes most of
    MTT's angular (lattice) bias, so a calm fire stays a circle and a wind-driven
    fire a smooth ellipse instead of a faceted polygon (Sethian & Vladimirsky 2003;
    Mirebeau 2014; cf. the Randers-Finsler formulation of Gahtan et al. 2026).

    ``backend`` (all give the same arrival field): ``"heap"`` is a heap-based
    narrow-band Ordered-Upwind solve that only touches the advancing front, so a
    finite ``max_time`` **prunes like Dijkstra** -- much cheaper for a small bounded
    fire and the right choice for burn probability; ``"numba"`` is Gauss-Seidel
    sweeping (good when the whole grid burns); ``"numpy"`` is a portable vectorized
    sweep (no numba). ``None`` picks ``"heap"`` when numba is available, else
    ``"numpy"``. Speed is taken from the cell being updated, so on a *uniform* field
    this matches the analytic ``distance / R(bearing)`` solution closely; on
    heterogeneous fields it differs from MTT's harmonic-mean-of-endpoints by
    construction. ``minimum_travel_time(method="fast_marching")`` dispatches here.

    Returns the arrival-time array (inf for unreached cells / past ``max_time``).
    """
    nrows, ncols = field.shape
    ros_max = np.ascontiguousarray(field.ros_max, dtype=np.float64)
    ecc = np.ascontiguousarray(field.eccentricity, dtype=np.float64)
    heading = np.ascontiguousarray(field.heading, dtype=np.float64)
    burnable = ros_max > 0.0

    arrival = np.full(field.shape, np.inf, dtype=np.float64)
    seeded = False
    for r, c in ignitions:
        if 0 <= r < nrows and 0 <= c < ncols and burnable[r, c]:
            arrival[r, c] = 0.0
            seeded = True
    if not seeded:
        return arrival

    dar, dac, dbr, dbc, als, dists, bears = _eikonal_samples(
        field.cellsize_x, field.cellsize_y, alpha_samples)
    if max_iter is None:
        max_iter = 4 * (nrows + ncols)

    if backend is None:
        backend = "heap" if _HAVE_NUMBA else "numpy"
    if backend in ("heap", "numba") and not _HAVE_NUMBA:
        raise ImportError(
            f"backend={backend!r} needs numba: pip install 'pyflam[accel]'")

    if backend == "heap":
        src = [(r, c) for r, c in ignitions
               if 0 <= r < nrows and 0 <= c < ncols and burnable[r, c]]
        src_r = np.array([s[0] for s in src], dtype=np.int64)
        src_c = np.array([s[1] for s in src], dtype=np.int64)
        # Pass 1: heap Fast Marching identifies + prunes the burned region (cheap,
        # first-order). Pass 2: Gauss-Seidel corrects it to the accurate fixed point
        # over only that region's bounding box -- sweep accuracy at pruned cost.
        _eikonal_heap_numba(arrival, ros_max, ecc, heading, dar, dac, dbr, dbc,
                            als, dists, bears, _NB8, src_r, src_c, float(max_time))
        rows, cols = np.where(np.isfinite(arrival))
        if rows.size:
            r0 = max(0, int(rows.min()) - 2); r1 = min(nrows, int(rows.max()) + 3)
            c0 = max(0, int(cols.min()) - 2); c1 = min(ncols, int(cols.max()) + 3)
            _eikonal_gauss_seidel(arrival, ros_max, ecc, heading, dar, dac, dbr, dbc,
                                  als, dists, bears, r0, r1, c0, c1, max_iter, tol)
    elif backend == "numba":
        _eikonal_gauss_seidel(arrival, ros_max, ecc, heading, dar, dac, dbr, dbc,
                              als, dists, bears, 0, nrows, 0, ncols, max_iter, tol)
    elif backend == "numpy":
        _eikonal_numpy(arrival, ros_max, ecc, heading, dar, dac, dbr, dbc,
                       als, dists, bears, max_iter, tol)
    else:
        raise ValueError("backend must be None, 'heap', 'numba' or 'numpy'")

    if math.isfinite(max_time):
        arrival[arrival > max_time] = np.inf
    return arrival


def _eikonal_numpy(arrival, ros_max, ecc, heading, dar, dac, dbr, dbc,
                   als, dists, bears, max_iter, tol):
    """Portable vectorized fallback for :func:`anisotropic_eikonal` (no numba)."""
    ns = als.shape[0]
    prev_finite = -1
    for it in range(max_iter):
        before = arrival.copy()
        idx = range(ns) if it % 2 == 0 else range(ns - 1, -1, -1)
        for s in idx:
            psi = np.radians(bears[s] - heading)
            speed = ros_max * (1.0 - ecc) / (1.0 - ecc * np.cos(psi))
            with np.errstate(divide="ignore", invalid="ignore"):
                step = np.where(speed > 0.0, dists[s] / speed, np.inf)
            ta = _shift_inf(arrival, dar[s], dac[s])
            tb = _shift_inf(arrival, dbr[s], dbc[s])
            al = als[s]
            if al == 0.0:
                interp = tb
            elif al == 1.0:
                interp = ta
            else:
                interp = al * ta + (1.0 - al) * tb
            np.minimum(arrival, interp + step, out=arrival)
        nfin = int(np.isfinite(arrival).sum())
        both = np.isfinite(arrival) & np.isfinite(before)
        delta = float(np.max(before[both] - arrival[both])) if both.any() else 0.0
        if nfin == prev_finite and delta < tol:
            break
        prev_finite = nfin
    return arrival


def minimum_travel_time(
    field: SpreadField,
    ignitions,
    *,
    max_time: float = math.inf,
    ring: int = 2,
    chunk_rows: int | None = None,
    method: str = "mtt",
) -> np.ndarray:
    """Fire arrival time (minutes) by Minimum Travel Time (Finney 2002).

    ``ignitions`` is an iterable of ``(row, col)`` cell indices that start at
    time 0. The time to cross a straight segment is its length divided by the
    harmonic mean of the elliptical spread rate at its two endpoints (the right
    average for travel time over the segment); the minimum-time path to every
    cell is then a shortest-path problem. Nonburnable cells (``ros_max == 0``)
    are barriers.

    ``method`` selects the propagation engine: ``"mtt"`` (default) is the
    Dijkstra-on-a-lattice solver below; ``"fast_marching"`` dispatches to
    :func:`anisotropic_eikonal`, a semi-Lagrangian anisotropic-Eikonal front solver
    that removes most of MTT's lattice (angular) bias at the cost of being a pure
    NumPy prototype rather than the C-level Dijkstra. ``ring``/``chunk_rows`` apply
    to the ``"mtt"`` path only.

    The graph is built vectorized (chunked for huge grids, see
    :func:`build_traveltime_graph` and ``chunk_rows``) and the shortest path is
    solved by SciPy's C-level multi-source Dijkstra, so this scales to large /
    high-resolution landscapes (many millions of cells). ``max_time`` bounds the
    search (cells that would arrive later are left ``inf``), which also makes
    bounded-duration runs much cheaper.

    Returns a float array (same shape as the field); unburned/unreachable cells
    and anything past ``max_time`` are ``inf``.
    """
    if method == "fast_marching":
        return anisotropic_eikonal(field, ignitions, max_time=max_time)
    if method != "mtt":
        raise ValueError("method must be 'mtt' or 'fast_marching'")

    from scipy.sparse.csgraph import dijkstra

    nrows, ncols = field.shape
    n = nrows * ncols
    burnable = np.asarray(field.ros_max) > 0.0

    sources = [
        r * ncols + c
        for r, c in ignitions
        if 0 <= r < nrows and 0 <= c < ncols and burnable[r, c]
    ]
    if not sources:
        return np.full(field.shape, math.inf, dtype=float)

    graph = build_traveltime_graph(field, ring=ring, chunk_rows=chunk_rows)
    dist = dijkstra(
        graph, directed=True, indices=np.asarray(sources, dtype=np.int64),
        min_only=True, limit=max_time,
    )
    return dist.reshape(nrows, ncols)


def _dijkstra_with_start_times(graph, n, sources, start_times, limit):
    """Multi-source Dijkstra where each source activates at its own start time.

    Uses a virtual super-source connected to each real source by an edge equal to
    that source's start time, so a single shortest-path solve yields
    ``min_s (start_time_s + travel_s->cell)`` -- exactly fire arrival time when
    several fronts (a primary ignition and spot fires) start at different times.
    """
    from scipy import sparse
    from scipy.sparse.csgraph import dijkstra

    sources = np.asarray(sources, dtype=np.int64)
    start_times = np.asarray(start_times, dtype=float)
    super_row = sparse.csr_matrix(
        (start_times, (np.zeros(sources.size, dtype=np.int64), sources)),
        shape=(1, n))
    top = sparse.hstack([graph, sparse.csr_matrix((n, 1))], format="csr")
    bottom = sparse.hstack([super_row, sparse.csr_matrix((1, 1))], format="csr")
    aug = sparse.vstack([top, bottom], format="csr")
    dist = dijkstra(aug, directed=True, indices=n, min_only=True, limit=limit)
    return dist[:n]


def spread_with_spotting(
    field: SpreadField,
    ignitions,
    *,
    max_time: float,
    wind_20ft: float,
    wind_direction: float,
    model,
    rng=None,
    ring: int = 2,
    chunk_rows: int | None = None,
    max_generations: int = 5,
    graph=None,
    fuel_moisture=None,
):
    """MTT fire growth with ember spotting (arrival time, minutes).

    Grows the primary fire, then repeatedly lofts embers from the burning area
    (``model.generate_spots``), adds the surviving spot landings as new ignitions
    at their landing times, and re-grows -- so the fire can jump fuel barriers.
    Stops when a generation adds no new burned cells or after ``max_generations``.

    ``model`` is a :class:`pyflam.spotting.SpottingModel` (parameterized) or
    :class:`pyflam.spotting.FirebrandPhysics` (stochastic, physics-based);
    ``wind_20ft`` (ft/min) and ``wind_direction`` (deg FROM) drive firebrand
    transport, and ``fuel_moisture`` (dead 1-h fraction; scalar or array) governs
    landing ignition in the physics model. Returns the arrival-time array.
    """
    if rng is None:
        rng = np.random.default_rng()
    nrows, ncols = field.shape
    n = nrows * ncols
    burnable = np.asarray(field.ros_max) > 0.0
    if graph is None:
        graph = build_traveltime_graph(field, ring=ring, chunk_rows=chunk_rows)

    sources, starts = [], []
    for r, c in ignitions:
        if 0 <= r < nrows and 0 <= c < ncols and burnable[r, c]:
            sources.append(r * ncols + c)
            starts.append(0.0)
    if not sources:
        return np.full(field.shape, math.inf, dtype=float)

    arrival = _dijkstra_with_start_times(
        graph, n, sources, starts, max_time).reshape(field.shape)
    burned = np.isfinite(arrival) & (arrival <= max_time)

    for _ in range(max_generations):
        spots = model.generate_spots(
            field, arrival, wind_20ft=wind_20ft, wind_direction=wind_direction,
            max_time=max_time, rng=rng, fuel_moisture=fuel_moisture,
        )
        if not spots:
            break
        for r, c, t in spots:
            sources.append(r * ncols + c)
            starts.append(t)
        arrival = _dijkstra_with_start_times(
            graph, n, sources, starts, max_time).reshape(field.shape)
        new_burned = np.isfinite(arrival) & (arrival <= max_time)
        if int(new_burned.sum()) == int(burned.sum()):
            break
        burned = new_burned

    return arrival


# FlamMap-style flame-length-probability (FLP) classes: 20 bins of 2 ft from 0 to
# 40 ft with an open-ended top class, matching the ``FIL1..FIL20`` columns of a
# FlamMap MTT random-ignition ``FLP_METRIC`` table. Override via
# ``flame_length_classes`` to match a run configured with different categories.
DEFAULT_FLAME_LENGTH_CLASSES = np.array(
    [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0,
     22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, math.inf])


@dataclass
class BurnProbabilityResult:
    """Burn probability plus the connected per-cell fire-behavior metrics.

    Returned by :func:`burn_probability` when ``return_metrics=True`` -- the
    pyflam analog of the full output set of a FlamMap MTT Random-Ignition run
    (``BURN_PROB``, ``FIRE_LINE_INT`` and ``FLP_METRIC``). All "conditional"
    rasters are *given the cell burned* (i.e. averaged only over the fires that
    reached the cell) and are ``nan`` where the cell never burned.
    """

    burn_prob: np.ndarray              # [0, 1] fraction of fires that reached each cell
    n_fires: int                       # fires actually ignited
    conditional_flame_length: np.ndarray   # ft, mean flame length given the cell burned
    conditional_intensity: np.ndarray      # Btu/ft/s, mean fireline intensity given burned
    flame_length_classes: np.ndarray       # class edges (ft), length n_classes + 1
    flp: np.ndarray                    # (n_classes, nrows, ncols) P(class | burned)
    fire_sizes: np.ndarray             # per-fire burned area (landscape units squared)

    @property
    def shape(self) -> tuple[int, int]:
        return self.burn_prob.shape

    def flame_length_class_centers(self) -> np.ndarray:
        """Representative flame length (ft) for each FLP class (top class = its low edge)."""
        edges = self.flame_length_classes
        centers = 0.5 * (edges[:-1] + edges[1:])
        if not math.isfinite(centers[-1]):
            centers[-1] = edges[-2]
        return centers


def _normalize_scenarios(field, graph):
    """Coerce ``field`` into a list of weather-scenario dicts.

    ``field`` is either one :class:`SpreadField` (a single deterministic weather)
    or a sequence describing a weather ensemble, whose entries are each a
    ``SpreadField``, a ``(weight, field)`` pair, a ``(weight, field, graph)``
    triple, or a ``dict`` with a ``"field"`` key and optional ``"weight"``,
    ``"graph"``, ``"wind_20ft"`` and ``"wind_direction"`` keys. The dict form is
    what lets a scenario carry its **own spotting wind** (ft/min and deg FROM); a
    scenario that omits them falls back to ``burn_probability``'s ``wind_20ft`` /
    ``wind_direction``. Each scenario keeps its own (lazily built) travel-time
    graph.
    """
    def make(weight, f, g, w20=None, wdir=None):
        return {"weight": float(weight), "field": f, "graph": g,
                "wind_20ft": w20, "wind_direction": wdir}

    if isinstance(field, SpreadField):
        return [make(1.0, field, graph)]
    scenarios = []
    for entry in field:
        if isinstance(entry, SpreadField):
            scenarios.append(make(1.0, entry, None))
        elif isinstance(entry, dict):
            if "field" not in entry:
                raise ValueError("scenario dict must have a 'field' key")
            scenarios.append(make(
                entry.get("weight", 1.0), entry["field"], entry.get("graph"),
                entry.get("wind_20ft"), entry.get("wind_direction")))
        else:
            entry = tuple(entry)
            if len(entry) == 2:
                scenarios.append(make(entry[0], entry[1], None))
            elif len(entry) == 3:
                scenarios.append(make(entry[0], entry[1], entry[2]))
            else:
                raise ValueError(
                    "scenario entries must be a SpreadField, a dict, "
                    "(weight, field) or (weight, field, graph)")
    if not scenarios:
        raise ValueError("no weather scenarios provided")
    return scenarios


def burn_probability(
    field,
    ignitions,
    *,
    max_time: float,
    ring: int = 2,
    chunk_rows: int | None = None,
    graph=None,
    spotting=None,
    wind_20ft: float = 0.0,
    wind_direction: float = 0.0,
    fuel_moisture=None,
    rng=None,
    return_metrics: bool = False,
    flame_length_classes=None,
    batch_size: int | None = None,
):
    """Burn probability + connected metrics from many MTT fires (FlamMap MTT BP).

    Each entry in ``ignitions`` is a *separate* fire (a ``(row, col)`` ignition);
    every fire grows for ``max_time`` minutes and the per-cell burn probability is
    the fraction of fires that reached the cell -- pyflam's analog of FlamMap's
    random-ignition burn-probability run.

    **Weather variation (accuracy).** ``field`` is either one :class:`SpreadField`
    (one deterministic weather, as before) or a *weather ensemble*: a sequence of
    ``SpreadField`` / ``(weight, field)`` / ``(weight, field, graph)`` scenarios,
    or scenario ``dict``\\s (``{"field":..., "weight":..., "graph":...,
    "wind_20ft":..., "wind_direction":...}``). Each fire is assigned a scenario
    drawn in proportion to its weight, so the fire-to-fire weather variation that
    makes FlamMap's flame-length distribution (FLP) meaningful is reproduced. One
    travel-time graph is built (and reused) per scenario, so an ensemble of a few
    scenarios stays cheap.

    **Spotting.** Pass a :class:`pyflam.spotting.SpottingModel` /
    :class:`pyflam.spotting.FirebrandPhysics` as ``spotting`` (with ``wind_20ft``
    in ft/min and ``wind_direction`` in deg FROM) to grow each fire *with ember
    spotting* via :func:`spread_with_spotting`, letting fires cross fuel barriers.
    In an ensemble, a scenario's ``dict`` may carry its own ``wind_20ft`` /
    ``wind_direction`` so ember transport uses that scenario's wind (otherwise the
    call-level values apply to every scenario).

    **Speed.** Without spotting the fires are solved in batches with one
    multi-source SciPy Dijkstra call per batch (set ``batch_size``; the default is
    chosen to bound peak memory), which is markedly faster than one solve per fire.

    Returns ``(burn_prob, n_fires)`` by default: a float array (cells in [0, 1])
    and the number of fires actually ignited (ignitions on nonburnable cells are
    skipped). With ``return_metrics=True`` returns a :class:`BurnProbabilityResult`
    that also carries conditional flame length, conditional fireline intensity and
    the per-class flame-length probabilities (``flame_length_classes`` overrides
    the default 20 FlamMap-style bins).
    """
    from scipy.sparse.csgraph import dijkstra

    if rng is None:
        rng = np.random.default_rng()

    scenarios = _normalize_scenarios(field, graph)
    nrows, ncols = scenarios[0]["field"].shape
    n = nrows * ncols
    base = scenarios[0]["field"]
    cellarea = base.cellsize_x * base.cellsize_y

    if flame_length_classes is None:
        edges = DEFAULT_FLAME_LENGTH_CLASSES
    else:
        edges = np.asarray(flame_length_classes, dtype=float)
    n_classes = edges.size - 1

    # Auto batch size keeps the dense (batch x n_cells) Dijkstra output near a
    # ~256 MB ceiling; tiny grids fall back to a modest cap.
    if batch_size is None:
        batch_size = max(1, min(256, int(256e6 / (8.0 * max(n, 1)))))

    # Per-scenario flattened rasters used to score burned cells.
    for sc in scenarios:
        f = sc["field"]
        ros = np.asarray(f.ros_max, dtype=float).reshape(-1)
        sc["burnable"] = ros > 0.0
        if return_metrics:
            fli = (np.zeros(n) if f.fireline_intensity is None
                   else np.asarray(f.fireline_intensity, dtype=float).reshape(-1))
            fl = np.where(fli > 0.0, 0.45 * fli ** 0.46, 0.0)
            sc["fli"] = fli
            sc["fl"] = fl
            sc["cls"] = np.clip(np.digitize(fl, edges) - 1, 0, n_classes - 1)

    def scenario_graph(sc):
        if sc["graph"] is None:
            sc["graph"] = build_traveltime_graph(
                sc["field"], ring=ring, chunk_rows=chunk_rows)
        return sc["graph"]

    # Assign each ignition to a weather scenario (drawn by weight).
    ignitions = list(ignitions)
    weights = np.array([sc["weight"] for sc in scenarios], dtype=float)
    weights /= weights.sum()
    if len(scenarios) == 1:
        scen_idx = np.zeros(len(ignitions), dtype=int)
    else:
        scen_idx = rng.choice(len(scenarios), size=len(ignitions), p=weights)

    count = np.zeros(n, dtype=np.int64)
    fli_sum = np.zeros(n) if return_metrics else None
    fl_sum = np.zeros(n) if return_metrics else None
    flp_count = np.zeros((n_classes, n), dtype=np.int64) if return_metrics else None
    fire_sizes: list[float] = []
    cell_ids = np.arange(n) if return_metrics else None

    def accumulate(sc, burned):
        """Fold a ``(k, n)`` boolean burned matrix (k fires) into the running totals."""
        per_cell = burned.sum(axis=0).astype(np.int64)
        count[:] += per_cell
        fire_sizes.extend((burned.sum(axis=1) * cellarea).tolist())
        if return_metrics:
            fli_sum[:] += per_cell * sc["fli"]
            fl_sum[:] += per_cell * sc["fl"]
            np.add.at(flp_count, (sc["cls"], cell_ids), per_cell)

    n_fires = 0
    for si, sc in enumerate(scenarios):
        burnable = sc["burnable"]
        sources = []
        for i, (r, c) in enumerate(ignitions):
            if scen_idx[i] != si:
                continue
            if 0 <= r < nrows and 0 <= c < ncols and burnable[r * ncols + c]:
                sources.append(r * ncols + c)
        if not sources:
            continue
        g = scenario_graph(sc)

        if spotting is not None:
            # Per-scenario spotting wind if the scenario set it, else the call's.
            sc_w20 = wind_20ft if sc["wind_20ft"] is None else sc["wind_20ft"]
            sc_wdir = (wind_direction if sc["wind_direction"] is None
                       else sc["wind_direction"])
            for src in sources:
                r, c = divmod(src, ncols)
                arrival = spread_with_spotting(
                    sc["field"], [(r, c)], max_time=max_time, wind_20ft=sc_w20,
                    wind_direction=sc_wdir, model=spotting, rng=rng,
                    graph=g, fuel_moisture=fuel_moisture)
                burned = (np.isfinite(arrival) & (arrival <= max_time)).reshape(1, n)
                accumulate(sc, burned)
                n_fires += 1
        else:
            srcarr = np.asarray(sources, dtype=np.int64)
            for start in range(0, srcarr.size, batch_size):
                chunk = srcarr[start:start + batch_size]
                dist = dijkstra(g, directed=True, indices=chunk,
                                min_only=False, limit=max_time)
                accumulate(sc, dist <= max_time)
                n_fires += chunk.size

    prob = (count / n_fires if n_fires else count.astype(float))
    prob = prob.reshape(nrows, ncols)
    if not return_metrics:
        return prob, n_fires

    with np.errstate(invalid="ignore", divide="ignore"):
        safe = np.where(count > 0, count, 1)
        cfl = np.where(count > 0, fl_sum / safe, np.nan).reshape(nrows, ncols)
        cint = np.where(count > 0, fli_sum / safe, np.nan).reshape(nrows, ncols)
        flp = np.where(count > 0, flp_count / safe, 0.0).reshape(n_classes, nrows, ncols)
    return BurnProbabilityResult(
        burn_prob=prob,
        n_fires=n_fires,
        conditional_flame_length=cfl,
        conditional_intensity=cint,
        flame_length_classes=edges,
        flp=flp,
        fire_sizes=np.asarray(fire_sizes, dtype=float),
    )


def _mtt_python(
    field: SpreadField,
    ignitions,
    *,
    max_time: float = math.inf,
    ring: int = 2,
) -> np.ndarray:
    """Pure-Python heap Dijkstra reference for :func:`minimum_travel_time`.

    Kept for cross-checking the vectorized/SciPy path in the tests; the public
    function above is the one to use (this one does not scale).
    """
    nrows, ncols = field.shape
    arrival = np.full(field.shape, math.inf, dtype=float)
    burnable = field.ros_max > 0.0

    ros_max = field.ros_max
    ecc = field.eccentricity
    heading = field.heading
    csx, csy = field.cellsize_x, field.cellsize_y
    geom = []
    for dr, dc in _template(ring):
        dx, dy = dc * csx, -dr * csy
        geom.append((dr, dc, math.hypot(dx, dy),
                     math.degrees(math.atan2(dx, dy)) % 360.0))

    heap: list[tuple[float, int, int]] = []
    for r, c in ignitions:
        if 0 <= r < nrows and 0 <= c < ncols and burnable[r, c]:
            arrival[r, c] = 0.0
            heapq.heappush(heap, (0.0, r, c))

    while heap:
        t, r, c = heapq.heappop(heap)
        if t > arrival[r, c] or t > max_time:
            continue
        ecc_s, head_s, rmax_s = ecc[r, c], heading[r, c], ros_max[r, c]
        for dr, dc, dist, az in geom:
            nr, nc = r + dr, c + dc
            if nr < 0 or nr >= nrows or nc < 0 or nc >= ncols:
                continue
            if not burnable[nr, nc]:
                continue
            psi_s = math.radians(az - head_s)
            r_src = rmax_s * (1.0 - ecc_s) / (1.0 - ecc_s * math.cos(psi_s))
            psi_t = math.radians(az - heading[nr, nc])
            e_t = ecc[nr, nc]
            r_tgt = ros_max[nr, nc] * (1.0 - e_t) / (1.0 - e_t * math.cos(psi_t))
            if r_src <= 0.0 or r_tgt <= 0.0:
                continue
            r_seg = 2.0 / (1.0 / r_src + 1.0 / r_tgt)
            nt = t + dist / r_seg
            if nt < arrival[nr, nc] and nt <= max_time:
                arrival[nr, nc] = nt
                heapq.heappush(heap, (nt, nr, nc))

    return arrival


def perimeter_mask(arrival_time: np.ndarray, time: float) -> np.ndarray:
    """Boolean burned area at ``time`` minutes (cells reached by then)."""
    return np.asarray(arrival_time) <= time


def burned_area(arrival_time: np.ndarray, time: float,
                cellsize_x: float, cellsize_y: float) -> float:
    """Burned area (in the landscape's linear units squared) at ``time``."""
    cells = int(np.count_nonzero(perimeter_mask(arrival_time, time)))
    return cells * cellsize_x * cellsize_y


def ignition_from_xy(ls, x: float, y: float) -> tuple[int, int]:
    """Convert world (x, y) to a ``(row, col)`` ignition index on ``ls``."""
    col = int((x - ls.west) / ls.cellsize_x)
    row = int((ls.north - y) / ls.cellsize_y)
    return row, col


def spread_perimeter(
    ls,
    ignitions,
    *,
    m_1h: float,
    m_10h: float,
    m_100h: float,
    m_live_herb: float = 0.0,
    m_live_woody: float = 0.0,
    wind_midflame=0.0,
    wind_direction=0.0,
    max_time: float = math.inf,
    ring: int = 2,
    chunk_rows: int | None = None,
    method: str = "mtt",
) -> dict[str, np.ndarray]:
    """One-call fire growth: build the spread field and run MTT over a landscape.

    ``ignitions`` is an iterable of ``(row, col)`` indices (use
    :func:`ignition_from_xy` for world coordinates). Returns a dict with the
    ``arrival_time`` raster (minutes), the ``spread_field`` (a
    :class:`SpreadField`) and the heading ``ros_max`` raster (ft/min).

    ``method`` chooses the propagation engine (``"mtt"`` Dijkstra, default, or
    ``"fast_marching"`` anisotropic-Eikonal; see :func:`minimum_travel_time`).
    ``chunk_rows`` controls the chunked graph build for huge landscapes (see
    :func:`build_traveltime_graph`); leave ``None`` to auto-chunk.
    """
    field = spread_field(
        ls, m_1h=m_1h, m_10h=m_10h, m_100h=m_100h,
        m_live_herb=m_live_herb, m_live_woody=m_live_woody,
        wind_midflame=wind_midflame, wind_direction=wind_direction,
    )
    arrival = minimum_travel_time(
        field, ignitions, max_time=max_time, ring=ring, chunk_rows=chunk_rows,
        method=method,
    )
    return {
        "arrival_time": arrival,
        "spread_field": field,
        "ros_max": field.ros_max,
    }
