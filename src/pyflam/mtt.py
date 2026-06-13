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


def minimum_travel_time(
    field: SpreadField,
    ignitions,
    *,
    max_time: float = math.inf,
    ring: int = 2,
    chunk_rows: int | None = None,
) -> np.ndarray:
    """Fire arrival time (minutes) by Minimum Travel Time (Finney 2002).

    ``ignitions`` is an iterable of ``(row, col)`` cell indices that start at
    time 0. The time to cross a straight segment is its length divided by the
    harmonic mean of the elliptical spread rate at its two endpoints (the right
    average for travel time over the segment); the minimum-time path to every
    cell is then a shortest-path problem. Nonburnable cells (``ros_max == 0``)
    are barriers.

    The graph is built vectorized (chunked for huge grids, see
    :func:`build_traveltime_graph` and ``chunk_rows``) and the shortest path is
    solved by SciPy's C-level multi-source Dijkstra, so this scales to large /
    high-resolution landscapes (many millions of cells). ``max_time`` bounds the
    search (cells that would arrive later are left ``inf``), which also makes
    bounded-duration runs much cheaper.

    Returns a float array (same shape as the field); unburned/unreachable cells
    and anything past ``max_time`` are ``inf``.
    """
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


def burn_probability(
    field: SpreadField,
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
):
    """Burn probability from many fixed-duration MTT fires (FlamMap-style MTT BP).

    Each entry in ``ignitions`` is a *separate* fire (a ``(row, col)`` ignition);
    every fire grows for ``max_time`` minutes and the per-cell burn probability is
    the fraction of fires that reached the cell. This is the deterministic-weather
    analog of FlamMap's random-ignition burn-probability run (fire-to-fire weather
    variation is not modelled).

    Pass a :class:`pyflam.spotting.SpottingModel` as ``spotting`` (with
    ``wind_20ft`` in ft/min and ``wind_direction`` in deg FROM) to grow each fire
    *with ember spotting* via :func:`spread_with_spotting` -- which is what lets
    fires cross fuel barriers and reproduces a spotting-on FlamMap run. Without it
    each fire is a plain contiguous MTT fire.

    The travel-time graph is built **once** and reused for every fire, so
    thousands of bounded fires are cheap. Pass a prebuilt ``graph`` to reuse it.

    Returns ``(burn_prob, n_fires)``: a float array (cells in [0, 1]) and the
    number of fires actually ignited (ignitions on nonburnable cells are skipped).
    """
    from scipy.sparse.csgraph import dijkstra

    nrows, ncols = field.shape
    n = nrows * ncols
    burnable = np.asarray(field.ros_max) > 0.0
    if graph is None:
        graph = build_traveltime_graph(field, ring=ring, chunk_rows=chunk_rows)
    if rng is None:
        rng = np.random.default_rng()

    count = np.zeros((nrows, ncols), dtype=np.int32)
    n_fires = 0
    for r, c in ignitions:
        if not (0 <= r < nrows and 0 <= c < ncols and burnable[r, c]):
            continue
        if spotting is not None:
            arrival = spread_with_spotting(
                field, [(r, c)], max_time=max_time, wind_20ft=wind_20ft,
                wind_direction=wind_direction, model=spotting, rng=rng,
                graph=graph, fuel_moisture=fuel_moisture)
            count += np.isfinite(arrival) & (arrival <= max_time)
        else:
            dist = dijkstra(graph, directed=True, indices=r * ncols + c,
                            min_only=True, limit=max_time)
            count += (dist <= max_time).reshape(nrows, ncols)
        n_fires += 1

    prob = (count / n_fires if n_fires else count).astype(float)
    return prob, n_fires


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
) -> dict[str, np.ndarray]:
    """One-call fire growth: build the spread field and run MTT over a landscape.

    ``ignitions`` is an iterable of ``(row, col)`` indices (use
    :func:`ignition_from_xy` for world coordinates). Returns a dict with the
    ``arrival_time`` raster (minutes), the ``spread_field`` (a
    :class:`SpreadField`) and the heading ``ros_max`` raster (ft/min).

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
    )
    return {
        "arrival_time": arrival,
        "spread_field": field,
        "ros_max": field.ros_max,
    }
