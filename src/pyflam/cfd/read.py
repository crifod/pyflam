"""Read an OpenFOAM result field back into a pyflam :class:`WindField`.

The case is built without ``renumberMesh``, so the solved ``U`` ``internalField``
is in our structured cell order (``cidx = k*ny*nx + j*nx + i``). We reshape it to
``(nz, ny, nx, 3)``, take each column's vertical profile, and interpolate the
horizontal wind to ``output_height`` above ground (log-law below the first cell
centre, matching :func:`pyflam.windsolver._extract_at_height`).
"""

from __future__ import annotations

import math
import os
import re

import numpy as np

from ..wind import WindField

_VEC_RE = re.compile(r"\(\s*([-\deE.+]+)\s+([-\deE.+]+)\s+([-\deE.+]+)\s*\)")


def latest_time_dir(case_dir: str) -> str:
    """Name of the highest numeric time directory (the converged result)."""
    times = []
    for entry in os.listdir(case_dir):
        if os.path.isdir(os.path.join(case_dir, entry)):
            try:
                times.append((float(entry), entry))
            except ValueError:
                pass
    if not times:
        raise FileNotFoundError(f"no time directories in {case_dir!r}")
    return max(times)[1]


def read_internal_vector(path: str, n_cells: int) -> np.ndarray:
    """Parse a volVectorField ``internalField`` (uniform or nonuniform), -> (N,3)."""
    with open(path) as fh:
        text = fh.read()
    m = re.search(r"internalField\s+(uniform|nonuniform)", text)
    if not m:
        raise ValueError(f"no internalField in {path!r}")
    if m.group(1) == "uniform":
        vec = _VEC_RE.search(text, m.end())
        v = [float(x) for x in vec.groups()]
        return np.tile(v, (n_cells, 1))
    # nonuniform: the data list follows "List<vector>".
    start = text.index("List<vector>", m.end())
    triples = _VEC_RE.findall(text, start)[:n_cells]
    arr = np.array(triples, dtype=float)
    if arr.shape != (n_cells, 3):
        raise ValueError(
            f"expected {n_cells} vectors in {path!r}, got {arr.shape[0]}")
    return arr


def read_wind_field(
    case_dir: str,
    mesh,
    *,
    output_height: float,
    z0,
    west: float = 0.0,
    north: float = 0.0,
    crs: object | None = None,
    time: str | None = None,
) -> WindField:
    """Build a :class:`WindField` from a solved case's ``U`` field.

    ``z0`` (m) is used for the near-ground log-law extrapolation when
    ``output_height`` is below the first cell centre; it may be a scalar or a 2D
    DEM-grid array (per-cell roughness).
    """
    time = time or latest_time_dir(case_dir)
    u = read_internal_vector(os.path.join(case_dir, time, "U"), mesh.n_cells)
    uf = u.reshape(mesh.nz, mesh.ny, mesh.nx, 3)
    z_agl = mesh.cell_agl                       # (nz,), ascending
    ue_all, vn_all = uf[..., 0], uf[..., 1]
    z0_arr = None if np.isscalar(z0) else np.asarray(z0, dtype=float)

    speed2d = np.zeros((mesh.ny, mesh.nx))
    dir2d = np.zeros((mesh.ny, mesh.nx))
    for j in range(mesh.ny):
        for i in range(mesh.nx):
            ue_col, vn_col = ue_all[:, j, i], vn_all[:, j, i]
            row, col = mesh.column_of(i, j)     # map mesh col -> DEM (row, col)
            if output_height < z_agl[0]:
                z0_local = float(z0) if z0_arr is None else float(z0_arr[row, col])
                factor = (math.log((output_height + z0_local) / z0_local)
                          / math.log((z_agl[0] + z0_local) / z0_local))
                ue, vn = ue_col[0] * factor, vn_col[0] * factor
            else:
                ue = np.interp(output_height, z_agl, ue_col)
                vn = np.interp(output_height, z_agl, vn_col)
            speed2d[row, col] = math.hypot(ue, vn)
            dir2d[row, col] = math.degrees(math.atan2(-ue, -vn)) % 360.0

    return WindField(
        speed=speed2d, direction=dir2d, cellsize=mesh.cellsize,
        west=west, north=north, speed_units="m/s",
        height=output_height, height_units="m", crs=crs,
    )
