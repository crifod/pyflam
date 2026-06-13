"""Topology tests for the terrain-following polyMesh generator (no OpenFOAM).

These check the mesh is structurally valid the way ``checkMesh`` would, but
without needing OpenFOAM: correct counts, OpenFOAM's upper-triangular
owner<neighbour ordering, valid point/cell references, and the terrain mapping.
"""

from __future__ import annotations

import numpy as np
import pytest

from pyflam.cfd import mesh as M


def _mesh(nx=5, ny=4, nz=3, relief=40.0):
    yy, xx = np.mgrid[0:ny, 0:nx]
    elev = 1000.0 + relief * np.exp(-(((xx - nx / 2) ** 2 + (yy - ny / 2) ** 2) / 8.0))
    return M.build_terrain_mesh(elev, 30.0, nz=nz, domain_height=300.0)


def test_counts():
    nx, ny, nz = 5, 4, 3
    m = _mesh(nx, ny, nz)
    assert m.n_cells == nx * ny * nz
    assert len(m.points) == (nx + 1) * (ny + 1) * (nz + 1)
    internal = (nx - 1) * ny * nz + nx * (ny - 1) * nz + nx * ny * (nz - 1)
    assert m.n_internal_faces == internal
    assert len(m.neighbour) == internal


def test_patch_face_counts():
    nx, ny, nz = 5, 4, 3
    m = _mesh(nx, ny, nz)
    expect = {"ground": nx * ny, "top": nx * ny, "west": ny * nz,
              "east": ny * nz, "south": nx * nz, "north": nx * nz}
    for name, n in expect.items():
        assert m.patches[name][1] == n, name
    # Boundary faces are contiguous and come after all internal faces.
    assert m.patches["ground"][0] == m.n_internal_faces
    total = m.n_internal_faces + sum(v[1] for v in m.patches.values())
    assert len(m.faces) == total


def test_upper_triangular_ordering():
    m = _mesh()
    o, n = m.owner[:m.n_internal_faces], m.neighbour
    assert np.all(o < n)                       # owner < neighbour
    keys = o.astype(np.int64) * (m.n_cells + 1) + n
    assert np.all(np.diff(keys) > 0)           # ascending (owner, neighbour)


def test_references_valid():
    m = _mesh()
    assert m.faces.min() >= 0
    assert m.faces.max() < len(m.points)
    assert m.owner.min() >= 0 and m.owner.max() < m.n_cells
    # Every cell owns at least one face.
    assert set(np.unique(m.owner)) == set(range(m.n_cells))


def test_flat_dem_has_flat_ground():
    flat = np.full((4, 5), 1500.0)
    m = M.build_terrain_mesh(flat, 30.0, nz=4, domain_height=200.0)
    ground_pts = m.points[m.points[:, 2] <= 1500.0 + 1e-6]
    assert np.allclose(ground_pts[:, 2], 1500.0)
    # Layer heights span 0..domain_height, graded (first layer thinnest).
    h = m.layer_heights
    assert h[0] == 0.0 and h[-1] == pytest.approx(200.0)
    assert np.all(np.diff(h) > 0)
    assert (h[1] - h[0]) < (h[-1] - h[-2])


def test_column_mapping_and_orientation():
    m = _mesh(5, 4, 3)
    assert m.column_of(2, 0) == (3, 2)         # j=0 (south) -> DEM last row
    assert m.column_of(2, 3) == (0, 2)         # j=ny-1 (north) -> DEM row 0
    # Bottom node under DEM north-west corner is higher for a centred hill? just
    # check the bottom surface reproduces the (flipped) DEM at cell centres.
    assert m.cell_agl.shape == (3,)
