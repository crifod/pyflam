"""Terrain-following structured hexahedral mesh -> OpenFOAM ``polyMesh``.

Builds a graded hex mesh whose bottom follows the DEM and whose layers are
geometrically refined near the ground (for the rough-wall functions), then writes
it as an OpenFOAM ``constant/polyMesh`` (``points``, ``faces``, ``owner``,
``neighbour``, ``boundary``).

Conventions used throughout:
    * local coordinates: x = east (DEM column ``i``), y = north (``j``), z = up.
    * the DEM ``elevation[row, col]`` has ``row 0 = north``; we flip rows so mesh
      ``j`` increases northward. Cell ``(i, j, k)`` maps back to DEM
      ``row = ny-1-j``, ``col = i`` (see :meth:`TerrainMesh.column_of`).
    * six boundary patches: ``ground`` (a ``wall``) and ``top``/``west``/``east``/
      ``south``/``north`` (generic ``patch`` — field files use ``slip`` on the
      non-inlet ones, so no geometric ``symmetry`` patch type is required).

Internal faces are emitted in OpenFOAM's required upper-triangular order
(ascending owner, then neighbour). Validated with ``checkMesh``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import dictwriter

GROUND, TOP, WEST, EAST, SOUTH, NORTH = (
    "ground", "top", "west", "east", "south", "north")
# Patch write order (boundary file + face blocks).
PATCHES = (GROUND, TOP, WEST, EAST, SOUTH, NORTH)
LATERAL = (WEST, EAST, SOUTH, NORTH)


@dataclass
class TerrainMesh:
    nx: int
    ny: int
    nz: int
    cellsize: float
    points: np.ndarray            # (Npts, 3) float
    faces: np.ndarray             # (Nfaces, 4) int (all faces are quads)
    owner: np.ndarray             # (Nfaces,) int
    neighbour: np.ndarray         # (nInternalFaces,) int
    patches: dict[str, tuple[int, int]]  # name -> (startFace, nFaces)
    layer_heights: np.ndarray     # (nz+1,) interface heights above ground (m)

    @property
    def cell_agl(self) -> np.ndarray:
        """Height above ground (m) of each vertical cell centre, shape (nz,)."""
        h = self.layer_heights
        return 0.5 * (h[:-1] + h[1:])

    @property
    def n_cells(self) -> int:
        return self.nx * self.ny * self.nz

    @property
    def n_internal_faces(self) -> int:
        return len(self.neighbour)

    def column_of(self, i: int, j: int) -> tuple[int, int]:
        """Map mesh column (i, j) back to DEM (row, col)."""
        return self.ny - 1 - j, i

    # --- writing ----------------------------------------------------------
    def write(self, case_dir: str) -> None:
        """Write ``<case_dir>/constant/polyMesh`` files."""
        import os
        pm = os.path.join(case_dir, "constant", "polyMesh")
        os.makedirs(pm, exist_ok=True)
        loc = "constant/polyMesh"

        n_pts, n_faces = len(self.points), len(self.faces)
        n_int = self.n_internal_faces
        note = (f'"nPoints:{n_pts} nCells:{self.n_cells} '
                f'nFaces:{n_faces} nInternalFaces:{n_int}"')

        with open(os.path.join(pm, "points"), "w") as fh:
            fh.write(dictwriter.header("vectorField", "points", loc))
            fh.write(f"\n{n_pts}\n(\n")
            fh.write("\n".join(f"({x:.6f} {y:.6f} {z:.6f})"
                               for x, y, z in self.points))
            fh.write("\n)\n")

        with open(os.path.join(pm, "faces"), "w") as fh:
            fh.write(dictwriter.header("faceList", "faces", loc))
            fh.write(f"\n{n_faces}\n(\n")
            fh.write("\n".join(f"4({a} {b} {c} {d})"
                               for a, b, c, d in self.faces))
            fh.write("\n)\n")

        for name, data in (("owner", self.owner), ("neighbour", self.neighbour)):
            with open(os.path.join(pm, name), "w") as fh:
                hdr = dictwriter.header("labelList", name, loc)
                # Insert the mesh note into the FoamFile block.
                hdr = hdr.replace("    object", f"    note        {note};\n    object")
                fh.write(hdr)
                fh.write(f"\n{len(data)}\n(\n")
                fh.write("\n".join(str(int(v)) for v in data))
                fh.write("\n)\n")

        boundary = {}
        for name in PATCHES:
            start, count = self.patches[name]
            boundary[name] = {
                "type": "wall" if name == GROUND else "patch",
                "nFaces": count,
                "startFace": start,
            }
        with open(os.path.join(pm, "boundary"), "w") as fh:
            fh.write(dictwriter.header("polyBoundaryMesh", "boundary", loc))
            fh.write(f"\n{len(PATCHES)}\n(\n")
            fh.write(dictwriter.render(boundary, level=0))
            fh.write("\n)\n")


def _node_elevations(elev_north_up: np.ndarray) -> np.ndarray:
    """Cell-centred DEM -> node-cornered elevations, shape (ny+1, nx+1)."""
    p = np.pad(elev_north_up, 1, mode="edge")
    return 0.25 * (p[:-1, :-1] + p[1:, :-1] + p[:-1, 1:] + p[1:, 1:])


def _layer_heights(domain_height: float, nz: int, expansion: float) -> np.ndarray:
    """Heights above ground of the nz+1 layer interfaces, geometric grading."""
    if abs(expansion - 1.0) < 1e-9:
        thick = np.ones(nz)
    else:
        thick = expansion ** np.arange(nz)
    h = np.concatenate([[0.0], np.cumsum(thick)])
    return h / h[-1] * domain_height


def build_terrain_mesh(
    elevation: np.ndarray,
    cellsize: float,
    *,
    nz: int = 20,
    domain_height: float = 1000.0,
    expansion_ratio: float = 1.2,
) -> TerrainMesh:
    """Build a terrain-following hex mesh from a DEM.

    ``elevation`` is the DEM in metres (row 0 = north). ``domain_height`` is the
    depth of the domain above the ground surface (m); ``nz`` vertical cells are
    geometrically graded by ``expansion_ratio`` (fine near the ground).
    """
    elev = np.asarray(elevation, dtype=float)
    ny, nx = elev.shape
    elev_nu = elev[::-1]  # flip so mesh j increases north
    node_elev = _node_elevations(elev_nu)            # (ny+1, nx+1)
    h = _layer_heights(domain_height, nz, expansion_ratio)  # (nz+1,)

    nPi, nPj = nx + 1, ny + 1

    # Points, ordered i-fastest then j then k (matches nidx below).
    kk, jj, ii = np.meshgrid(np.arange(nz + 1), np.arange(nPj), np.arange(nPi),
                             indexing="ij")
    x = ii * cellsize
    y = jj * cellsize
    z = node_elev[jj, ii] + h[kk]
    points = np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)

    def nidx(i, j, k):
        return k * nPj * nPi + j * nPi + i

    def cidx(i, j, k):
        return k * ny * nx + j * nx + i

    int_owner, int_neigh, int_quads = [], [], []
    patch_quads: dict[str, list] = {p: [] for p in PATCHES}
    patch_owner: dict[str, list] = {p: [] for p in PATCHES}

    # --- x-normal faces (planes I = 0..nx) --------------------------------
    J, K = np.meshgrid(np.arange(ny), np.arange(nz), indexing="ij")
    J, K = J.ravel(), K.ravel()
    for I in range(nx + 1):
        quad = np.stack([nidx(I, J, K), nidx(I, J + 1, K),
                         nidx(I, J + 1, K + 1), nidx(I, J, K + 1)], axis=1)
        if I == 0:                       # west boundary, outward -x
            patch_quads[WEST].append(quad[:, ::-1])
            patch_owner[WEST].append(cidx(0, J, K))
        elif I == nx:                    # east boundary, outward +x
            patch_quads[EAST].append(quad)
            patch_owner[EAST].append(cidx(nx - 1, J, K))
        else:
            int_quads.append(quad)
            int_owner.append(cidx(I - 1, J, K))
            int_neigh.append(cidx(I, J, K))

    # --- y-normal faces (planes J = 0..ny) --------------------------------
    I, K = np.meshgrid(np.arange(nx), np.arange(nz), indexing="ij")
    I, K = I.ravel(), K.ravel()
    for Jp in range(ny + 1):
        quad = np.stack([nidx(I, Jp, K), nidx(I, Jp, K + 1),
                         nidx(I + 1, Jp, K + 1), nidx(I + 1, Jp, K)], axis=1)
        if Jp == 0:                      # south boundary, outward -y
            patch_quads[SOUTH].append(quad[:, ::-1])
            patch_owner[SOUTH].append(cidx(I, 0, K))
        elif Jp == ny:                   # north boundary, outward +y
            patch_quads[NORTH].append(quad)
            patch_owner[NORTH].append(cidx(I, ny - 1, K))
        else:
            int_quads.append(quad)
            int_owner.append(cidx(I, Jp - 1, K))
            int_neigh.append(cidx(I, Jp, K))

    # --- z-normal faces (planes K = 0..nz) --------------------------------
    I, J = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    I, J = I.ravel(), J.ravel()
    for Kp in range(nz + 1):
        quad = np.stack([nidx(I, J, Kp), nidx(I + 1, J, Kp),
                         nidx(I + 1, J + 1, Kp), nidx(I, J + 1, Kp)], axis=1)
        if Kp == 0:                      # ground, outward -z
            patch_quads[GROUND].append(quad[:, ::-1])
            patch_owner[GROUND].append(cidx(I, J, 0))
        elif Kp == nz:                   # top, outward +z
            patch_quads[TOP].append(quad)
            patch_owner[TOP].append(cidx(I, J, nz - 1))
        else:
            int_quads.append(quad)
            int_owner.append(cidx(I, J, Kp - 1))
            int_neigh.append(cidx(I, J, Kp))

    int_quads = np.concatenate(int_quads)
    int_owner = np.concatenate(int_owner)
    int_neigh = np.concatenate(int_neigh)
    order = np.lexsort((int_neigh, int_owner))  # upper-triangular ordering
    int_quads, int_owner, int_neigh = (
        int_quads[order], int_owner[order], int_neigh[order])

    faces = [int_quads]
    owners = [int_owner]
    patches: dict[str, tuple[int, int]] = {}
    start = len(int_owner)
    for name in PATCHES:
        q = np.concatenate(patch_quads[name])
        o = np.concatenate(patch_owner[name])
        faces.append(q)
        owners.append(o)
        patches[name] = (start, len(o))
        start += len(o)

    return TerrainMesh(
        nx=nx, ny=ny, nz=nz, cellsize=cellsize,
        points=points,
        faces=np.concatenate(faces),
        owner=np.concatenate(owners),
        neighbour=int_neigh,
        patches=patches,
        layer_heights=h,
    )
