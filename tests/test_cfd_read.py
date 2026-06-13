"""Tests for the OpenFOAM result reader (no OpenFOAM; fabricated fields)."""

from __future__ import annotations

import numpy as np
import pytest

from pyflam.cfd import mesh as M
from pyflam.cfd import read as R
from pyflam.cfd.dictwriter import header


def _write_u(path, vectors):
    with open(path, "w") as fh:
        fh.write(header("volVectorField", "U", str(path.parent.name)))
        fh.write(f"\ndimensions      [0 1 -1 0 0 0 0];\n")
        fh.write(f"internalField   nonuniform List<vector>\n{len(vectors)}\n(\n")
        fh.write("\n".join(f"({x} {y} {z})" for x, y, z in vectors))
        fh.write("\n)\n;\n")


def test_read_internal_vector_uniform(tmp_path):
    p = tmp_path / "U"
    p.write_text(header("volVectorField", "U") +
                 "\ndimensions [0 1 -1 0 0 0 0];\ninternalField uniform (3 4 0);\n")
    arr = R.read_internal_vector(str(p), 5)
    assert arr.shape == (5, 3)
    assert np.allclose(arr, [3, 4, 0])


def test_read_internal_vector_nonuniform(tmp_path):
    p = tmp_path / "U"
    _write_u(p, [(1, 0, 0), (2, 0, 0), (3, 0, 0)])
    arr = R.read_internal_vector(str(p), 3)
    assert np.allclose(arr[:, 0], [1, 2, 3])


def test_read_wind_field_uniform_flow(tmp_path):
    # 2x2x2 mesh; uniform east wind everywhere -> direction "from west" (270).
    m = M.build_terrain_mesh(np.full((2, 2), 1000.0), 30.0, nz=2,
                             domain_height=100.0)
    t = tmp_path / "100"
    t.mkdir()
    _write_u(t / "U", [(5.0, 0.0, 0.0)] * m.n_cells)
    wf = R.read_wind_field(str(tmp_path), m, output_height=m.cell_agl[0],
                           z0=0.1, west=0.0, north=60.0)
    assert wf.shape == (2, 2)
    assert np.allclose(wf.speed, 5.0)
    assert np.allclose(wf.direction, 270.0)
    assert wf.speed_units == "m/s"


def test_latest_time_dir(tmp_path):
    for t in ("0", "100", "2000", "constant"):
        (tmp_path / t).mkdir()
    assert R.latest_time_dir(str(tmp_path)) == "2000"
