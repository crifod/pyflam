"""Tests for the pure-Python LCP reader/writer.

The strongest check here is cross-validation against GDAL's own LCP driver
(via rasterio): we write an LCP with pyflam and confirm GDAL reads back the same
bands and geometry. That tests our writer against an independent implementation,
not just against our own reader.
"""

from __future__ import annotations

import numpy as np
import pytest

from pyflam import io_lcp
from pyflam.io_lcp import HEADER_SIZE, LcpData, read_lcp, write_lcp


def _sample(nrows=4, ncols=5, crown=False):
    rng = np.random.default_rng(0)
    bands = {
        "elevation": rng.integers(800, 2000, (nrows, ncols)).astype("<i2"),
        "slope": rng.integers(0, 45, (nrows, ncols)).astype("<i2"),
        "aspect": rng.integers(0, 360, (nrows, ncols)).astype("<i2"),
        "fuel_model": rng.choice([1, 2, 101, 102, 165, 91], (nrows, ncols)).astype("<i2"),
        "canopy_cover": rng.integers(0, 100, (nrows, ncols)).astype("<i2"),
    }
    if crown:
        bands["canopy_height"] = rng.integers(0, 400, (nrows, ncols)).astype("<i2")
        bands["canopy_base_height"] = rng.integers(0, 100, (nrows, ncols)).astype("<i2")
        bands["canopy_bulk_density"] = rng.integers(0, 40, (nrows, ncols)).astype("<i2")
    return LcpData(
        bands=bands, ncols=ncols, nrows=nrows,
        west=500000.0, north=5000000.0,
        east=500000.0 + ncols * 30.0, south=5000000.0 - nrows * 30.0,
        cellsize_x=30.0, cellsize_y=30.0, slope_units="degrees", latitude=46,
    )


def test_header_size_and_file_length(tmp_path):
    lcp = _sample()
    path = tmp_path / "x.lcp"
    write_lcp(str(path), lcp)
    expected = HEADER_SIZE + lcp.nrows * lcp.ncols * 5 * 2
    assert path.stat().st_size == expected


@pytest.mark.parametrize("crown", [False, True])
def test_roundtrip_bands_and_geometry(tmp_path, crown):
    lcp = _sample(crown=crown)
    path = tmp_path / "x.lcp"
    write_lcp(str(path), lcp)
    back = read_lcp(str(path))

    assert set(back.bands) == set(lcp.bands)
    for name, arr in lcp.bands.items():
        assert np.array_equal(back.bands[name], arr), name
    assert (back.ncols, back.nrows) == (lcp.ncols, lcp.nrows)
    assert back.cellsize_x == lcp.cellsize_x
    assert back.west == lcp.west and back.north == lcp.north
    assert back.slope_units == "degrees"


def test_band_count_flag(tmp_path):
    """Crown bands must round-trip via the crown-fuels presence flag."""
    base = tmp_path / "base.lcp"
    crowned = tmp_path / "crown.lcp"
    write_lcp(str(base), _sample(crown=False))
    write_lcp(str(crowned), _sample(crown=True))
    assert len(read_lcp(str(base)).bands) == 5
    assert len(read_lcp(str(crowned)).bands) == 8


def test_truncated_file_raises(tmp_path):
    path = tmp_path / "short.lcp"
    path.write_bytes(b"\x00" * 100)
    with pytest.raises(ValueError):
        read_lcp(str(path))


# --- Cross-validation against GDAL's LCP driver -------------------------------

def test_gdal_reads_our_lcp(tmp_path):
    rasterio = pytest.importorskip("rasterio")
    lcp = _sample(nrows=6, ncols=7, crown=True)
    path = tmp_path / "x.lcp"
    write_lcp(str(path), lcp)

    with rasterio.open(str(path)) as ds:
        assert ds.driver == "LCP"
        assert (ds.width, ds.height) == (lcp.ncols, lcp.nrows)
        assert ds.count == 8  # 5 base + 3 crown
        order = list(io_lcp.BASE_BANDS) + list(io_lcp.CROWN_BANDS)
        for i, name in enumerate(order, start=1):
            assert np.array_equal(
                ds.read(i).astype("<i2"), lcp.bands[name]
            ), f"GDAL band {i} ({name}) mismatch"
