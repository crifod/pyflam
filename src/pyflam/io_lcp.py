"""Read and write FlamMap/FARSITE landscape (``.lcp``) files.

An LCP is a flat binary raster stack: a fixed 7316-byte header followed by the
cell data as little-endian ``int16``, band-interleaved-by-pixel (all bands for a
cell stored consecutively), row-major from the north-west corner.

Bands, in file order:
    0 elevation, 1 slope, 2 aspect, 3 fuel model, 4 canopy cover,
    then (if crown fuels present) 5 canopy height, 6 canopy base height,
    7 canopy bulk density,
    then (if ground fuels present) duff, coarse woody debris.

The header advertises crown/ground fuel presence with two flags (20 = absent,
21 = present), which is how the band count is determined.

This is a pure-Python implementation (no GDAL needed). It is validated in the
test suite against GDAL's own LCP driver. Only the fields needed to reconstruct
geometry and the band stack are written; the per-band value dictionaries are
filled minimally (lo/hi + ``num = -1``) and the trailing file-name/description
fields are left zeroed, which FlamMap and GDAL both accept.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass

import numpy as np

HEADER_SIZE = 7316
_CROWN_PRESENT = 21
_GROUND_PRESENT = 21
_CROWN_ABSENT = 20
_GROUND_ABSENT = 20

# Header field offsets (bytes). Only the geometry-critical fields are named;
# everything else in the header is left at its default (zero) on write.
_OFF_CROWN = 0      # int32
_OFF_GROUND = 4     # int32
_OFF_LATITUDE = 8   # int32
_OFF_VALUE_BLOCKS = 44   # 10 blocks x (lo,hi,num int32 + 100 int32 values) = 4120
_OFF_NUMEAST = 4164   # int32 (columns)
_OFF_NUMNORTH = 4168  # int32 (rows)
_OFF_EASTUTM = 4172   # double (max x)
_OFF_WESTUTM = 4180   # double (min x)
_OFF_NORTHUTM = 4188  # double (max y)
_OFF_SOUTHUTM = 4196  # double (min y)
_OFF_GRIDUNITS = 4204  # double (0 = metres, 1 = feet, 2 = km)
_OFF_XRESOL = 4212    # double
_OFF_YRESOL = 4220    # double
_OFF_SUNITS = 4230    # int16 (slope units: 0 = degrees, 1 = percent)

# Standard band names in file order.
BASE_BANDS = ("elevation", "slope", "aspect", "fuel_model", "canopy_cover")
CROWN_BANDS = ("canopy_height", "canopy_base_height", "canopy_bulk_density")
GROUND_BANDS = ("duff", "coarse_woody")


@dataclass
class LcpData:
    """Raw contents of an LCP file."""

    bands: dict[str, np.ndarray]  # name -> 2D int16 array (north-up, row 0 = north)
    ncols: int
    nrows: int
    west: float
    north: float
    east: float
    south: float
    cellsize_x: float
    cellsize_y: float
    slope_units: str  # "degrees" or "percent"
    latitude: int


def read_lcp(path: str) -> LcpData:
    """Read an ``.lcp`` file into its bands and geometry."""
    with open(path, "rb") as fh:
        raw = fh.read()
    if len(raw) < HEADER_SIZE:
        raise ValueError(f"{path}: too short to be an LCP ({len(raw)} bytes)")

    crown = struct.unpack_from("<i", raw, _OFF_CROWN)[0]
    ground = struct.unpack_from("<i", raw, _OFF_GROUND)[0]
    n_band_names = list(BASE_BANDS)
    if crown == _CROWN_PRESENT:
        n_band_names += list(CROWN_BANDS)
    if ground == _GROUND_PRESENT:
        n_band_names += list(GROUND_BANDS)
    nbands = len(n_band_names)

    ncols = struct.unpack_from("<i", raw, _OFF_NUMEAST)[0]
    nrows = struct.unpack_from("<i", raw, _OFF_NUMNORTH)[0]
    east = struct.unpack_from("<d", raw, _OFF_EASTUTM)[0]
    west = struct.unpack_from("<d", raw, _OFF_WESTUTM)[0]
    north = struct.unpack_from("<d", raw, _OFF_NORTHUTM)[0]
    south = struct.unpack_from("<d", raw, _OFF_SOUTHUTM)[0]
    xres = struct.unpack_from("<d", raw, _OFF_XRESOL)[0]
    yres = struct.unpack_from("<d", raw, _OFF_YRESOL)[0]
    sunits = struct.unpack_from("<h", raw, _OFF_SUNITS)[0]

    expected = nrows * ncols * nbands * 2
    body = raw[HEADER_SIZE:HEADER_SIZE + expected]
    if len(body) < expected:
        raise ValueError(
            f"{path}: truncated data (need {expected} bytes, have {len(body)})"
        )

    stack = np.frombuffer(body, dtype="<i2").reshape(nrows, ncols, nbands)
    bands = {name: stack[:, :, i].copy() for i, name in enumerate(n_band_names)}

    return LcpData(
        bands=bands, ncols=ncols, nrows=nrows,
        west=west, north=north, east=east, south=south,
        cellsize_x=xres, cellsize_y=yres,
        slope_units="percent" if sunits == 1 else "degrees",
        latitude=struct.unpack_from("<i", raw, _OFF_LATITUDE)[0],
    )


def write_lcp(path: str, lcp: LcpData) -> None:
    """Write an :class:`LcpData` to an ``.lcp`` file."""
    names = list(lcp.bands.keys())
    has_crown = all(b in lcp.bands for b in CROWN_BANDS)
    has_ground = all(b in lcp.bands for b in GROUND_BANDS)

    # Bands must appear in canonical file order.
    order = list(BASE_BANDS)
    if has_crown:
        order += list(CROWN_BANDS)
    if has_ground:
        order += list(GROUND_BANDS)
    missing = [b for b in BASE_BANDS if b not in lcp.bands]
    if missing:
        raise ValueError(f"LCP requires base bands; missing {missing}")

    header = bytearray(HEADER_SIZE)
    struct.pack_into("<i", header, _OFF_CROWN,
                     _CROWN_PRESENT if has_crown else _CROWN_ABSENT)
    struct.pack_into("<i", header, _OFF_GROUND,
                     _GROUND_PRESENT if has_ground else _GROUND_ABSENT)
    struct.pack_into("<i", header, _OFF_LATITUDE, int(lcp.latitude))

    # Per-band value dictionaries: lo, hi, num=-1, 100 zeros.
    for i, name in enumerate(order):
        arr = lcp.bands[name]
        lo = int(arr.min())
        hi = int(arr.max())
        base = _OFF_VALUE_BLOCKS + i * 412
        struct.pack_into("<iii", header, base, lo, hi, -1)

    struct.pack_into("<i", header, _OFF_NUMEAST, lcp.ncols)
    struct.pack_into("<i", header, _OFF_NUMNORTH, lcp.nrows)
    struct.pack_into("<d", header, _OFF_EASTUTM, lcp.east)
    struct.pack_into("<d", header, _OFF_WESTUTM, lcp.west)
    struct.pack_into("<d", header, _OFF_NORTHUTM, lcp.north)
    struct.pack_into("<d", header, _OFF_SOUTHUTM, lcp.south)
    struct.pack_into("<d", header, _OFF_GRIDUNITS, 0.0)  # metres
    struct.pack_into("<d", header, _OFF_XRESOL, lcp.cellsize_x)
    struct.pack_into("<d", header, _OFF_YRESOL, lcp.cellsize_y)
    struct.pack_into("<h", header, _OFF_SUNITS,
                     1 if lcp.slope_units == "percent" else 0)

    stack = np.stack(
        [lcp.bands[name].astype("<i2") for name in order], axis=-1
    )
    with open(path, "wb") as fh:
        fh.write(header)
        fh.write(stack.tobytes())
