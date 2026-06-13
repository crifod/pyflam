"""Spatially-varying wind fields.

A :class:`WindField` is the shared contract between pyflam's wind models and the
fire model: gridded wind speed + direction that feeds
:func:`pyflam.landscape.basic_fire_behavior` through a per-cell ``wind_midflame``
array. Two solvers produce one:

* :mod:`pyflam.windsolver` — fast diagnostic mass-consistent model (no deps);
* :mod:`pyflam.cfd` — momentum/RANS model via OpenFOAM (terrain + stability +
  diurnal slope flows).

This module also provides a generic ESRI/Arc ASCII grid reader
(:func:`read_esri_ascii`), handy for ingesting gridded data from other tools.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# Wind-speed unit -> ft/min (midflame wind feeds Rothermel in ft/min).
_SPEED_TO_FT_MIN = {
    "mph": 88.0,
    "mps": 1.0 / 0.3048 * 60.0,   # m/s
    "m/s": 1.0 / 0.3048 * 60.0,
    "kph": 1000.0 / 0.3048 / 60.0,
    "km/h": 1000.0 / 0.3048 / 60.0,
    "kmph": 1000.0 / 0.3048 / 60.0,
    "kts": 6076.12 / 60.0,        # knots
    "knots": 6076.12 / 60.0,
    "ftmin": 1.0,
    "ft/min": 1.0,
}


@dataclass
class EsriGrid:
    """A parsed ESRI/Arc ASCII raster (row 0 = north)."""

    data: np.ndarray
    cellsize: float
    xllcorner: float
    yllcorner: float
    nodata: float

    @property
    def west(self) -> float:
        return self.xllcorner

    @property
    def north(self) -> float:
        return self.yllcorner + self.data.shape[0] * self.cellsize


def read_esri_ascii(path: str) -> EsriGrid:
    """Read an ESRI ASCII (``.asc``) raster."""
    header: dict[str, float] = {}
    keys = {"ncols", "nrows", "xllcorner", "yllcorner",
            "xllcenter", "yllcenter", "cellsize", "nodata_value"}
    with open(path) as fh:
        values: list[float] = []
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            key = parts[0].lower()
            if key in keys and len(parts) >= 2:
                header[key] = float(parts[1])
            else:
                values.extend(float(v) for v in parts)

    nrows = int(header["nrows"])
    ncols = int(header["ncols"])
    cellsize = header["cellsize"]
    if "xllcorner" in header:
        xll, yll = header["xllcorner"], header["yllcorner"]
    else:
        xll = header["xllcenter"] - cellsize / 2.0
        yll = header["yllcenter"] - cellsize / 2.0
    data = np.array(values, dtype=float).reshape(nrows, ncols)
    return EsriGrid(data=data, cellsize=cellsize, xllcorner=xll, yllcorner=yll,
                    nodata=header.get("nodata_value", -9999.0))


@dataclass
class WindField:
    """A gridded wind field (north-up, row 0 = north).

    ``speed`` is in ``speed_units``; ``direction`` is degrees clockwise from
    north giving the direction the wind blows *from* (meteorological convention,
    matching FlamMap's wind-direction input).
    """

    speed: np.ndarray
    direction: np.ndarray  # degrees FROM, met convention
    cellsize: float
    west: float
    north: float
    speed_units: str = "mph"
    height: float = 20.0
    height_units: str = "ft"
    crs: object | None = None

    @property
    def shape(self) -> tuple[int, int]:
        return self.speed.shape

    @property
    def direction_toward(self) -> np.ndarray:
        """Direction the wind blows *toward* (degrees from north)."""
        return (self.direction + 180.0) % 360.0

    def speed_ft_per_min(self) -> np.ndarray:
        try:
            factor = _SPEED_TO_FT_MIN[self.speed_units.lower()]
        except KeyError:
            raise ValueError(
                f"unknown wind speed unit {self.speed_units!r}; "
                f"expected one of {sorted(_SPEED_TO_FT_MIN)}"
            ) from None
        return np.asarray(self.speed, dtype=float) * factor

    def to_landscape(self, ls) -> "WindField":
        """Nearest-neighbour resample onto a landscape's grid.

        Angles can't be linearly averaged, and wind grids are usually coarser
        than the fuel grid, so nearest-neighbour is the safe, simple choice.
        """
        nrows, ncols = ls.shape
        cols = np.arange(ncols)
        rows = np.arange(nrows)
        x = ls.west + (cols + 0.5) * ls.cellsize_x          # cell-centre x
        y = ls.north - (rows + 0.5) * ls.cellsize_y          # cell-centre y
        wc = np.clip(((x - self.west) / self.cellsize).astype(int),
                     0, self.shape[1] - 1)
        wr = np.clip(((self.north - y) / self.cellsize).astype(int),
                     0, self.shape[0] - 1)
        rr, cc = np.meshgrid(wr, wc, indexing="ij")
        return WindField(
            speed=self.speed[rr, cc], direction=self.direction[rr, cc],
            cellsize=ls.cellsize_x, west=ls.west, north=ls.north,
            speed_units=self.speed_units, height=self.height,
            height_units=self.height_units, crs=ls.crs,
        )

    def to_midflame(self, ls=None, *, wind_reduction_factor: float = 0.4) -> np.ndarray:
        """Midflame wind speed (ft/min) for the surface model.

        Multiplies this field's wind by a wind reduction factor (WRF) to get
        midflame wind. The default 0.4 is a reasonable unsheltered value; a
        canopy-derived per-cell WRF is a later refinement. If ``ls`` is given,
        the field is first resampled onto its grid.
        """
        wf = self.to_landscape(ls) if ls is not None else self
        return wf.speed_ft_per_min() * wind_reduction_factor
