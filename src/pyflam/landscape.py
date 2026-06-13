"""Gridded landscapes and whole-landscape surface fire behavior.

A :class:`Landscape` is the in-memory band stack (fuel model, slope, and the
other LANDFIRE/LCP layers) plus georeferencing. It can be built from an ``.lcp``
file (pure Python, see :mod:`pyflam.io_lcp`) or from GeoTIFFs (via rasterio).

:func:`basic_fire_behavior` runs the Rothermel surface model over every cell,
the per-cell equivalent of FlamMap's "Basic Fire Behavior" outputs.

Step-2 scope and caveats:
    * Fuel moisture and wind are uniform across the landscape; the fuel model and
      slope vary per cell.
    * Wind and slope are combined as Rothermel's scalar ``1 + phi_w + phi_s``
      (i.e. wind blowing upslope, the maximum-spread case). Directional spread
      using aspect + wind vector is step 3.
    * Outputs are in English units (ft/min, Btu/ft/s, ft); convert with
      :mod:`pyflam.units`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import fuel_models, io_lcp
from .rothermel import kernel_param_groups, surface_kernel


@dataclass
class Landscape:
    """A georeferenced raster stack of fuel/terrain layers.

    Arrays are 2D, north-up (row 0 is the northern edge). Only ``fuel_model`` and
    ``slope`` are required for surface fire behavior; the rest are carried for
    later roadmap steps (crown fire, wind reduction) and for round-tripping LCPs.
    """

    fuel_model: np.ndarray
    slope: np.ndarray
    cellsize_x: float
    cellsize_y: float
    west: float
    north: float
    elevation: np.ndarray | None = None
    aspect: np.ndarray | None = None
    canopy_cover: np.ndarray | None = None
    canopy_height: np.ndarray | None = None
    canopy_base_height: np.ndarray | None = None
    canopy_bulk_density: np.ndarray | None = None
    slope_units: str = "degrees"  # "degrees" | "percent" | "tangent"
    crs: object | None = None

    @property
    def shape(self) -> tuple[int, int]:
        return self.fuel_model.shape

    @property
    def transform_gdal(self) -> tuple[float, float, float, float, float, float]:
        """GDAL-style geotransform (origin top-left, north-up)."""
        return (self.west, self.cellsize_x, 0.0,
                self.north, 0.0, -self.cellsize_y)

    @property
    def slope_tangent(self) -> np.ndarray:
        """Slope as rise/run (the form the Rothermel slope factor needs)."""
        s = np.asarray(self.slope, dtype=float)
        if self.slope_units == "degrees":
            return np.tan(np.radians(np.clip(s, 0.0, 89.9)))
        if self.slope_units == "percent":
            return np.clip(s, 0.0, None) / 100.0
        return np.clip(s, 0.0, None)

    # --- LCP I/O -----------------------------------------------------------
    @classmethod
    def from_lcp(cls, path: str) -> "Landscape":
        return cls.from_lcpdata(io_lcp.read_lcp(path))

    @classmethod
    def from_lcpdata(cls, lcp: io_lcp.LcpData) -> "Landscape":
        b = lcp.bands
        return cls(
            fuel_model=b["fuel_model"], slope=b["slope"],
            elevation=b.get("elevation"), aspect=b.get("aspect"),
            canopy_cover=b.get("canopy_cover"),
            canopy_height=b.get("canopy_height"),
            canopy_base_height=b.get("canopy_base_height"),
            canopy_bulk_density=b.get("canopy_bulk_density"),
            cellsize_x=lcp.cellsize_x, cellsize_y=lcp.cellsize_y,
            west=lcp.west, north=lcp.north, slope_units=lcp.slope_units,
        )

    def to_lcpdata(self, latitude: int = 0) -> io_lcp.LcpData:
        nrows, ncols = self.shape
        zeros = np.zeros(self.shape, dtype="<i2")
        bands = {
            "elevation": _as_i16(self.elevation, zeros),
            "slope": _as_i16(self.slope, zeros),
            "aspect": _as_i16(self.aspect, zeros),
            "fuel_model": _as_i16(self.fuel_model, zeros),
            "canopy_cover": _as_i16(self.canopy_cover, zeros),
        }
        crown = (self.canopy_height, self.canopy_base_height,
                 self.canopy_bulk_density)
        if all(c is not None for c in crown):
            bands["canopy_height"] = _as_i16(self.canopy_height, zeros)
            bands["canopy_base_height"] = _as_i16(self.canopy_base_height, zeros)
            bands["canopy_bulk_density"] = _as_i16(self.canopy_bulk_density, zeros)
        return io_lcp.LcpData(
            bands=bands, ncols=ncols, nrows=nrows,
            west=self.west, north=self.north,
            east=self.west + ncols * self.cellsize_x,
            south=self.north - nrows * self.cellsize_y,
            cellsize_x=self.cellsize_x, cellsize_y=self.cellsize_y,
            slope_units=("percent" if self.slope_units == "percent"
                         else "degrees"),
            latitude=latitude,
        )

    def to_lcp(self, path: str, latitude: int = 0) -> None:
        io_lcp.write_lcp(path, self.to_lcpdata(latitude=latitude))

    # --- GeoTIFF I/O (rasterio) -------------------------------------------
    @classmethod
    def from_geotiffs(cls, paths: dict[str, str], *,
                      slope_units: str = "degrees") -> "Landscape":
        """Build a landscape from one GeoTIFF per band.

        ``paths`` maps band names (``"fuel_model"``, ``"slope"``, ...) to file
        paths. ``fuel_model`` and ``slope`` are required; all rasters must share
        the same grid (shape and transform).
        """
        rasterio = _require_rasterio()
        for required in ("fuel_model", "slope"):
            if required not in paths:
                raise ValueError(f"from_geotiffs needs a {required!r} band")

        arrays: dict[str, np.ndarray] = {}
        transform = crs = ref_shape = None
        for name, path in paths.items():
            with rasterio.open(path) as ds:
                arrays[name] = ds.read(1)
                if ref_shape is None:
                    ref_shape, transform, crs = arrays[name].shape, ds.transform, ds.crs
                elif arrays[name].shape != ref_shape:
                    raise ValueError(
                        f"band {name!r} shape {arrays[name].shape} != {ref_shape}"
                    )
        return cls(
            fuel_model=arrays["fuel_model"], slope=arrays["slope"],
            elevation=arrays.get("elevation"), aspect=arrays.get("aspect"),
            canopy_cover=arrays.get("canopy_cover"),
            canopy_height=arrays.get("canopy_height"),
            canopy_base_height=arrays.get("canopy_base_height"),
            canopy_bulk_density=arrays.get("canopy_bulk_density"),
            cellsize_x=transform.a, cellsize_y=-transform.e,
            west=transform.c, north=transform.f,
            slope_units=slope_units, crs=crs,
        )

    def to_geotiff(self, path: str, array: np.ndarray, *,
                   nodata: float = np.nan, dtype: str = "float32") -> None:
        """Write one 2D array on this landscape's grid to a GeoTIFF."""
        rasterio = _require_rasterio()
        from rasterio.transform import Affine

        if array.shape != self.shape:
            raise ValueError(f"array shape {array.shape} != landscape {self.shape}")
        transform = Affine(self.cellsize_x, 0.0, self.west,
                           0.0, -self.cellsize_y, self.north)
        with rasterio.open(
            path, "w", driver="GTiff",
            height=self.shape[0], width=self.shape[1], count=1,
            dtype=dtype, crs=self.crs, transform=transform, nodata=nodata,
        ) as ds:
            ds.write(array.astype(dtype), 1)


def basic_fire_behavior(
    ls: Landscape,
    *,
    m_1h: float,
    m_10h: float,
    m_100h: float,
    m_live_herb: float = 0.0,
    m_live_woody: float = 0.0,
    wind_midflame: float | np.ndarray = 0.0,
    load_factor: float = 1.0,
    nodata: float = np.nan,
) -> dict[str, np.ndarray]:
    """Per-cell surface fire behavior over a whole landscape.

    ``wind_midflame`` (ft/min) may be a single value or a 2D array the same shape
    as the landscape — pass a spatially-varying field (from
    :mod:`pyflam.windsolver` or :mod:`pyflam.cfd`) to capture terrain-driven wind
    across the grid.

    Returns float arrays (same shape as the landscape) for ``rate_of_spread``
    (ft/min), ``fireline_intensity`` (Btu/ft/s), ``flame_length`` (ft) and
    ``reaction_intensity`` (Btu/ft^2/min). Nonburnable cells are 0; cells with an
    unrecognized fuel-model code are set to ``nodata``.

    Efficiency: the expensive Rothermel terms are computed once per *unique* fuel
    model (there are at most a few dozen), then the cheap per-cell wind and slope
    factors are applied vectorized across all cells of that model.
    """
    shape = ls.shape
    ros = np.full(shape, nodata, dtype=float)
    fli = np.full(shape, nodata, dtype=float)
    fl = np.full(shape, nodata, dtype=float)
    ri = np.full(shape, nodata, dtype=float)

    tan_slope = ls.slope_tangent
    fuel = np.asarray(ls.fuel_model)
    # Wind may be a single value or a per-cell field (from the wind solvers).
    wind = np.broadcast_to(np.asarray(wind_midflame, dtype=float), shape)
    moist = dict(m_1h=m_1h, m_10h=m_10h, m_100h=m_100h,
                 m_live_herb=m_live_herb, m_live_woody=m_live_woody)

    for num in np.unique(fuel):
        num = int(num)
        mask = fuel == num
        try:
            fm = fuel_models.get(num)
        except KeyError:
            continue  # leave as nodata (e.g. -9999 fill or undefined code)

        if not fm.is_burnable:
            ros[mask] = fli[mask] = fl[mask] = ri[mask] = 0.0
            continue

        # Heavy fuel/moisture terms computed once per (fuel, load factor, moisture
        # bin); wind & slope applied per cell. load_factor and the moistures may be
        # scalars or per-cell fields (e.g. spatial weather).
        lf = load_factor if not isinstance(load_factor, dict) \
            else load_factor.get(num, 1.0)
        for sub, p in kernel_param_groups(mask, {"load_factor": lf, **moist}):
            kernel = surface_kernel(fm, **p)
            r = kernel.rate_of_spread(wind[sub], tan_slope[sub])
            i_byram = kernel.heat_per_unit_area * r / 60.0
            ros[sub] = r
            ri[sub] = kernel.reaction_intensity
            fli[sub] = i_byram
            fl[sub] = np.where(i_byram > 0.0, 0.45 * i_byram ** 0.46, 0.0)

    return {
        "rate_of_spread": ros,
        "fireline_intensity": fli,
        "flame_length": fl,
        "reaction_intensity": ri,
    }


def _as_i16(arr, default):
    return default if arr is None else np.asarray(arr).astype("<i2")


def _require_rasterio():
    try:
        import rasterio
    except ImportError as exc:  # pragma: no cover - exercised only without rasterio
        raise ImportError(
            "GeoTIFF I/O needs rasterio; install with: pip install 'pyflam[geo]'"
        ) from exc
    return rasterio
