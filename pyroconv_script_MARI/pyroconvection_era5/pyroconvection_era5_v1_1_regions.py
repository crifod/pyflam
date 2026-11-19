#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pyroconvection_era5.py
======================

Mappa la predisposizione atmosferica alla piroconvezione usando ERA5,
seguendo il quadro concettuale di Castellnou et al. (2022):

    0 = nessuna classe / condizioni non sufficienti
    1 = Convective plume
    2 = Overshooting PyroCu potenziale
    3 = Resilient PyroCu potenziale
    4 = PyroCb / deep PyroCu potenziale

Il programma:
  1. scarica ERA5 single levels e pressure levels tramite CDS API;
  2. ricostruisce i profili verticali in quota AGL;
  3. stima ABL con Bulk Richardson number (fallback: BLH ERA5);
  4. calcola LCL, stabilità ABL, gradiente di theta sopra ABL,
     quota del massimo wind shear e RH al top ABL;
  5. classifica ogni pixel;
  6. esporta NetCDF, GeoTIFF orari, PNG orari e tavola multi-pannello.
     I PNG possono includere i confini amministrativi regionali.

IMPORTANTE
----------
ERA5 è una rianalisi, non una previsione operativa. Il prodotto ottenuto è
quindi una diagnosi/hindcast della predisposizione atmosferica. Per una vera
previsione la stessa diagnostica va applicata a campi NWP previsionali.

La classificazione del paper è prototipale e non fornisce una tabella di
soglie universali. Le soglie numeriche qui adottate sono quindi iniziali,
"paper-informed", e devono essere validate/calibrate su casi locali.

Installazione minima:
    pip install "cdsapi>=0.7.7" xarray netCDF4 numpy pandas matplotlib rasterio

Cartopy è opzionale, ma migliora le mappe:
    conda install -c conda-forge cartopy

Configurazione CDS API (~/.cdsapirc):
    url: https://cds.climate.copernicus.eu/api
    key: <PERSONAL-ACCESS-TOKEN>

Esempio:
    python pyroconvection_era5.py \
        --start 2026-06-20 --end 2026-06-24 \
        --hours 00 03 06 09 12 15 18 21 \
        --area 45.2 8.8 41.4 13.5 \
        --outdir output_piroconvezione

Area CDS: NORD OVEST SUD EST.
"""

from __future__ import annotations

import argparse
import logging
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
import xarray as xr


# -----------------------------------------------------------------------------
# Costanti fisiche
# -----------------------------------------------------------------------------

G = 9.80665                 # m s-2
RD = 287.05                 # J kg-1 K-1
CP = 1004.0                 # J kg-1 K-1
EPSILON = 0.622             # Rd/Rv circa
KAPPA = RD / CP
P0 = 100000.0               # Pa
RIB_CRITICAL = 0.33
NODATA_CLASS = 255


# -----------------------------------------------------------------------------
# Soglie operative iniziali, derivate dal quadro concettuale del paper
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class Thresholds:
    # Il paper distingue ABL neutra/instabile da stabile intorno a 1e-3 K m-1.
    theta_grad_stable: float = 1.0e-3

    # Overshooting può avvenire anche in ABL leggermente stabile.
    theta_grad_overshoot_max: float = 3.0e-3

    # Nel paper 3.9e-3 favorisce deepening e 5.1e-3 lo inibisce.
    # Non sono soglie universali: 4.2e-3 è una soglia operativa iniziale.
    gamma_theta_weak_cap: float = 4.2e-3
    gamma_theta_strong_cap: float = 4.8e-3

    # Rapporti geometrici e prossimità normalizzate rispetto ad ABL.
    lcl_ratio_overshoot_max: float = 1.60
    lcl_ratio_deep_max: float = 1.10
    lcl_ratio_resilient_max: float = 1.00
    shear_distance_deep: float = 0.30

    # Umidità relativa media attorno al top ABL.
    rh_top_moist: float = 80.0

    # Limiti di qualità/ricerca.
    min_abl_m: float = 150.0
    max_abl_m: float = 4500.0
    max_profile_m: float = 8000.0
    min_land_fraction: float = 0.50


DEFAULT_THRESHOLDS = Thresholds()


@dataclass(frozen=True)
class Domain:
    north: float
    west: float
    south: float
    east: float

    def as_cds_area(self) -> list[float]:
        return [self.north, self.west, self.south, self.east]


@dataclass(frozen=True)
class RunConfig:
    start: str
    end: str
    hours: tuple[str, ...]
    area: Domain
    outdir: Path
    pressure_levels: tuple[int, ...]
    download: bool = True
    overwrite_downloads: bool = False
    mask_sea: bool = True
    make_individual_png: bool = True
    make_panel_png: bool = True
    panel_cols: int = 8
    use_cartopy: bool = True
    draw_regional_boundaries: bool = False
    regions_shapefile: Path | None = None


@dataclass
class ProfileMetrics:
    abl_m: float = np.nan
    lcl_m: float = np.nan
    lcl_ratio: float = np.nan
    z_shear_m: float = np.nan
    shear_distance_ratio: float = np.nan
    theta_grad_abl: float = np.nan
    gamma_theta_ft: float = np.nan
    rh_top_abl: float = np.nan
    score: float = np.nan
    pyro_class: int = 0


CLASS_NAMES = {
    0: "Nessuna classe",
    1: "Convective plumes",
    2: "Overshooting PyroCu",
    3: "Resilient PyroCu",
    4: "PyroCb",
}

# Colori coerenti con la legenda allegata.
CLASS_COLORS = {
    0: "#FFFFFF",
    1: "#99E699",
    2: "#FFF32B",
    3: "#F7C83E",
    4: "#8B0000",
}


# -----------------------------------------------------------------------------
# Logging
# -----------------------------------------------------------------------------

def setup_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


# -----------------------------------------------------------------------------
# CDS API: costruzione richiesta e download
# -----------------------------------------------------------------------------

def normalize_hours(hours: Iterable[str]) -> tuple[str, ...]:
    normalized: list[str] = []
    for raw in hours:
        h = str(raw).strip()
        if ":" not in h:
            h = f"{int(h):02d}:00"
        else:
            hh, mm = h.split(":", 1)
            h = f"{int(hh):02d}:{int(mm):02d}"
        if h not in normalized:
            normalized.append(h)
    return tuple(normalized)


def split_dates(start: str, end: str) -> tuple[list[str], list[str], list[str]]:
    dates = pd.date_range(start=start, end=end, freq="D")
    if len(dates) == 0:
        raise ValueError("Intervallo temporale vuoto.")
    years = sorted({f"{d.year:04d}" for d in dates})
    months = sorted({f"{d.month:02d}" for d in dates})
    days = sorted({f"{d.day:02d}" for d in dates})
    return years, months, days


def request_is_single_calendar_rectangle(start: str, end: str) -> bool:
    """
    I campi year/month/day del CDS formano un prodotto cartesiano. È sicuro
    usarli in una sola richiesta solo quando l'intervallo ricade nello stesso
    mese; altrimenti si scarica mese per mese.
    """
    a = pd.Timestamp(start)
    b = pd.Timestamp(end)
    return a.year == b.year and a.month == b.month


def month_chunks(start: str, end: str) -> list[tuple[pd.Timestamp, pd.Timestamp]]:
    start_ts = pd.Timestamp(start).normalize()
    end_ts = pd.Timestamp(end).normalize()
    chunks: list[tuple[pd.Timestamp, pd.Timestamp]] = []
    current = start_ts
    while current <= end_ts:
        month_end = (current + pd.offsets.MonthEnd(0)).normalize()
        chunk_end = min(month_end, end_ts)
        chunks.append((current, chunk_end))
        current = chunk_end + pd.Timedelta(days=1)
    return chunks


def cds_request_common(
    chunk_start: pd.Timestamp,
    chunk_end: pd.Timestamp,
    hours: Sequence[str],
    area: Domain,
) -> dict:
    dates = pd.date_range(chunk_start, chunk_end, freq="D")
    return {
        "product_type": ["reanalysis"],
        "year": sorted({f"{d.year:04d}" for d in dates}),
        "month": sorted({f"{d.month:02d}" for d in dates}),
        "day": [f"{d.day:02d}" for d in dates],
        "time": list(hours),
        "area": area.as_cds_area(),
        "data_format": "netcdf",
        "download_format": "unarchived",
    }


def download_one_request(dataset: str, request: dict, target: Path) -> None:
    try:
        import cdsapi
    except ImportError as exc:
        raise RuntimeError(
            'Manca cdsapi. Installare con: pip install "cdsapi>=0.7.7"'
        ) from exc

    target.parent.mkdir(parents=True, exist_ok=True)
    logging.info("CDS download %s -> %s", dataset, target)
    client = cdsapi.Client()
    client.retrieve(dataset, request, str(target))

    if not target.exists() or target.stat().st_size == 0:
        raise RuntimeError(f"Download non riuscito o file vuoto: {target}")


def download_era5(cfg: RunConfig) -> tuple[list[Path], list[Path]]:
    """Scarica ERA5 mese per mese per evitare prodotti cartesiani sulle date."""
    single_files: list[Path] = []
    pressure_files: list[Path] = []
    raw_dir = cfg.outdir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    for chunk_start, chunk_end in month_chunks(cfg.start, cfg.end):
        tag = f"{chunk_start:%Y%m%d}_{chunk_end:%Y%m%d}"
        sl_path = raw_dir / f"era5_single_{tag}.nc"
        pl_path = raw_dir / f"era5_pressure_{tag}.nc"

        common = cds_request_common(chunk_start, chunk_end, cfg.hours, cfg.area)

        sl_request = {
            **common,
            "variable": [
                "2m_temperature",
                "2m_dewpoint_temperature",
                "10m_u_component_of_wind",
                "10m_v_component_of_wind",
                "surface_pressure",
                "boundary_layer_height",
                "geopotential",
                "land_sea_mask",
            ],
        }

        pl_request = {
            **common,
            "variable": [
                "geopotential",
                "temperature",
                "specific_humidity",
                "u_component_of_wind",
                "v_component_of_wind",
            ],
            "pressure_level": [str(v) for v in cfg.pressure_levels],
        }

        if cfg.overwrite_downloads or not sl_path.exists():
            download_one_request("reanalysis-era5-single-levels", sl_request, sl_path)
        else:
            logging.info("File già presente, salto download: %s", sl_path)

        if cfg.overwrite_downloads or not pl_path.exists():
            download_one_request("reanalysis-era5-pressure-levels", pl_request, pl_path)
        else:
            logging.info("File già presente, salto download: %s", pl_path)

        single_files.append(sl_path)
        pressure_files.append(pl_path)

    return single_files, pressure_files


# -----------------------------------------------------------------------------
# Apertura e normalizzazione dataset ERA5
# -----------------------------------------------------------------------------

def _rename_common(ds: xr.Dataset) -> xr.Dataset:
    rename: dict[str, str] = {}
    candidates = {
        "valid_time": "time",
        "pressure_level": "level",
        "isobaricInhPa": "level",
        "lat": "latitude",
        "lon": "longitude",
    }
    for old, new in candidates.items():
        if old in ds.dims or old in ds.coords:
            if new not in ds.dims and new not in ds.coords:
                rename[old] = new
    if rename:
        ds = ds.rename(rename)
    return ds


def _drop_or_combine_expver(ds: xr.Dataset) -> xr.Dataset:
    if "expver" not in ds.dims:
        return ds
    # ERA5 può avere due versioni complementari. combine_first mantiene i dati validi.
    pieces = [ds.isel(expver=i, drop=True) for i in range(ds.sizes["expver"])]
    out = pieces[0]
    for p in pieces[1:]:
        out = out.combine_first(p)
    return out


def open_many_netcdf(paths: Sequence[Path]) -> xr.Dataset:
    missing = [str(p) for p in paths if not p.exists()]
    if missing:
        raise FileNotFoundError("File mancanti:\n" + "\n".join(missing))

    datasets: list[xr.Dataset] = []
    for p in paths:
        logging.info("Apro %s", p)
        ds = xr.open_dataset(p, decode_times=True)
        ds = _rename_common(ds)
        ds = _drop_or_combine_expver(ds)
        datasets.append(ds)

    if len(datasets) == 1:
        out = datasets[0]
    else:
        out = xr.concat(datasets, dim="time", data_vars="minimal", coords="minimal", compat="override")

    if "time" not in out.coords:
        raise ValueError("Coordinata temporale non trovata nei NetCDF ERA5.")
    out = out.sortby("time")
    _, unique_idx = np.unique(out["time"].values, return_index=True)
    out = out.isel(time=np.sort(unique_idx))
    return out


def select_requested_time(ds: xr.Dataset, start: str, end: str, hours: Sequence[str]) -> xr.Dataset:
    start_ts = pd.Timestamp(start)
    end_ts = pd.Timestamp(end) + pd.Timedelta(days=1) - pd.Timedelta(microseconds=1)
    ds = ds.sel(time=slice(start_ts, end_ts))
    wanted_hours = {int(h.split(":")[0]) for h in hours}
    mask = xr.DataArray(
        np.array([pd.Timestamp(t).hour in wanted_hours for t in ds.time.values]),
        dims=("time",),
        coords={"time": ds.time},
    )
    return ds.sel(time=mask)


def find_var(ds: xr.Dataset, candidates: Sequence[str]) -> xr.DataArray:
    for name in candidates:
        if name in ds.data_vars:
            return ds[name]
    raise KeyError(
        f"Nessuna variabile tra {list(candidates)}. Presenti: {list(ds.data_vars)}"
    )


def ensure_lat_lon_order(ds: xr.Dataset) -> xr.Dataset:
    if ds.latitude.values[0] < ds.latitude.values[-1]:
        ds = ds.sortby("latitude", ascending=False)
    if ds.longitude.values[0] > ds.longitude.values[-1]:
        ds = ds.sortby("longitude", ascending=True)
    return ds


def align_datasets(ds_sl: xr.Dataset, ds_pl: xr.Dataset) -> tuple[xr.Dataset, xr.Dataset]:
    ds_sl = ensure_lat_lon_order(ds_sl)
    ds_pl = ensure_lat_lon_order(ds_pl)

    # Pressione ordinata dal basso verso l'alto: 1000 -> 100 hPa.
    if ds_pl.level.values[0] < ds_pl.level.values[-1]:
        ds_pl = ds_pl.sortby("level", ascending=False)

    common_time = np.intersect1d(ds_sl.time.values, ds_pl.time.values)
    if len(common_time) == 0:
        raise ValueError("Nessun timestamp comune fra single levels e pressure levels.")

    ds_sl = ds_sl.sel(time=common_time)
    ds_pl = ds_pl.sel(time=common_time)

    # In caso di minime differenze di griglia, interpola i pressure levels sulla griglia single.
    same_lat = np.array_equal(ds_sl.latitude.values, ds_pl.latitude.values)
    same_lon = np.array_equal(ds_sl.longitude.values, ds_pl.longitude.values)
    if not (same_lat and same_lon):
        logging.warning("Griglie non coincidenti: interpolo pressure levels su single levels.")
        ds_pl = ds_pl.interp(latitude=ds_sl.latitude, longitude=ds_sl.longitude)

    return ds_sl, ds_pl


# -----------------------------------------------------------------------------
# Termodinamica
# -----------------------------------------------------------------------------

def saturation_vapor_pressure_pa(t_k: np.ndarray | float) -> np.ndarray:
    """Magnus-Tetens su acqua, sufficiente nel range operativo."""
    t_c = np.asarray(t_k, dtype=float) - 273.15
    return 611.2 * np.exp(17.67 * t_c / (t_c + 243.5))


def specific_humidity_from_dewpoint(td_k: float, p_pa: float) -> float:
    e = float(saturation_vapor_pressure_pa(td_k))
    denominator = p_pa - (1.0 - EPSILON) * e
    if denominator <= 0:
        return np.nan
    return EPSILON * e / denominator


def vapor_pressure_from_specific_humidity(q: np.ndarray, p_pa: np.ndarray) -> np.ndarray:
    return q * p_pa / (EPSILON + (1.0 - EPSILON) * q)


def relative_humidity_from_q_t_p(q: np.ndarray, t_k: np.ndarray, p_pa: np.ndarray) -> np.ndarray:
    e = vapor_pressure_from_specific_humidity(q, p_pa)
    es = saturation_vapor_pressure_pa(t_k)
    return np.clip(100.0 * e / es, 0.0, 100.0)


def potential_temperature(t_k: np.ndarray | float, p_pa: np.ndarray | float) -> np.ndarray:
    return np.asarray(t_k, dtype=float) * (P0 / np.asarray(p_pa, dtype=float)) ** KAPPA


def virtual_potential_temperature(theta_k: np.ndarray, q: np.ndarray) -> np.ndarray:
    return theta_k * (1.0 + 0.61 * q)


def lcl_height_bolton_m(t_k: float, td_k: float, p_pa: float) -> float:
    """
    LCL con formula di Bolton per Tlcl e conversione ipsometrica della pressione
    LCL in quota AGL. Restituisce metri sopra il punto di partenza.
    """
    if not np.isfinite([t_k, td_k, p_pa]).all() or p_pa <= 0:
        return np.nan
    td_k = min(td_k, t_k)
    if t_k - td_k < 0.02:
        return 0.0

    tlcl = 1.0 / (1.0 / (td_k - 56.0) + math.log(t_k / td_k) / 800.0) + 56.0
    p_lcl = p_pa * (tlcl / t_k) ** (CP / RD)

    q_sfc = specific_humidity_from_dewpoint(td_k, p_pa)
    tv_mean = 0.5 * (t_k + tlcl) * (1.0 + 0.61 * q_sfc)
    z_lcl = RD * tv_mean / G * math.log(p_pa / p_lcl)
    return max(0.0, float(z_lcl))


# -----------------------------------------------------------------------------
# Utilità profilo verticale
# -----------------------------------------------------------------------------

def clean_profile(z: np.ndarray, *values: np.ndarray) -> tuple[np.ndarray, ...]:
    arrays = [np.asarray(z, dtype=float)] + [np.asarray(v, dtype=float) for v in values]
    mask = np.isfinite(arrays[0])
    for arr in arrays[1:]:
        mask &= np.isfinite(arr)
    arrays = [arr[mask] for arr in arrays]
    if len(arrays[0]) < 2:
        return tuple(arrays)

    order = np.argsort(arrays[0])
    arrays = [arr[order] for arr in arrays]
    unique = np.concatenate(([True], np.diff(arrays[0]) > 0.5))
    arrays = [arr[unique] for arr in arrays]
    return tuple(arrays)


def interp_to_grid(z: np.ndarray, x: np.ndarray, grid: np.ndarray) -> np.ndarray:
    zc, xc = clean_profile(z, x)
    if len(zc) < 2:
        return np.full_like(grid, np.nan, dtype=float)
    out = np.interp(grid, zc, xc)
    out[(grid < zc.min()) | (grid > zc.max())] = np.nan
    return out


def linear_gradient(z: np.ndarray, x: np.ndarray, zmin: float, zmax: float) -> float:
    mask = np.isfinite(z) & np.isfinite(x) & (z >= zmin) & (z <= zmax)
    if mask.sum() < 3 or zmax <= zmin:
        return np.nan
    coeff = np.polyfit(z[mask], x[mask], deg=1)
    return float(coeff[0])


def mean_in_layer(z: np.ndarray, x: np.ndarray, zmin: float, zmax: float) -> float:
    mask = np.isfinite(z) & np.isfinite(x) & (z >= zmin) & (z <= zmax)
    if mask.sum() == 0:
        return np.nan
    return float(np.nanmean(x[mask]))


def first_threshold_crossing(z: np.ndarray, y: np.ndarray, threshold: float, zmin: float) -> float:
    mask = np.isfinite(z) & np.isfinite(y) & (z >= zmin)
    zz = z[mask]
    yy = y[mask]
    if len(zz) < 2:
        return np.nan

    idx = np.where(yy >= threshold)[0]
    if len(idx) == 0:
        return np.nan
    i = int(idx[0])
    if i == 0:
        return float(zz[0])

    z0, z1 = zz[i - 1], zz[i]
    y0, y1 = yy[i - 1], yy[i]
    if not np.isfinite(y0) or y1 == y0:
        return float(z1)
    frac = (threshold - y0) / (y1 - y0)
    return float(z0 + np.clip(frac, 0.0, 1.0) * (z1 - z0))


def estimate_abl_rib(
    z_agl: np.ndarray,
    theta_v: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    theta_v_surface: float,
    u10: float,
    v10: float,
    fallback_blh: float,
    thresholds: Thresholds,
) -> float:
    """
    Stima ABL dal primo superamento Rib >= 0.33 a partire da 200 m AGL.
    Per il profilo ambientale ERA5 si usa il riferimento di superficie; la
    scelta di 400 m nel paper riguarda soprattutto le sonde dentro il plume,
    per escludere l'indraft del fuoco.
    """
    if not np.isfinite([theta_v_surface, u10, v10]).all():
        return float(fallback_blh) if np.isfinite(fallback_blh) else np.nan

    wind2 = (u - u10) ** 2 + (v - v10) ** 2
    wind2 = np.maximum(wind2, 0.5)
    rib = (G / theta_v_surface) * (theta_v - theta_v_surface) * z_agl / wind2
    abl = first_threshold_crossing(z_agl, rib, RIB_CRITICAL, zmin=200.0)

    if not np.isfinite(abl) or not (thresholds.min_abl_m <= abl <= thresholds.max_abl_m):
        abl = float(fallback_blh) if np.isfinite(fallback_blh) else np.nan

    if np.isfinite(abl):
        abl = float(np.clip(abl, thresholds.min_abl_m, thresholds.max_abl_m))
    return abl


def maximum_shear_height(
    z_grid: np.ndarray,
    u_grid: np.ndarray,
    v_grid: np.ndarray,
    zmin: float = 200.0,
    zmax: float = 6000.0,
    half_window_m: float = 100.0,
) -> float:
    """Quota del massimo shear vettoriale calcolato su finestra di 200 m."""
    dz = float(np.nanmedian(np.diff(z_grid)))
    half_steps = max(1, int(round(half_window_m / dz)))
    shear = np.full_like(z_grid, np.nan, dtype=float)

    for i in range(half_steps, len(z_grid) - half_steps):
        if not np.isfinite(
            [u_grid[i - half_steps], u_grid[i + half_steps],
             v_grid[i - half_steps], v_grid[i + half_steps]]
        ).all():
            continue
        span = z_grid[i + half_steps] - z_grid[i - half_steps]
        du = u_grid[i + half_steps] - u_grid[i - half_steps]
        dv = v_grid[i + half_steps] - v_grid[i - half_steps]
        shear[i] = math.hypot(du, dv) / span

    mask = np.isfinite(shear) & (z_grid >= zmin) & (z_grid <= zmax)
    if not mask.any():
        return np.nan
    local_indices = np.where(mask)[0]
    return float(z_grid[local_indices[np.nanargmax(shear[mask])]])


# -----------------------------------------------------------------------------
# Classificazione
# -----------------------------------------------------------------------------

def continuous_score(
    lcl_ratio: float,
    theta_grad_abl: float,
    gamma_theta_ft: float,
    shear_distance_ratio: float,
    rh_top_abl: float,
    th: Thresholds,
) -> float:
    if not np.isfinite(
        [lcl_ratio, theta_grad_abl, gamma_theta_ft, shear_distance_ratio, rh_top_abl]
    ).all():
        return np.nan

    # 0-100: score diagnostico, non una probabilità statistica.
    lcl_term = np.clip(1.0 - abs(lcl_ratio - 1.0) / 1.0, 0.0, 1.0)
    stability_term = np.clip(
        (th.theta_grad_overshoot_max - theta_grad_abl)
        / (th.theta_grad_overshoot_max + 1.0e-3),
        0.0,
        1.0,
    )
    cap_term = np.clip(
        (th.gamma_theta_strong_cap - gamma_theta_ft)
        / (th.gamma_theta_strong_cap - 2.5e-3),
        0.0,
        1.0,
    )
    shear_term = np.clip(1.0 - shear_distance_ratio / 1.0, 0.0, 1.0)
    moisture_term = np.clip((rh_top_abl - 55.0) / 40.0, 0.0, 1.0)

    return float(
        np.clip(
            30.0 * lcl_term
            + 25.0 * stability_term
            + 20.0 * cap_term
            + 15.0 * shear_term
            + 10.0 * moisture_term,
            0.0,
            100.0,
        )
    )


def classify_pyroconvection(
    abl_m: float,
    lcl_m: float,
    z_shear_m: float,
    theta_grad_abl: float,
    gamma_theta_ft: float,
    rh_top_abl: float,
    th: Thresholds,
) -> tuple[int, float, float, float]:
    required = [abl_m, lcl_m, z_shear_m, theta_grad_abl, gamma_theta_ft, rh_top_abl]
    if not np.isfinite(required).all() or abl_m <= 0:
        return 0, np.nan, np.nan, np.nan

    lcl_ratio = lcl_m / abl_m
    shear_distance_ratio = min(
        abs(z_shear_m - abl_m),
        abs(z_shear_m - lcl_m),
    ) / abl_m

    score = continuous_score(
        lcl_ratio,
        theta_grad_abl,
        gamma_theta_ft,
        shear_distance_ratio,
        rh_top_abl,
        th,
    )

    stable = theta_grad_abl > th.theta_grad_stable
    neutral_or_unstable = theta_grad_abl <= th.theta_grad_stable
    slightly_stable_or_better = theta_grad_abl <= th.theta_grad_overshoot_max
    weak_cap = gamma_theta_ft <= th.gamma_theta_weak_cap
    strong_cap = gamma_theta_ft >= th.gamma_theta_strong_cap
    shear_close = shear_distance_ratio <= th.shear_distance_deep
    moist_top = rh_top_abl >= th.rh_top_moist

    # 4. Deep pyroCu / PyroCb potenziale:
    # LCL e ABL vicini, ABL neutra/instabile, cap debole e shear vicino.
    if (
        lcl_ratio <= th.lcl_ratio_deep_max
        and neutral_or_unstable
        and weak_cap
        and shear_close
        and moist_top
    ):
        return 4, score, lcl_ratio, shear_distance_ratio

    # 3. Resilient PyroCu:
    # LCL dentro ABL, condensazione persistente, ma crescita verticale limitata
    # da cap e/o shear lontano.
    if (
        lcl_ratio < th.lcl_ratio_resilient_max
        and neutral_or_unstable
        and moist_top
        and (strong_cap or not shear_close or not weak_cap)
    ):
        return 3, score, lcl_ratio, shear_distance_ratio

    # 2. Overshooting PyroCu potenziale:
    # LCL sopra ABL ma non troppo lontano; richiede un impulso di calore del fuoco,
    # non ricavabile dal solo ERA5. La classe indica quindi predisposizione.
    if (
        1.0 <= lcl_ratio <= th.lcl_ratio_overshoot_max
        and slightly_stable_or_better
    ):
        return 2, score, lcl_ratio, shear_distance_ratio

    # 1. Convective plume senza nube persistente: LCL sopra ABL.
    if lcl_ratio > 1.0:
        return 1, score, lcl_ratio, shear_distance_ratio

    # Caso residuo con LCL sotto ABL ma profilo non abbastanza coerente/umido.
    # Si mantiene una classe prudente solo se l'ambiente è almeno ben rimescolato.
    if neutral_or_unstable and score >= 35.0:
        return 2, score, lcl_ratio, shear_distance_ratio

    return 0, score, lcl_ratio, shear_distance_ratio


# -----------------------------------------------------------------------------
# Elaborazione di un profilo
# -----------------------------------------------------------------------------

def process_profile(
    pressure_hpa: np.ndarray,
    geopotential: np.ndarray,
    temperature_k: np.ndarray,
    q: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    surface_geopotential: float,
    surface_pressure_pa: float,
    t2m_k: float,
    td2m_k: float,
    u10: float,
    v10: float,
    blh_fallback_m: float,
    thresholds: Thresholds,
) -> ProfileMetrics:
    metrics = ProfileMetrics()

    scalar_values = [
        surface_geopotential,
        surface_pressure_pa,
        t2m_k,
        td2m_k,
        u10,
        v10,
    ]
    if not np.isfinite(scalar_values).all() or surface_pressure_pa <= 0:
        return metrics

    p_pa = np.asarray(pressure_hpa, dtype=float) * 100.0
    z_asl = np.asarray(geopotential, dtype=float) / G
    orog_m = surface_geopotential / G
    z_agl = z_asl - orog_m

    # Esclude livelli sotto il terreno e oltre il dominio verticale utile.
    valid = (
        np.isfinite(z_agl)
        & np.isfinite(temperature_k)
        & np.isfinite(q)
        & np.isfinite(u)
        & np.isfinite(v)
        & (z_agl >= 20.0)
        & (z_agl <= thresholds.max_profile_m)
        & (p_pa <= surface_pressure_pa + 100.0)
    )
    if valid.sum() < 4:
        return metrics

    z_agl = z_agl[valid]
    p_pa = p_pa[valid]
    temperature_k = np.asarray(temperature_k, dtype=float)[valid]
    q = np.asarray(q, dtype=float)[valid]
    u = np.asarray(u, dtype=float)[valid]
    v = np.asarray(v, dtype=float)[valid]

    z_agl, p_pa, temperature_k, q, u, v = clean_profile(
        z_agl, p_pa, temperature_k, q, u, v
    )
    if len(z_agl) < 4:
        return metrics

    theta = potential_temperature(temperature_k, p_pa)
    theta_v = virtual_potential_temperature(theta, q)

    q_sfc = specific_humidity_from_dewpoint(td2m_k, surface_pressure_pa)
    theta_sfc = float(potential_temperature(t2m_k, surface_pressure_pa))
    theta_v_sfc = theta_sfc * (1.0 + 0.61 * q_sfc)

    abl_m = estimate_abl_rib(
        z_agl=z_agl,
        theta_v=theta_v,
        u=u,
        v=v,
        theta_v_surface=theta_v_sfc,
        u10=u10,
        v10=v10,
        fallback_blh=blh_fallback_m,
        thresholds=thresholds,
    )
    if not np.isfinite(abl_m):
        return metrics

    lcl_m = lcl_height_bolton_m(t2m_k, td2m_k, surface_pressure_pa)

    # Griglia AGL regolare per gradienti e shear.
    z_top = min(thresholds.max_profile_m, float(np.nanmax(z_agl)))
    z_grid = np.arange(50.0, z_top + 25.0, 50.0)
    if len(z_grid) < 10:
        return metrics

    theta_grid = interp_to_grid(z_agl, theta, z_grid)
    t_grid = interp_to_grid(z_agl, temperature_k, z_grid)
    q_grid = interp_to_grid(z_agl, q, z_grid)
    p_grid = interp_to_grid(z_agl, p_pa, z_grid)
    u_grid = interp_to_grid(z_agl, u, z_grid)
    v_grid = interp_to_grid(z_agl, v, z_grid)
    rh_grid = relative_humidity_from_q_t_p(q_grid, t_grid, p_grid)

    z_shear_m = maximum_shear_height(z_grid, u_grid, v_grid)

    # Stabilità media nell'ABL: parte centrale per evitare superficie e inversione top.
    abl_low = max(100.0, 0.20 * abl_m)
    abl_high = max(abl_low + 150.0, 0.80 * abl_m)
    abl_high = min(abl_high, abl_m - 50.0 if abl_m > 250.0 else abl_m)
    theta_grad_abl = linear_gradient(z_grid, theta_grid, abl_low, abl_high)

    # Free-troposphere theta lapse rate sopra ABL.
    ft_low = abl_m + 200.0
    ft_high = min(abl_m + 1200.0, z_top)
    gamma_theta_ft = linear_gradient(z_grid, theta_grid, ft_low, ft_high)

    # RH media in ±150 m attorno al top ABL.
    rh_top_abl = mean_in_layer(z_grid, rh_grid, max(50.0, abl_m - 150.0), abl_m + 150.0)

    pyro_class, score, lcl_ratio, shear_distance_ratio = classify_pyroconvection(
        abl_m=abl_m,
        lcl_m=lcl_m,
        z_shear_m=z_shear_m,
        theta_grad_abl=theta_grad_abl,
        gamma_theta_ft=gamma_theta_ft,
        rh_top_abl=rh_top_abl,
        th=thresholds,
    )

    metrics.abl_m = abl_m
    metrics.lcl_m = lcl_m
    metrics.lcl_ratio = lcl_ratio
    metrics.z_shear_m = z_shear_m
    metrics.shear_distance_ratio = shear_distance_ratio
    metrics.theta_grad_abl = theta_grad_abl
    metrics.gamma_theta_ft = gamma_theta_ft
    metrics.rh_top_abl = rh_top_abl
    metrics.score = score
    metrics.pyro_class = pyro_class
    return metrics


# -----------------------------------------------------------------------------
# Elaborazione dominio
# -----------------------------------------------------------------------------

def process_domain(
    ds_sl: xr.Dataset,
    ds_pl: xr.Dataset,
    thresholds: Thresholds,
    mask_sea: bool = True,
) -> xr.Dataset:
    # Variabili con alias tipici dei NetCDF CDS.
    t2m = find_var(ds_sl, ["t2m", "2m_temperature"])
    d2m = find_var(ds_sl, ["d2m", "2m_dewpoint_temperature"])
    u10_da = find_var(ds_sl, ["u10", "10m_u_component_of_wind"])
    v10_da = find_var(ds_sl, ["v10", "10m_v_component_of_wind"])
    sp = find_var(ds_sl, ["sp", "surface_pressure"])
    blh = find_var(ds_sl, ["blh", "boundary_layer_height"])
    z_sfc = find_var(ds_sl, ["z", "geopotential"])
    lsm = find_var(ds_sl, ["lsm", "land_sea_mask"])

    z_pl = find_var(ds_pl, ["z", "geopotential"])
    t_pl = find_var(ds_pl, ["t", "temperature"])
    q_pl = find_var(ds_pl, ["q", "specific_humidity"])
    u_pl = find_var(ds_pl, ["u", "u_component_of_wind"])
    v_pl = find_var(ds_pl, ["v", "v_component_of_wind"])

    # Garantisce ordine dimensionale stabile.
    pl_order = ("time", "level", "latitude", "longitude")
    sl_order = ("time", "latitude", "longitude")
    z_pl = z_pl.transpose(*pl_order)
    t_pl = t_pl.transpose(*pl_order)
    q_pl = q_pl.transpose(*pl_order)
    u_pl = u_pl.transpose(*pl_order)
    v_pl = v_pl.transpose(*pl_order)

    def broadcast_sl(da: xr.DataArray) -> xr.DataArray:
        # Geopotenziale/LSM possono essere statici e privi di time.
        if "time" not in da.dims:
            da = da.expand_dims(time=ds_sl.time)
        da = da.broadcast_like(t2m)
        return da.transpose(*sl_order)

    t2m = broadcast_sl(t2m)
    d2m = broadcast_sl(d2m)
    u10_da = broadcast_sl(u10_da)
    v10_da = broadcast_sl(v10_da)
    sp = broadcast_sl(sp)
    blh = broadcast_sl(blh)
    z_sfc = broadcast_sl(z_sfc)
    lsm = broadcast_sl(lsm)

    times = ds_sl.time.values
    lats = ds_sl.latitude.values
    lons = ds_sl.longitude.values
    levels = ds_pl.level.values.astype(float)

    nt, ny, nx = len(times), len(lats), len(lons)
    shape = (nt, ny, nx)

    pyro_class = np.zeros(shape, dtype=np.uint8)
    score = np.full(shape, np.nan, dtype=np.float32)
    abl_m = np.full(shape, np.nan, dtype=np.float32)
    lcl_m = np.full(shape, np.nan, dtype=np.float32)
    lcl_ratio = np.full(shape, np.nan, dtype=np.float32)
    z_shear_m = np.full(shape, np.nan, dtype=np.float32)
    shear_distance_ratio = np.full(shape, np.nan, dtype=np.float32)
    theta_grad_abl = np.full(shape, np.nan, dtype=np.float32)
    gamma_theta_ft = np.full(shape, np.nan, dtype=np.float32)
    rh_top_abl = np.full(shape, np.nan, dtype=np.float32)

    z_pl_np = np.asarray(z_pl.values)
    t_pl_np = np.asarray(t_pl.values)
    q_pl_np = np.asarray(q_pl.values)
    u_pl_np = np.asarray(u_pl.values)
    v_pl_np = np.asarray(v_pl.values)

    t2m_np = np.asarray(t2m.values)
    d2m_np = np.asarray(d2m.values)
    u10_np = np.asarray(u10_da.values)
    v10_np = np.asarray(v10_da.values)
    sp_np = np.asarray(sp.values)
    blh_np = np.asarray(blh.values)
    z_sfc_np = np.asarray(z_sfc.values)
    lsm_np = np.asarray(lsm.values)

    total = nt * ny * nx
    done = 0
    for it, timestamp in enumerate(times):
        logging.info("Elaboro %s (%d/%d)", pd.Timestamp(timestamp), it + 1, nt)
        for iy in range(ny):
            for ix in range(nx):
                done += 1
                if mask_sea and lsm_np[it, iy, ix] < thresholds.min_land_fraction:
                    continue

                m = process_profile(
                    pressure_hpa=levels,
                    geopotential=z_pl_np[it, :, iy, ix],
                    temperature_k=t_pl_np[it, :, iy, ix],
                    q=q_pl_np[it, :, iy, ix],
                    u=u_pl_np[it, :, iy, ix],
                    v=v_pl_np[it, :, iy, ix],
                    surface_geopotential=float(z_sfc_np[it, iy, ix]),
                    surface_pressure_pa=float(sp_np[it, iy, ix]),
                    t2m_k=float(t2m_np[it, iy, ix]),
                    td2m_k=float(d2m_np[it, iy, ix]),
                    u10=float(u10_np[it, iy, ix]),
                    v10=float(v10_np[it, iy, ix]),
                    blh_fallback_m=float(blh_np[it, iy, ix]),
                    thresholds=thresholds,
                )

                pyro_class[it, iy, ix] = np.uint8(m.pyro_class)
                score[it, iy, ix] = m.score
                abl_m[it, iy, ix] = m.abl_m
                lcl_m[it, iy, ix] = m.lcl_m
                lcl_ratio[it, iy, ix] = m.lcl_ratio
                z_shear_m[it, iy, ix] = m.z_shear_m
                shear_distance_ratio[it, iy, ix] = m.shear_distance_ratio
                theta_grad_abl[it, iy, ix] = m.theta_grad_abl
                gamma_theta_ft[it, iy, ix] = m.gamma_theta_ft
                rh_top_abl[it, iy, ix] = m.rh_top_abl

        logging.debug("Avanzamento %.1f%%", 100.0 * done / total)

    coords = {
        "time": times,
        "latitude": lats,
        "longitude": lons,
    }
    dims = ("time", "latitude", "longitude")

    out = xr.Dataset(
        data_vars={
            "pyro_class": (dims, pyro_class),
            "pyro_score": (dims, score),
            "abl_height": (dims, abl_m),
            "lcl_height": (dims, lcl_m),
            "lcl_abl_ratio": (dims, lcl_ratio),
            "maximum_shear_height": (dims, z_shear_m),
            "shear_distance_ratio": (dims, shear_distance_ratio),
            "theta_gradient_abl": (dims, theta_grad_abl),
            "gamma_theta_free_troposphere": (dims, gamma_theta_ft),
            "relative_humidity_abl_top": (dims, rh_top_abl),
            "land_sea_mask": (dims, lsm_np.astype(np.float32)),
        },
        coords=coords,
        attrs={
            "title": "ERA5 atmospheric predisposition to pyroconvection",
            "method": "Paper-informed diagnostic based on ABL stability, LCL/ABL, shear proximity, free-tropospheric theta gradient and RH at ABL top",
            "reference": "Castellnou et al. (2022), JGR Atmospheres, doi:10.1029/2022JD036920",
            "warning": "ERA5 is a reanalysis. Classes indicate atmospheric predisposition under an intense-fire forcing, not occurrence probability.",
            "class_0": CLASS_NAMES[0],
            "class_1": CLASS_NAMES[1],
            "class_2": CLASS_NAMES[2],
            "class_3": CLASS_NAMES[3],
            "class_4": CLASS_NAMES[4],
        },
    )

    out["pyro_class"].attrs.update(
        long_name="pyroconvection atmospheric class",
        flag_values=np.array([0, 1, 2, 3, 4], dtype=np.uint8),
        flag_meanings="none convective_plume overshooting_pyrocu resilient_pyrocu pyrocb",
    )
    out["pyro_score"].attrs.update(
        long_name="continuous atmospheric predisposition score",
        units="1",
        comment="0-100 diagnostic score; not a calibrated probability",
    )
    for name in ["abl_height", "lcl_height", "maximum_shear_height"]:
        out[name].attrs["units"] = "m"
    out["theta_gradient_abl"].attrs["units"] = "K m-1"
    out["gamma_theta_free_troposphere"].attrs["units"] = "K m-1"
    out["relative_humidity_abl_top"].attrs["units"] = "%"

    return out


# -----------------------------------------------------------------------------
# Esportazione
# -----------------------------------------------------------------------------

def save_netcdf(ds: xr.Dataset, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    encoding: dict[str, dict] = {}
    for name, da in ds.data_vars.items():
        if name == "pyro_class":
            encoding[name] = {"zlib": True, "complevel": 4, "dtype": "u1", "_FillValue": NODATA_CLASS}
        else:
            encoding[name] = {"zlib": True, "complevel": 4, "dtype": "f4"}
    try:
        ds.to_netcdf(path, encoding=encoding)
    except (ValueError, TypeError) as exc:
        logging.warning("Compressione NetCDF non disponibile (%s). Salvo senza encoding.", exc)
        ds.to_netcdf(path)
    logging.info("Salvato NetCDF: %s", path)


def prepare_raster_array(da: xr.DataArray, nodata: int = NODATA_CLASS) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    arr = np.asarray(da.values)
    lat = np.asarray(da.latitude.values)
    lon = np.asarray(da.longitude.values)

    if lat[0] < lat[-1]:
        lat = lat[::-1]
        arr = arr[::-1, :]
    if lon[0] > lon[-1]:
        lon = lon[::-1]
        arr = arr[:, ::-1]

    arr = np.where(np.isfinite(arr), arr, nodata).astype(np.uint8)
    return arr, lat, lon


def save_geotiffs(ds: xr.Dataset, outdir: Path) -> list[Path]:
    try:
        import rasterio
        from rasterio.transform import from_bounds
    except ImportError as exc:
        logging.warning("rasterio non installato: GeoTIFF non prodotti (%s)", exc)
        return []

    tif_dir = outdir / "geotiff"
    tif_dir.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []

    for timestamp in ds.time.values:
        da = ds["pyro_class"].sel(time=timestamp)
        arr, lat, lon = prepare_raster_array(da)

        dx = float(np.median(np.diff(lon)))
        dy = abs(float(np.median(np.diff(lat))))
        west = float(lon.min() - dx / 2.0)
        east = float(lon.max() + dx / 2.0)
        south = float(lat.min() - dy / 2.0)
        north = float(lat.max() + dy / 2.0)
        transform = from_bounds(west, south, east, north, len(lon), len(lat))

        tag = pd.Timestamp(timestamp).strftime("%Y%m%dT%H%MZ")
        path = tif_dir / f"pyro_class_{tag}.tif"
        with rasterio.open(
            path,
            "w",
            driver="GTiff",
            height=arr.shape[0],
            width=arr.shape[1],
            count=1,
            dtype="uint8",
            crs="EPSG:4326",
            transform=transform,
            nodata=NODATA_CLASS,
            compress="DEFLATE",
            predictor=2,
        ) as dst:
            dst.write(arr, 1)
            dst.update_tags(
                title="Atmospheric predisposition to pyroconvection",
                class_0=CLASS_NAMES[0],
                class_1=CLASS_NAMES[1],
                class_2=CLASS_NAMES[2],
                class_3=CLASS_NAMES[3],
                class_4=CLASS_NAMES[4],
                valid_time_utc=str(pd.Timestamp(timestamp)),
            )
        outputs.append(path)
        logging.info("Salvato GeoTIFF: %s", path)
    return outputs


def matplotlib_objects():
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.patches import Patch

    cmap = ListedColormap([CLASS_COLORS[i] for i in range(5)])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], cmap.N)
    handles = [
        Patch(facecolor=CLASS_COLORS[i], edgecolor="0.35", label=CLASS_NAMES[i])
        for i in range(1, 5)
    ]
    return plt, cmap, norm, handles


def get_cartopy(use_cartopy: bool):
    if not use_cartopy:
        return None
    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        return ccrs, cfeature
    except Exception as exc:
        logging.warning("Cartopy non disponibile: uso mappa lat/lon semplice (%s)", exc)
        return None


def add_regional_boundaries(
    ax,
    cartopy,
    *,
    enabled: bool = False,
    shapefile: Path | None = None,
    linewidth: float = 0.55,
) -> None:
    """Disegna i confini amministrativi regionali sulle mappe Cartopy.

    Se ``shapefile`` è valorizzato usa il file locale indicato (formato ESRI
    Shapefile, geometrie in EPSG:4326 o con coordinate geografiche). In caso
    contrario usa il layer Natural Earth 1° livello amministrativo
    ``admin_1_states_provinces_lines`` alla scala 1:10m.

    Il layer Natural Earth viene scaricato e memorizzato nella cache Cartopy
    al primo utilizzo.
    """
    if not enabled and shapefile is None:
        return
    if cartopy is None:
        logging.warning(
            "Confini regionali richiesti ma Cartopy non è disponibile; "
            "l'opzione viene ignorata."
        )
        return

    ccrs, cfeature = cartopy
    try:
        if shapefile is not None:
            from cartopy.io import shapereader

            shp_path = Path(shapefile)
            if not shp_path.exists():
                raise FileNotFoundError(f"Shapefile dei confini non trovato: {shp_path}")
            geometries = list(shapereader.Reader(str(shp_path)).geometries())
            ax.add_geometries(
                geometries,
                crs=ccrs.PlateCarree(),
                facecolor="none",
                edgecolor="0.15",
                linewidth=linewidth,
                zorder=8,
            )
        else:
            regions = cfeature.NaturalEarthFeature(
                category="cultural",
                name="admin_1_states_provinces_lines",
                scale="10m",
                facecolor="none",
            )
            ax.add_feature(
                regions,
                edgecolor="0.15",
                linewidth=linewidth,
                zorder=8,
            )
    except Exception as exc:
        logging.warning("Impossibile disegnare i confini regionali: %s", exc)


def plot_one_map(
    da: xr.DataArray,
    path: Path,
    use_cartopy: bool = True,
    draw_regional_boundaries: bool = False,
    regions_shapefile: Path | None = None,
    title_prefix: str = "Tipus de piroconvecció - ERA5",
) -> None:
    plt, cmap, norm, handles = matplotlib_objects()
    cartopy = get_cartopy(use_cartopy)

    timestamp = pd.Timestamp(da.time.values)
    lon = da.longitude.values
    lat = da.latitude.values
    values = da.values

    if cartopy:
        ccrs, cfeature = cartopy
        fig = plt.figure(figsize=(9.2, 7.0))
        ax = plt.axes(projection=ccrs.PlateCarree())
        mesh = ax.pcolormesh(lon, lat, values, cmap=cmap, norm=norm, shading="nearest", transform=ccrs.PlateCarree())
        ax.coastlines(resolution="10m", linewidth=0.8)
        ax.add_feature(cfeature.BORDERS.with_scale("10m"), linewidth=0.45)
        add_regional_boundaries(
            ax,
            cartopy,
            enabled=draw_regional_boundaries,
            shapefile=regions_shapefile,
            linewidth=0.65,
        )
        ax.add_feature(cfeature.LAKES.with_scale("10m"), facecolor="#B7DDE8", edgecolor="0.5", linewidth=0.3)
        ax.add_feature(cfeature.OCEAN.with_scale("10m"), facecolor="#B7DDE8", zorder=-2)
        ax.set_extent([float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max())])
        gl = ax.gridlines(draw_labels=True, linewidth=0.25, alpha=0.45)
        gl.top_labels = False
        gl.right_labels = False
    else:
        fig, ax = plt.subplots(figsize=(9.2, 7.0))
        mesh = ax.pcolormesh(lon, lat, values, cmap=cmap, norm=norm, shading="nearest")
        ax.set_xlabel("Longitudine")
        ax.set_ylabel("Latitudine")
        ax.grid(True, linewidth=0.25, alpha=0.45)

    ax.set_title(f"{title_prefix}\n{timestamp:%Y-%m-%d %H:%M UTC}", fontsize=13)
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.01))
    fig.text(
        0.5,
        0.055,
        "Predisposizione atmosferica sotto forcing di incendio intenso; ERA5 = rianalisi",
        ha="center",
        fontsize=8,
    )
    fig.subplots_adjust(bottom=0.14, top=0.88, left=0.08, right=0.96)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    logging.info("Salvato PNG: %s", path)


def save_individual_pngs(
    ds: xr.Dataset,
    outdir: Path,
    use_cartopy: bool,
    draw_regional_boundaries: bool = False,
    regions_shapefile: Path | None = None,
) -> list[Path]:
    png_dir = outdir / "png"
    png_dir.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    for timestamp in ds.time.values:
        tag = pd.Timestamp(timestamp).strftime("%Y%m%dT%H%MZ")
        path = png_dir / f"pyro_class_{tag}.png"
        plot_one_map(
            ds["pyro_class"].sel(time=timestamp),
            path,
            use_cartopy=use_cartopy,
            draw_regional_boundaries=draw_regional_boundaries,
            regions_shapefile=regions_shapefile,
        )
        outputs.append(path)
    return outputs


def save_panel_png(
    ds: xr.Dataset,
    path: Path,
    cols: int = 8,
    use_cartopy: bool = True,
    draw_regional_boundaries: bool = False,
    regions_shapefile: Path | None = None,
) -> Path:
    plt, cmap, norm, handles = matplotlib_objects()
    cartopy = get_cartopy(use_cartopy)

    n = ds.sizes["time"]
    cols = max(1, min(cols, n))
    rows = int(math.ceil(n / cols))
    width = max(12.0, cols * 2.15)
    height = max(4.8, rows * 2.05 + 1.0)

    subplot_kw = {"projection": cartopy[0].PlateCarree()} if cartopy else {}
    fig, axes = plt.subplots(rows, cols, figsize=(width, height), squeeze=False, subplot_kw=subplot_kw)

    lon = ds.longitude.values
    lat = ds.latitude.values
    for idx, ax in enumerate(axes.flat):
        if idx >= n:
            ax.set_visible(False)
            continue

        timestamp = ds.time.values[idx]
        values = ds["pyro_class"].isel(time=idx).values
        if cartopy:
            ccrs, cfeature = cartopy
            ax.pcolormesh(lon, lat, values, cmap=cmap, norm=norm, shading="nearest", transform=ccrs.PlateCarree())
            ax.coastlines(resolution="10m", linewidth=0.35)
            ax.add_feature(cfeature.BORDERS.with_scale("10m"), linewidth=0.2)
            add_regional_boundaries(
                ax,
                cartopy,
                enabled=draw_regional_boundaries,
                shapefile=regions_shapefile,
                linewidth=0.35,
            )
            ax.add_feature(cfeature.OCEAN.with_scale("10m"), facecolor="#B7DDE8", zorder=-2)
            ax.set_extent([float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max())])
        else:
            ax.pcolormesh(lon, lat, values, cmap=cmap, norm=norm, shading="nearest")
            ax.set_xticks([])
            ax.set_yticks([])

        ax.set_title(pd.Timestamp(timestamp).strftime("%Y-%m-%d %HZ"), fontsize=7.2, pad=2)

    fig.suptitle("Tipus de piroconvecció - ERA5 diagnostic", fontsize=13, y=0.985)
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.015), fontsize=8)
    fig.text(
        0.5,
        0.055,
        "Classificazione atmosferica paper-informed; le classi 2-4 presuppongono un incendio sufficientemente intenso",
        ha="center",
        fontsize=7.5,
    )
    fig.subplots_adjust(left=0.015, right=0.995, top=0.945, bottom=0.11, wspace=0.04, hspace=0.16)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    logging.info("Salvata tavola multi-pannello: %s", path)
    return path


def save_class_csv(ds: xr.Dataset, path: Path) -> None:
    rows = []
    for value in range(5):
        rows.append(
            {
                "value": value,
                "label": CLASS_NAMES[value],
                "hex_color": CLASS_COLORS[value],
            }
        )
    pd.DataFrame(rows).to_csv(path, index=False)
    logging.info("Salvata tabella classi: %s", path)


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def parse_pressure_levels(text_values: Sequence[str]) -> tuple[int, ...]:
    levels: list[int] = []
    for item in text_values:
        for token in item.replace(",", " ").split():
            levels.append(int(token))
    unique = sorted(set(levels), reverse=True)
    if not unique:
        raise ValueError("Lista pressure levels vuota.")
    return tuple(unique)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Scarica ERA5 e produce mappe di predisposizione atmosferica alla piroconvezione.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--start", required=True, help="Data iniziale YYYY-MM-DD")
    parser.add_argument("--end", required=True, help="Data finale YYYY-MM-DD")
    parser.add_argument(
        "--hours",
        nargs="+",
        default=["00", "03", "06", "09", "12", "15", "18", "21"],
        help="Ore UTC, ad esempio: 00 03 06 09 12 15 18 21",
    )
    parser.add_argument(
        "--area",
        nargs=4,
        type=float,
        metavar=("NORTH", "WEST", "SOUTH", "EAST"),
        default=[45.2, 8.8, 41.4, 13.5],
        help="Bounding box CDS: nord ovest sud est",
    )
    parser.add_argument("--outdir", type=Path, default=Path("output_piroconvezione"))
    parser.add_argument(
        "--pressure-levels",
        nargs="+",
        default=[
            "1000", "975", "950", "925", "900", "875", "850", "825", "800",
            "775", "750", "700", "650", "600", "550", "500", "450", "400",
            "350", "300", "250", "200",
        ],
        help="Livelli di pressione hPa, separati da spazi o virgole",
    )
    parser.add_argument(
        "--no-download",
        action="store_true",
        help="Non scarica: usa i file ERA5 già presenti in OUTDIR/raw",
    )
    parser.add_argument(
        "--single-files",
        nargs="*",
        type=Path,
        help="File NetCDF single levels già scaricati; implica --no-download",
    )
    parser.add_argument(
        "--pressure-files",
        nargs="*",
        type=Path,
        help="File NetCDF pressure levels già scaricati; implica --no-download",
    )
    parser.add_argument("--overwrite-downloads", action="store_true")
    parser.add_argument("--keep-sea", action="store_true", help="Non maschera i punti marini")
    parser.add_argument("--no-individual-png", action="store_true")
    parser.add_argument("--no-panel-png", action="store_true")
    parser.add_argument("--panel-cols", type=int, default=8)
    parser.add_argument("--no-cartopy", action="store_true")
    parser.add_argument(
        "--regional-boundaries",
        action="store_true",
        help=(
            "Disegna i confini amministrativi regionali (Natural Earth admin-1) "
            "nei PNG; richiede Cartopy"
        ),
    )
    parser.add_argument(
        "--regions-shapefile",
        type=Path,
        help=(
            "Shapefile .shp alternativo per i confini regionali. "
            "Se specificato abilita automaticamente i confini"
        ),
    )
    parser.add_argument("--verbose", action="store_true")
    return parser


def discover_existing_files(outdir: Path) -> tuple[list[Path], list[Path]]:
    raw_dir = outdir / "raw"
    single = sorted(raw_dir.glob("era5_single_*.nc"))
    pressure = sorted(raw_dir.glob("era5_pressure_*.nc"))
    if not single or not pressure:
        raise FileNotFoundError(
            f"Non trovo file ERA5 in {raw_dir}. Specificare --single-files e --pressure-files oppure rimuovere --no-download."
        )
    return single, pressure


def validate_config(cfg: RunConfig) -> None:
    start = pd.Timestamp(cfg.start)
    end = pd.Timestamp(cfg.end)
    if end < start:
        raise ValueError("--end precede --start")
    if not (-90 <= cfg.area.south < cfg.area.north <= 90):
        raise ValueError("Latitudini area non valide.")
    if not (-180 <= cfg.area.west < cfg.area.east <= 360):
        raise ValueError("Longitudini area non valide.")
    if cfg.panel_cols < 1:
        raise ValueError("--panel-cols deve essere >= 1")
    if cfg.regions_shapefile is not None and not cfg.regions_shapefile.exists():
        raise FileNotFoundError(
            f"Shapefile dei confini regionali non trovato: {cfg.regions_shapefile}"
        )


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    setup_logging(args.verbose)

    hours = normalize_hours(args.hours)
    pressure_levels = parse_pressure_levels(args.pressure_levels)
    area = Domain(*args.area)
    cfg = RunConfig(
        start=args.start,
        end=args.end,
        hours=hours,
        area=area,
        outdir=args.outdir,
        pressure_levels=pressure_levels,
        download=not args.no_download and not args.single_files and not args.pressure_files,
        overwrite_downloads=args.overwrite_downloads,
        mask_sea=not args.keep_sea,
        make_individual_png=not args.no_individual_png,
        make_panel_png=not args.no_panel_png,
        panel_cols=args.panel_cols,
        use_cartopy=not args.no_cartopy,
        draw_regional_boundaries=(
            args.regional_boundaries or args.regions_shapefile is not None
        ),
        regions_shapefile=args.regions_shapefile,
    )
    validate_config(cfg)
    cfg.outdir.mkdir(parents=True, exist_ok=True)

    if args.single_files or args.pressure_files:
        if not args.single_files or not args.pressure_files:
            parser.error("Specificare insieme --single-files e --pressure-files")
        single_files = list(args.single_files)
        pressure_files = list(args.pressure_files)
    elif cfg.download:
        single_files, pressure_files = download_era5(cfg)
    else:
        single_files, pressure_files = discover_existing_files(cfg.outdir)

    ds_sl = open_many_netcdf(single_files)
    ds_pl = open_many_netcdf(pressure_files)
    ds_sl = select_requested_time(ds_sl, cfg.start, cfg.end, cfg.hours)
    ds_pl = select_requested_time(ds_pl, cfg.start, cfg.end, cfg.hours)
    ds_sl, ds_pl = align_datasets(ds_sl, ds_pl)

    logging.info(
        "Dominio: %d tempi, %d latitudini, %d longitudini, %d livelli",
        ds_sl.sizes["time"],
        ds_sl.sizes["latitude"],
        ds_sl.sizes["longitude"],
        ds_pl.sizes["level"],
    )

    result = process_domain(
        ds_sl=ds_sl,
        ds_pl=ds_pl,
        thresholds=DEFAULT_THRESHOLDS,
        mask_sea=cfg.mask_sea,
    )

    nc_path = cfg.outdir / "pyroconvection_era5.nc"
    save_netcdf(result, nc_path)
    save_geotiffs(result, cfg.outdir)
    save_class_csv(result, cfg.outdir / "classi_piroconvezione.csv")

    if cfg.make_individual_png:
        save_individual_pngs(
            result,
            cfg.outdir,
            use_cartopy=cfg.use_cartopy,
            draw_regional_boundaries=cfg.draw_regional_boundaries,
            regions_shapefile=cfg.regions_shapefile,
        )
    if cfg.make_panel_png:
        save_panel_png(
            result,
            cfg.outdir / "pyroconvection_panel.png",
            cols=cfg.panel_cols,
            use_cartopy=cfg.use_cartopy,
            draw_regional_boundaries=cfg.draw_regional_boundaries,
            regions_shapefile=cfg.regions_shapefile,
        )

    logging.info("Elaborazione completata in %s", cfg.outdir.resolve())
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        logging.error("Interrotto dall'utente.")
        raise SystemExit(130)
    except Exception as exc:
        logging.exception("Errore: %s", exc)
        raise SystemExit(1)
