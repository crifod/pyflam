#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pyroconvection_wrf15_v1_7_4.py
===========================

Mappa la predisposizione atmosferica alla piroconvezione usando il WRF operativo ItaliaMeteo-CIMA-LaMMA a 1.5 km,
seguendo il quadro concettuale di Castellnou et al. (2022):

    0 = Surface plume; 255 = non classificabile / NoData
    1 = Convective plume
    2 = Overshooting PyroCu potenziale
    3 = Resilient PyroCu potenziale
    4 = PyroCb / deep PyroCu potenziale

Il programma:
  1. legge file WRF NetCDF/GRIB scaricati da MeteoHub (surface e pressure levels),
     applicando bbox e stride ai singoli file prima del merge;
  2. ricostruisce i profili verticali in quota AGL;
  3. stima ABL con Bulk Richardson number (fallback: PBLH del WRF); quando il primo livello
     disponibile supera già la soglia e PBLH non è valida, usa soltanto un vincolo
     superiore conservativo per riconoscere una ABL certamente inferiore a 600 m;
  4. valuta la risoluzione verticale reale dei livelli originali e sceglie
     automaticamente tra ladder completa, ladder senza shear e proxy di mixing;
  5. calcola LCL, stabilità ABL, gradiente di theta sopra ABL,
     quota del massimo wind shear e RH al top ABL;
  6. classifica ogni pixel con un indice di qualità e una motivazione esplicita;
  7. esporta NetCDF, GeoTIFF orari, PNG orari e tavola multi-pannello.
     I PNG includono per impostazione predefinita i confini provinciali italiani e possono includere anche i confini regionali.

IMPORTANTE
----------
Il dataset WRF è una previsione deterministica operativa a 1.5 km, con corse 00 e 12 UTC e orizzonte +72 h.
I file vanno scaricati dal portale MeteoHub, dataset ID ``WRF``, selezionando campi di superficie e livelli isobarici.
La classificazione resta prototipale e deve essere validata/calibrata su casi locali.

Installazione minima:
    pip install xarray netCDF4 numpy pandas matplotlib rasterio cfgrib eccodes shapely

Esempio con ultima corsa e intervallo temporale automatico:
    python pyroconvection_wrf15_v1_7_4.py \
        --download-wrf --wrf-run latest \
        --hours 00 03 06 09 12 15 18 21 \
        --outdir output_piroconvezione_wrf

Se ``--start`` e ``--end`` sono omessi, vengono elaborate tutte le scadenze
disponibili della corsa selezionata, limitatamente alle ore indicate da ``--hours``.
"""

from __future__ import annotations

import argparse
import gc
import json
import hashlib
import logging
import math
import os
import re
import sys
import time
from html.parser import HTMLParser
from urllib.parse import urljoin, urlparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence
from functools import lru_cache

import requests

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

OPENPOLIS_PROVINCES_URL = (
    "https://raw.githubusercontent.com/openpolis/geojson-italy/master/"
    "geojson/limits_IT_provinces.geojson"
)
DEFAULT_BOUNDARIES_CACHE_DIR = Path(".admin_boundaries")

# Livelli sufficienti per ricostruire il profilo fino a circa 8-9 km,
# evitando di aprire e scaricare livelli stratosferici non usati dal calcolo.
DEFAULT_PRESSURE_LEVELS = (
    1000, 975, 950, 925, 900, 850, 800, 750, 700, 650,
    600, 550, 500, 450, 400, 350, 300,
)
PRESSURE_VARIABLES = frozenset({"gh", "t", "q", "r", "u", "v"})


# -----------------------------------------------------------------------------
# Soglie operative iniziali, derivate dal quadro concettuale del paper
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class Thresholds:
    # Gradiente ABL: usato solo quando il layer è campionato da livelli originali.
    theta_grad_stable: float = 1.0e-3
    theta_grad_overshoot_max: float = 3.0e-3

    # Proxy quando il gradiente ABL non è realmente risolto.
    mixing_depth_proxy_m: float = 455.0

    # Cap sopra ABL.
    gamma_theta_weak_cap: float = 4.2e-3
    gamma_theta_strong_cap: float = 4.8e-3

    # Rapporti geometrici e prossimità normalizzate rispetto ad ABL.
    lcl_ratio_overshoot_max: float = 1.60
    lcl_ratio_deep_max: float = 1.10
    lcl_ratio_resilient_max: float = 1.00
    shear_distance_deep: float = 0.30

    # Umidità relativa media attorno al top ABL.
    rh_top_moist: float = 80.0

    # Limiti ABL. Sotto 600 m la colonna resta Surface plume: una ABL
    # notturna o scarsamente sviluppata non viene promossa a Convective plume.
    min_abl_m: float = 150.0
    min_classifiable_abl_m: float = 600.0
    low_confidence_abl_m: float = 600.0
    max_abl_m: float = 4500.0
    max_profile_m: float = 8000.0
    min_land_fraction: float = 0.50

    # Requisiti di risoluzione verticale sui livelli ORIGINALI, non interpolati.
    min_levels_abl_gradient: int = 3
    min_span_abl_gradient_m: float = 250.0
    min_levels_cap_gradient: int = 2
    min_span_cap_gradient_m: float = 350.0
    min_levels_shear: int = 5
    max_shear_layer_gap_m: float = 1000.0
    max_rh_bracket_gap_m: float = 900.0


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
    start: str | None
    end: str | None
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
    draw_province_boundaries: bool = True
    provinces_file: Path | None = None
    boundaries_cache_dir: Path = DEFAULT_BOUNDARIES_CACHE_DIR
    stride: int = 1


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
    pyro_class: int = NODATA_CLASS
    classification_mode: int = 0
    classification_quality: int = 0
    classification_reason: int = 0
    abl_source: int = 0
    abl_upper_bound_m: float = np.nan
    rib_crossing_resolved: int = 0
    stability_method: int = 0
    n_levels_valid: int = 0
    n_levels_abl: int = 0
    n_levels_cap: int = 0
    n_levels_shear: int = 0
    cap_resolved: int = 0
    rh_resolved: int = 0
    shear_resolved: int = 0


# classification_mode: 0 non classificabile, 1 proxy_noshear,
# 2 resolved_noshear, 3 full, 4 shallow_abl_upper_bound.
# classification_quality resta 0-3: la modalità 4 ha qualità bassa (=1) perché
# determina soltanto che ABL < 600 m, senza inventarne l'altezza esatta.
CLASSIFICATION_MODE_NAMES = {
    0: "invalid",
    1: "proxy_noshear",
    2: "resolved_noshear",
    3: "full",
    4: "shallow_abl_upper_bound",
}

# Motivo determinante della classe finale. Il codice 0 è riservato a celle
# non classificabili/NoData; gli altri codici rendono auditabile la ladder.
CLASSIFICATION_REASON_NAMES = {
    0: "invalid_or_nodata",
    1: "shallow_abl_surface_plume",
    2: "stable_or_poorly_mixed_convective_plume",
    3: "lcl_too_high_convective_plume",
    4: "overshooting_pyrocu",
    5: "insufficient_deep_mixing_convective_plume",
    6: "dry_abl_top_convective_plume",
    7: "rh_top_unresolved_resilient_pyrocu",
    8: "cap_unresolved_resilient_pyrocu",
    9: "cap_not_weak_resilient_pyrocu",
    10: "shear_not_close_resilient_pyrocu",
    11: "deep_pyrocu_noshear",
    12: "deep_pyrocu_full",
    13: "shallow_abl_upper_bound_surface_plume",
}

CLASS_NAMES = {
    0: "Surface plume",
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
# Download diretto WRF da MeteoHub (open data, senza login)
# -----------------------------------------------------------------------------

METEOHUB_WRF_BASE = "https://meteohub.agenziaitaliameteo.it/nwp/WRF/"
DEFAULT_WRF_VARIABLES = (
    "2t", "2d", "10u", "10v", "hpbl", "orog", "prmsl", "lsm",
    "gh", "t", "q", "u", "v",
)


class _HrefParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.hrefs: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() != "a":
            return
        for key, value in attrs:
            if key.lower() == "href" and value:
                self.hrefs.append(value)


def list_http_directory(url: str, session: requests.Session, timeout: int = 60) -> list[str]:
    response = session.get(url, timeout=timeout)
    response.raise_for_status()
    parser = _HrefParser()
    parser.feed(response.text)
    return [urljoin(url, href) for href in parser.hrefs if href not in {"../", "./"}]


def discover_wrf_runs(session: requests.Session, base_url: str = METEOHUB_WRF_BASE) -> list[str]:
    entries = list_http_directory(base_url, session)
    runs: list[str] = []
    for entry in entries:
        name = urlparse(entry).path.rstrip("/").split("/")[-1]
        if re.fullmatch(r"\d{10}", name):
            runs.append(name)
    return sorted(set(runs))


def resolve_wrf_run(requested: str, session: requests.Session, base_url: str = METEOHUB_WRF_BASE) -> str:
    if requested.lower() != "latest":
        if not re.fullmatch(r"\d{10}", requested):
            raise ValueError("--wrf-run deve essere 'latest' oppure YYYYMMDDHH")
        return requested
    runs = discover_wrf_runs(session, base_url)
    if not runs:
        raise RuntimeError(f"Nessuna corsa WRF trovata in {base_url}")
    return runs[-1]


def _download_file(
    url: str, destination: Path, session: requests.Session, overwrite: bool = False,
    timeout: int = 120, retries: int = 3,
) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0 and not overwrite:
        logging.info("File già presente, salto: %s", destination)
        return destination

    part = destination.with_suffix(destination.suffix + ".part")
    for attempt in range(1, retries + 1):
        try:
            headers: dict[str, str] = {}
            mode = "wb"
            existing = part.stat().st_size if part.exists() else 0
            if existing > 0 and not overwrite:
                headers["Range"] = f"bytes={existing}-"
                mode = "ab"
            elif overwrite and part.exists():
                part.unlink()

            with session.get(url, headers=headers, stream=True, timeout=timeout) as response:
                if response.status_code == 416 and part.exists():
                    part.replace(destination)
                    return destination
                response.raise_for_status()
                if existing > 0 and response.status_code != 206:
                    mode = "wb"
                    existing = 0
                total = response.headers.get("Content-Length")
                total_bytes = int(total) + existing if total and total.isdigit() else None
                downloaded = existing
                last_log = time.monotonic()
                with part.open(mode) as fh:
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if not chunk:
                            continue
                        fh.write(chunk)
                        downloaded += len(chunk)
                        if time.monotonic() - last_log >= 10:
                            if total_bytes:
                                logging.info("Download %s: %.1f%%", destination.name, 100.0 * downloaded / total_bytes)
                            else:
                                logging.info("Download %s: %.1f MB", destination.name, downloaded / 1024**2)
                            last_log = time.monotonic()
            part.replace(destination)
            logging.info("Scaricato: %s", destination)
            return destination
        except Exception as exc:
            if attempt >= retries:
                raise RuntimeError(f"Download fallito per {url}: {exc}") from exc
            logging.warning("Tentativo %d/%d fallito per %s: %s", attempt, retries, url, exc)
            time.sleep(2 ** attempt)
    raise RuntimeError(f"Download fallito per {url}")


def download_wrf_open_data(
    run: str, variables: Sequence[str], download_dir: Path, overwrite: bool = False,
    base_url: str = METEOHUB_WRF_BASE,
    pressure_levels: Sequence[int] | None = None,
) -> tuple[list[Path], str]:
    session = requests.Session()
    session.headers.update({"User-Agent": "pyroconvection-wrf15/1.7.4"})
    resolved_run = resolve_wrf_run(run, session, base_url)
    run_url = urljoin(base_url, resolved_run + "/")
    available_entries = list_http_directory(run_url, session)
    available = {urlparse(x).path.rstrip("/").split("/")[-1]: x for x in available_entries if x.endswith("/")}

    output_files: list[Path] = []
    missing_vars: list[str] = []
    for variable in variables:
        variable = variable.strip()
        if not variable:
            continue
        var_url = available.get(variable)
        if not var_url:
            missing_vars.append(variable)
            continue
        entries = list_http_directory(var_url, session)
        files = [u for u in entries if urlparse(u).path.lower().endswith((".grib", ".grb", ".grib2", ".grb2", ".nc", ".netcdf"))]
        if variable.lower() in PRESSURE_VARIABLES and pressure_levels:
            wanted = np.asarray(list(pressure_levels), dtype=float)
            filtered: list[str] = []
            for candidate in files:
                level = _pressure_level_from_path(Path(urlparse(candidate).path))
                if level is None or np.isclose(level, wanted, rtol=0.0, atol=0.05).any():
                    filtered.append(candidate)
            files = filtered
        if not files:
            logging.warning("Nessun file trovato per la variabile %s nella corsa %s", variable, resolved_run)
            continue
        for url in files:
            filename = urlparse(url).path.split("/")[-1]
            target = download_dir / resolved_run / variable / filename
            output_files.append(_download_file(url, target, session, overwrite=overwrite))

    if missing_vars:
        logging.warning("Variabili non disponibili nella corsa %s: %s", resolved_run, ", ".join(missing_vars))
    if not output_files:
        raise RuntimeError("Il download WRF non ha prodotto alcun file")
    manifest = download_dir / resolved_run / "download_manifest.txt"
    manifest.write_text("\n".join(str(x) for x in output_files) + "\n", encoding="utf-8")
    return output_files, resolved_run


# -----------------------------------------------------------------------------
# Input WRF da MeteoHub
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


def _cfgrib_index_path(path: Path, index_dir: Path | None) -> str:
    """Restituisce un percorso stabile per l'indice cfgrib/eccodes.

    Nella versione precedente ``indexpath=""`` disabilitava completamente gli
    indici. Ogni esecuzione doveva quindi rileggere e indicizzare da zero tutti
    i messaggi GRIB. Un indice persistente riduce sensibilmente i tempi delle
    esecuzioni successive.
    """
    if index_dir is None:
        return ""
    index_dir.mkdir(parents=True, exist_ok=True)
    signature = hashlib.sha1(str(path.resolve()).encode("utf-8")).hexdigest()[:12]
    return str(index_dir / f"{path.name}.{signature}.idx")


def open_one_dataset(
    path: Path,
    cfgrib_index_dir: Path | None = None,
) -> xr.Dataset:
    suffix = path.suffix.lower()
    if suffix in {".grib", ".grb", ".grib2", ".grb2"}:
        try:
            indexpath = _cfgrib_index_path(path, cfgrib_index_dir)
            return xr.open_dataset(
                path,
                engine="cfgrib",
                cache=False,
                backend_kwargs={"indexpath": indexpath},
            )
        except Exception as exc:
            raise RuntimeError(
                f"Impossibile aprire {path} con cfgrib. Installare cfgrib/eccodes o usare NetCDF: {exc}"
            ) from exc
    return xr.open_dataset(path, decode_times=True, cache=False)

# -----------------------------------------------------------------------------
# Apertura e normalizzazione dataset WRF
# -----------------------------------------------------------------------------

def _rename_common(ds: xr.Dataset) -> xr.Dataset:
    rename: dict[str, str] = {}
    candidates = {
        "valid_time": "time",
        "forecast_time": "time",
        "pressure_level": "level",
        "pressureLevel": "level",
        "isobaricInhPa": "level",
        "isobaricInPa": "level",
        "isobaric": "level",
        "plev": "level",
        "lat": "latitude",
        "lon": "longitude",
        "XLAT": "latitude",
        "XLONG": "longitude",
        "XLAT_M": "latitude",
        "XLONG_M": "longitude",
    }
    for old, new in candidates.items():
        if old in ds.dims or old in ds.coords or old in ds.data_vars:
            if new not in ds.dims and new not in ds.coords and new not in ds.data_vars:
                rename[old] = new
    if rename:
        ds = ds.rename(rename)

    # Nei NetCDF WRF XLAT/XLONG possono essere normali variabili anziché
    # coordinate. Le promuoviamo a coordinate per conservarle durante merge,
    # transpose ed esportazione.
    to_coords = [name for name in ("latitude", "longitude") if name in ds.data_vars]
    if to_coords:
        ds = ds.set_coords(to_coords)
    return ds


_ISOBARIC_FILENAME_RE = re.compile(
    r"(?:^|[_-])isobaricIn(?:hPa|Pa)[-_](?P<level>\d+(?:\.\d+)?)",
    flags=re.IGNORECASE,
)


def _pressure_level_from_path(path: Path) -> float | None:
    """Ricava il livello isobarico in hPa dal nome del file MeteoHub."""
    match = _ISOBARIC_FILENAME_RE.search(path.name)
    if not match:
        return None
    value = float(match.group("level"))
    # I nomi isobaricInPa usano Pa, quelli isobaricInhPa usano hPa.
    if "isobaricinpa" in path.name.lower() and "isobaricinhpa" not in path.name.lower():
        value /= 100.0
    return value


def _make_known_static_fields_timeless(ds: xr.Dataset, path: Path) -> xr.Dataset:
    """Rimuove il tempo soltanto dai campi WRF statici.

    Orografia e land mask possono essere distribuite con una sola scadenza. Nei
    NetCDF multidimensionali possono però convivere con T2, vento e altri campi
    dinamici: in tal caso non deve essere eliminata la dimensione temporale
    dell'intero dataset.
    """
    static_names = {
        "orog", "hgt", "terrain_height", "surface_geopotential",
        "z_sfc", "lsm", "landmask", "land_sea_mask",
    }
    parent_name = path.parent.name.lower()
    data_names = list(ds.data_vars)
    static_vars = [name for name in data_names if str(name).lower() in static_names]

    # Un file collocato nella directory di una variabile statica è normalmente
    # monovariabile; in caso contrario si interviene comunque solo sui campi noti.
    if parent_name in static_names and len(data_names) == 1:
        static_vars = data_names

    for name in static_vars:
        da = ds[name]
        if "time" in da.dims:
            timeless = da.isel(time=0, drop=True)
            ds = ds.drop_vars(name)
            ds[name] = timeless

    # Se nessuna variabile usa più time, elimina la coordinata residua.
    if "time" in ds.coords and not any("time" in ds[name].dims for name in ds.data_vars):
        ds = ds.drop_vars("time")
    return ds


def _normalize_pressure_level(ds: xr.Dataset, path: Path) -> xr.Dataset:
    """Trasforma il livello isobarico scalare in una vera dimensione ``level``.

    I file WRF di MeteoHub sono separati per variabile e livello. Cfgrib
    espone normalmente ``isobaricInhPa`` come coordinata scalare; senza questa
    promozione un merge tra 1000, 925, 850 hPa ecc. sovrascrive la coordinata
    invece di costruire il profilo verticale.
    """
    level_from_name = _pressure_level_from_path(path)

    if "level" in ds.coords or "level" in ds.data_vars:
        level_da = ds["level"]
        values = np.asarray(level_da.values, dtype=float).reshape(-1)
        if values.size == 0:
            return ds

        units = str(level_da.attrs.get("units", "")).lower().replace(" ", "")
        values_hpa = values / 100.0 if units in {"pa", "pascal", "pascals"} else values
        # Euristica aggiuntiva: un livello > 2000 è quasi certamente espresso in Pa.
        values_hpa = np.where(values_hpa > 2000.0, values_hpa / 100.0, values_hpa)

        if "level" in ds.dims:
            if ds.sizes["level"] != values_hpa.size:
                raise ValueError(
                    f"Dimensione level incoerente in {path}: "
                    f"size={ds.sizes['level']}, valori={values_hpa.size}."
                )
            ds = ds.assign_coords(level=("level", values_hpa.astype(float)))
        elif values_hpa.size == 1:
            # Rimuove la coordinata scalare prima di creare l'omonima dimensione.
            ds = ds.drop_vars("level")
            ds = ds.expand_dims(level=[float(values_hpa[0])])
        else:
            raise ValueError(
                f"La coordinata level in {path} contiene {values_hpa.size} valori "
                "ma non è una dimensione."
            )
        ds["level"].attrs.update({"units": "hPa", "positive": "down"})
        return ds

    if level_from_name is not None:
        ds = ds.expand_dims(level=[float(level_from_name)])
        ds["level"].attrs.update({"units": "hPa", "positive": "down"})
    return ds


def _merge_same_level(datasets: Sequence[xr.Dataset], level: float) -> xr.Dataset:
    """Unisce variabili diverse appartenenti allo stesso livello isobarico."""
    pieces: list[xr.Dataset] = []
    for ds in datasets:
        if "level" in ds.dims and ds.sizes["level"] == 1:
            pieces.append(ds.isel(level=0, drop=True))
        else:
            pieces.append(ds)
    try:
        merged = xr.merge(
            pieces,
            compat="override",
            join="exact",
            combine_attrs="override",
        )
    except ValueError as exc:
        logging.warning(
            "Coordinate non perfettamente coincidenti al livello %g hPa (%s); "
            "uso join esterno", level, exc
        )
        merged = xr.merge(
            pieces,
            compat="override",
            join="outer",
            combine_attrs="override",
        )
    merged = merged.expand_dims(level=[float(level)])
    merged["level"].attrs.update({"units": "hPa", "positive": "down"})
    return merged


def _combine_opened_datasets(datasets: Sequence[xr.Dataset]) -> xr.Dataset:
    """Combina dataset WRF senza perdere i livelli isobarici scalari."""
    if len(datasets) == 1:
        return datasets[0]

    with_level = [ds for ds in datasets if "level" in ds.dims]
    without_level = [ds for ds in datasets if "level" not in ds.dims]

    combined_parts: list[xr.Dataset] = []
    if with_level:
        if all(ds.sizes.get("level", 0) == 1 for ds in with_level):
            groups: dict[float, list[xr.Dataset]] = {}
            for ds in with_level:
                level = float(np.asarray(ds["level"].values).reshape(-1)[0])
                groups.setdefault(level, []).append(ds)
            per_level = [
                _merge_same_level(groups[level], level)
                for level in sorted(groups, reverse=True)
            ]
            try:
                pressure = xr.concat(
                    per_level,
                    dim="level",
                    data_vars="minimal",
                    coords="minimal",
                    compat="override",
                    join="exact",
                    combine_attrs="override",
                )
            except ValueError as exc:
                logging.warning(
                    "Coordinate dei livelli non perfettamente coincidenti (%s); "
                    "uso join esterno", exc
                )
                pressure = xr.concat(
                    per_level,
                    dim="level",
                    data_vars="minimal",
                    coords="minimal",
                    compat="override",
                    join="outer",
                    combine_attrs="override",
                )
        else:
            try:
                pressure = xr.combine_by_coords(
                    with_level,
                    data_vars="all",
                    coords="minimal",
                    compat="override",
                    join="outer",
                    combine_attrs="override",
                )
            except Exception as exc:
                logging.warning(
                    "Combinazione dei livelli isobarici non riuscita (%s); "
                    "provo merge esterno", exc
                )
                pressure = xr.merge(
                    with_level, compat="override", join="outer",
                    combine_attrs="override"
                )
        combined_parts.append(pressure)

    if without_level:
        # I file di superficie contengono normalmente variabili differenti sulle
        # stesse coordinate. Non devono essere concatenati: un merge diretto è
        # più rapido e non genera l'avviso inevitabile di combine_by_coords.
        try:
            surface = xr.merge(
                without_level, compat="override", join="exact",
                combine_attrs="override"
            )
        except ValueError as exc:
            logging.warning(
                "Coordinate dei campi di superficie non perfettamente coincidenti "
                "(%s); uso join esterno", exc
            )
            surface = xr.merge(
                without_level, compat="override", join="outer",
                combine_attrs="override"
            )
        combined_parts.append(surface)

    if len(combined_parts) == 1:
        return combined_parts[0]
    return xr.merge(
        combined_parts, compat="override", join="outer",
        combine_attrs="override"
    )


def split_surface_pressure_files(paths: Sequence[Path]) -> tuple[list[Path], list[Path]]:
    """Separa i file MeteoHub di superficie da quelli isobarici.

    Per un singolo NetCDF/GRIB già multidimensionale, o quando il nome non
    consente la classificazione, mantiene il comportamento compatibile usando
    lo stesso insieme per entrambi i dataset.
    """
    all_paths = list(paths)
    pressure = [p for p in all_paths if _pressure_level_from_path(p) is not None]
    surface = [p for p in all_paths if _pressure_level_from_path(p) is None]
    if pressure and surface:
        logging.info(
            "File WRF separati automaticamente: %d superficie, %d isobarici",
            len(surface), len(pressure)
        )
        return surface, pressure
    return all_paths, all_paths


def filter_pressure_files(
    paths: Sequence[Path],
    requested_levels: Sequence[int] | None,
) -> list[Path]:
    """Mantiene solo i livelli isobarici necessari al diagnostico.

    I file multidimensionali, per i quali il livello non è ricavabile dal nome,
    vengono mantenuti e saranno selezionati successivamente sul dataset.
    """
    all_paths = list(paths)
    if not requested_levels:
        return all_paths
    wanted = np.asarray(sorted(set(float(x) for x in requested_levels)), dtype=float)
    kept: list[Path] = []
    skipped = 0
    for path in all_paths:
        level = _pressure_level_from_path(path)
        if level is None or np.isclose(level, wanted, rtol=0.0, atol=0.05).any():
            kept.append(path)
        else:
            skipped += 1
    if skipped:
        logging.info(
            "Filtro livelli isobarici: apro %d file e ne escludo %d non necessari",
            len(kept), skipped
        )
    return kept


def select_pressure_levels(
    ds: xr.Dataset,
    requested_levels: Sequence[int] | None,
) -> xr.Dataset:
    if not requested_levels or "level" not in ds.coords:
        return ds
    available = np.asarray(ds["level"].values, dtype=float)
    wanted = np.asarray(sorted(set(float(x) for x in requested_levels), reverse=True))
    selected = [
        float(level) for level in available
        if np.isclose(level, wanted, rtol=0.0, atol=0.05).any()
    ]
    if not selected:
        raise ValueError(
            f"Nessuno dei livelli richiesti {list(requested_levels)} è presente. "
            f"Disponibili: {available.tolist()}"
        )
    return ds.sel(level=selected)


def _drop_or_combine_expver(ds: xr.Dataset) -> xr.Dataset:
    if "expver" not in ds.dims:
        return ds
    # ERA5 può avere due versioni complementari. combine_first mantiene i dati validi.
    pieces = [ds.isel(expver=i, drop=True) for i in range(ds.sizes["expver"])]
    out = pieces[0]
    for p in pieces[1:]:
        out = out.combine_first(p)
    return out


def _normalize_forecast_time(ds: xr.Dataset) -> xr.Dataset:
    """Normalizza i tempi GRIB cfgrib usando il tempo valido come dimensione ``time``.

    Nei GRIB previsionali cfgrib espone spesso:
    - ``time``: tempo iniziale della corsa (coordinata scalare);
    - ``step``: scadenza previsionale (dimensione);
    - ``valid_time``: tempo valido, generalmente dipendente da ``step``.

    Rinominare direttamente ``valid_time`` in ``time`` genera un conflitto con la
    coordinata scalare ``time``. La coordinata della corsa viene quindi conservata
    come ``forecast_reference_time`` e il tempo valido diventa la vera dimensione
    temporale del dataset.
    """
    if "valid_time" in ds.coords:
        valid = ds["valid_time"]

        # Conserva il tempo di inizializzazione della corsa senza occupare il nome time.
        if "time" in ds.coords and "time" not in ds.dims:
            if "forecast_reference_time" not in ds.coords:
                ds = ds.assign_coords(forecast_reference_time=ds["time"])
            ds = ds.drop_vars("time")

        if valid.ndim == 0:
            # Singola scadenza: crea esplicitamente la dimensione temporale.
            value = np.asarray(valid.values).reshape(-1)[0]
            ds = ds.drop_vars("valid_time")
            if "time" in ds.dims:
                return ds
            return ds.expand_dims(time=[value])

        valid_dim = valid.dims[0]
        if valid_dim != "valid_time":
            ds = ds.swap_dims({valid_dim: "valid_time"})
            if valid_dim in ds.coords:
                ds = ds.drop_vars(valid_dim)
        if "time" in ds.coords and "time" not in ds.dims:
            ds = ds.drop_vars("time")
        return ds.rename({"valid_time": "time"})

    # Fallback per GRIB che non espongono valid_time ma hanno time + step.
    if "step" in ds.coords and "time" in ds.coords:
        reference = ds["time"]
        valid = reference + ds["step"]
        if "time" not in ds.dims:
            if "forecast_reference_time" not in ds.coords:
                ds = ds.assign_coords(forecast_reference_time=reference)
            ds = ds.drop_vars("time")
        ds = ds.assign_coords(valid_time=valid)
        if "step" in ds.dims:
            ds = ds.swap_dims({"step": "valid_time"})
        if "step" in ds.coords:
            ds = ds.drop_vars("step")
        return ds.rename({"valid_time": "time"})

    # NetCDF già normalizzato o singolo tempo scalare.
    if "time" in ds.coords and "time" not in ds.dims:
        value = np.asarray(ds["time"].values).reshape(-1)[0]
        ds = ds.drop_vars("time")
        ds = ds.expand_dims(time=[value])
    return ds


def open_many_netcdf(
    paths: Sequence[Path],
    *,
    start: str | None = None,
    end: str | None = None,
    hours: Sequence[str] | None = None,
    domain: Domain | None = None,
    stride: int = 1,
    cfgrib_index_dir: Path | None = None,
) -> xr.Dataset:
    missing = [str(p) for p in paths if not p.exists()]
    if missing:
        raise FileNotFoundError("File mancanti:\n" + "\n".join(missing))

    if stride < 1:
        raise ValueError("stride deve essere >= 1")

    datasets: list[xr.Dataset] = []
    spatial_indexers: dict[str, slice] | None = None
    spatial_grid_signature: tuple[tuple[str, int], ...] | None = None

    for p in paths:
        logging.info("Apro %s", p)
        ds = open_one_dataset(p, cfgrib_index_dir=cfgrib_index_dir)
        try:
            ds = _rename_common(ds)
            ds = _drop_or_combine_expver(ds)
            ds = _normalize_pressure_level(ds, p)
            ds = _normalize_forecast_time(ds)
            ds = _make_known_static_fields_timeless(ds, p)

            # Riduce il grafo xarray e la memoria prima del merge/concat: vengono
            # mantenute soltanto le scadenze realmente richieste.
            if hours is not None and "time" in ds.coords:
                ds = select_requested_time(
                    ds, start=start, end=end, hours=hours, allow_empty=True
                )
                if ds.sizes.get("time", 0) == 0:
                    logging.debug("Nessuna scadenza richiesta in %s; salto", p)
                    ds.close()
                    continue

            # Ottimizzazione principale: bbox e stride vengono applicati a ogni
            # singolo file PRIMA che entri nel concat/merge. In questo modo le
            # operazioni xarray lavorano già sul solo dominio richiesto.
            if domain is not None or stride > 1:
                effective_domain = domain or Domain(90.0, -180.0, -90.0, 180.0)
                # Ricava dimensioni e forma senza materializzare latitude/longitude.
                # I valori geografici vengono letti soltanto quando è necessario
                # calcolare un nuovo ritaglio per una griglia differente.
                lat_da = _spatial_coord(ds["latitude"])
                lon_da = _spatial_coord(ds["longitude"])
                if lat_da.ndim == 1 and lon_da.ndim == 1:
                    y_dim, x_dim = lat_da.dims[0], lon_da.dims[0]
                elif lat_da.ndim == 2 and lon_da.ndim == 2:
                    y_dim, x_dim = lat_da.dims
                else:
                    raise ValueError(
                        "Latitude e longitude devono essere entrambe 1D oppure "
                        f"entrambe 2D; forme trovate: {lat_da.shape} e {lon_da.shape}."
                    )
                current_signature = (
                    (y_dim, int(ds.sizes[y_dim])),
                    (x_dim, int(ds.sizes[x_dim])),
                )
                if spatial_indexers is None or current_signature != spatial_grid_signature:
                    spatial_indexers = geographic_subset_indexers(
                        ds, effective_domain, stride=stride
                    )
                    spatial_grid_signature = current_signature
                    before_shape = (int(ds.sizes[y_dim]), int(ds.sizes[x_dim]))
                    after_shape = (
                        len(range(*spatial_indexers[y_dim].indices(before_shape[0]))),
                        len(range(*spatial_indexers[x_dim].indices(before_shape[1]))),
                    )
                    logging.info(
                        "Ritaglio pre-merge: griglia %dx%d -> %dx%d, stride=%d, "
                        "bbox N=%.4f W=%.4f S=%.4f E=%.4f",
                        before_shape[0], before_shape[1],
                        after_shape[0], after_shape[1], stride,
                        effective_domain.north, effective_domain.west,
                        effective_domain.south, effective_domain.east,
                    )
                ds = ds.isel(spatial_indexers)

            datasets.append(ds)
        except Exception:
            ds.close()
            raise

    if not datasets:
        raise ValueError("Nessun dataset contiene le scadenze temporali richieste.")

    out = _combine_opened_datasets(datasets)

    if "time" not in out.coords:
        raise ValueError("Coordinata temporale non trovata nei file WRF.")
    if "time" not in out.dims:
        value = np.asarray(out["time"].values).reshape(-1)[0]
        out = out.drop_vars("time").expand_dims(time=[value])

    out = out.sortby("time")
    _, unique_idx = np.unique(out["time"].values, return_index=True)
    out = out.isel(time=np.sort(unique_idx))
    return out

def select_requested_time(
    ds: xr.Dataset,
    start: str | None,
    end: str | None,
    hours: Sequence[str] | None,
    allow_empty: bool = False,
) -> xr.Dataset:
    """Seleziona intervallo e ore valide, consentendo un range automatico.

    Quando ``start`` e ``end`` sono entrambi ``None`` non viene applicato alcun
    limite di data: restano tutte le scadenze disponibili, eventualmente
    filtrate dalle ore UTC richieste. La CLI impone che start/end siano forniti
    insieme, ma la funzione gestisce anche limiti singoli per robustezza.
    """
    selected = ds

    start_ts = pd.Timestamp(start) if start is not None else None
    end_ts = (
        pd.Timestamp(end) + pd.Timedelta(days=1) - pd.Timedelta(microseconds=1)
        if end is not None else None
    )
    if start_ts is not None or end_ts is not None:
        selected = selected.sel(time=slice(start_ts, end_ts))

    if hours is not None:
        wanted_hours = {int(h.split(":")[0]) for h in hours}
        mask = xr.DataArray(
            np.array([pd.Timestamp(t).hour in wanted_hours for t in selected.time.values]),
            dims=("time",),
            coords={"time": selected.time},
        )
        selected = selected.sel(time=mask)

    if selected.sizes.get("time", 0) == 0 and not allow_empty:
        if start is None and end is None:
            detail = f"nelle ore UTC {list(hours or [])}"
        else:
            detail = f"tra {start or '-∞'} e {end or '+∞'} alle ore UTC {list(hours or [])}"
        raise ValueError(f"Nessun tempo WRF selezionato {detail}.")
    return selected


def find_var(ds: xr.Dataset, candidates: Sequence[str]) -> xr.DataArray:
    for name in candidates:
        if name in ds.data_vars:
            return ds[name]
    raise KeyError(
        f"Nessuna variabile tra {list(candidates)}. Presenti: {list(ds.data_vars)}"
    )




def find_optional_var(ds: xr.Dataset, candidates: Sequence[str]) -> xr.DataArray | None:
    for name in candidates:
        if name in ds.data_vars:
            return ds[name]
    return None

def convert_temperature_to_k(da: xr.DataArray) -> xr.DataArray:
    units = str(da.attrs.get("units", "")).lower()
    if units in {"c", "degc", "degree_celsius", "degrees_celsius", "°c"}:
        return da + 273.15
    return da

def convert_pressure_to_pa(da: xr.DataArray) -> xr.DataArray:
    units = str(da.attrs.get("units", "")).lower()
    if units in {"hpa", "mb", "mbar"}:
        return da * 100.0
    return da

def height_or_geopotential_to_geopotential(da: xr.DataArray) -> xr.DataArray:
    units = str(da.attrs.get("units", "")).lower().replace(" ", "")
    name = str(da.name or "").lower()
    if units in {"m", "meter", "metre", "meters", "metres", "gpm"} or "height" in name or name in {"gh", "hgt", "orog", "ter", "terrain"}:
        return da * G
    return da

def dewpoint_from_t_rh(t_k: xr.DataArray, rh_percent: xr.DataArray) -> xr.DataArray:
    tc = t_k - 273.15
    rh = xr.where(rh_percent > 1.5, rh_percent / 100.0, rh_percent).clip(min=1e-4, max=1.0)
    a, b = 17.625, 243.04
    gamma = np.log(rh) + a * tc / (b + tc)
    return (b * gamma / (a - gamma)) + 273.15

def specific_humidity_from_rh_da(rh_percent: xr.DataArray, t_k: xr.DataArray, level_hpa: xr.DataArray) -> xr.DataArray:
    rh = xr.where(rh_percent > 1.5, rh_percent / 100.0, rh_percent).clip(min=0.0, max=1.0)
    es = 611.2 * np.exp(17.67 * (t_k - 273.15) / (t_k - 29.65))
    p = level_hpa * 100.0
    e = rh * es
    return EPSILON * e / (p - (1.0 - EPSILON) * e)


def _spatial_coord(da: xr.DataArray) -> xr.DataArray:
    """Riduce una coordinata geografica alle sole dimensioni orizzontali.

    Alcuni NetCDF WRF memorizzano XLAT/XLONG anche con una dimensione temporale
    o con ulteriori dimensioni singleton. La griglia è statica: in questi casi
    si usa il primo elemento.
    """
    out = da
    while out.ndim > 2:
        removable = None
        for dim in out.dims:
            if dim == "time" or out.sizes[dim] == 1:
                removable = dim
                break
        if removable is None:
            raise ValueError(
                f"La coordinata {da.name!r} ha dimensioni non gestibili {da.dims}; "
                "sono supportate coordinate geografiche 1D o 2D statiche."
            )
        out = out.isel({removable: 0}, drop=True)
    return out


def horizontal_grid_info(
    obj: xr.Dataset | xr.DataArray,
) -> tuple[str, str, np.ndarray, np.ndarray, bool]:
    """Restituisce dimensioni Y/X e coordinate lat/lon della griglia.

    L'ultimo valore indica se la griglia è curvilinea (coordinate 2D).
    """
    if "latitude" not in obj.coords or "longitude" not in obj.coords:
        raise ValueError(
            "Coordinate latitude/longitude non trovate. "
            f"Coordinate disponibili: {list(obj.coords)}"
        )

    lat_da = _spatial_coord(obj["latitude"])
    lon_da = _spatial_coord(obj["longitude"])

    if lat_da.ndim == 1 and lon_da.ndim == 1:
        return (
            lat_da.dims[0],
            lon_da.dims[0],
            np.asarray(lat_da.values),
            np.asarray(lon_da.values),
            False,
        )

    if lat_da.ndim == 2 and lon_da.ndim == 2:
        if set(lat_da.dims) != set(lon_da.dims):
            raise ValueError(
                "Latitude e longitude 2D usano dimensioni incompatibili: "
                f"{lat_da.dims} e {lon_da.dims}."
            )
        if lon_da.dims != lat_da.dims:
            lon_da = lon_da.transpose(*lat_da.dims)
        y_dim, x_dim = lat_da.dims
        return (
            y_dim,
            x_dim,
            np.asarray(lat_da.values),
            np.asarray(lon_da.values),
            True,
        )

    raise ValueError(
        "Latitude e longitude devono essere entrambe 1D oppure entrambe 2D; "
        f"forme trovate: {lat_da.shape} e {lon_da.shape}."
    )


def _finite_edge_mean(values: np.ndarray, axis: int, index: int) -> float:
    edge = np.take(np.asarray(values, dtype=float), index, axis=axis)
    if not np.isfinite(edge).any():
        return np.nan
    return float(np.nanmean(edge))


def geographic_subset_indexers(
    ds: xr.Dataset,
    domain: Domain,
    *,
    stride: int = 1,
) -> dict[str, slice]:
    """Calcola gli indexer spaziali per bbox e stride senza modificare il dataset.

    Lo stride viene applicato *dopo* il ritaglio geografico. Gli indexer possono
    essere riutilizzati per tutti i file della stessa griglia, evitando di
    scandire ripetutamente latitude/longitude prima del merge.
    """
    if stride < 1:
        raise ValueError("stride deve essere >= 1")

    y_dim, x_dim, lat, lon, curvilinear = horizontal_grid_info(ds)
    global_domain = (
        domain.north >= 89.999 and domain.south <= -89.999
        and domain.west <= -179.999 and domain.east >= 179.999
    )

    if global_domain:
        return {
            y_dim: slice(0, int(ds.sizes[y_dim]), stride),
            x_dim: slice(0, int(ds.sizes[x_dim]), stride),
        }

    if curvilinear:
        mask = (
            np.isfinite(lat) & np.isfinite(lon)
            & (lat >= domain.south) & (lat <= domain.north)
            & (lon >= domain.west) & (lon <= domain.east)
        )
        if not mask.any():
            raise ValueError(
                "Il bounding box richiesto non interseca la griglia WRF: "
                f"N={domain.north}, W={domain.west}, S={domain.south}, E={domain.east}."
            )
        rows, cols = np.where(mask)
        return {
            y_dim: slice(int(rows.min()), int(rows.max()) + 1, stride),
            x_dim: slice(int(cols.min()), int(cols.max()) + 1, stride),
        }

    y_idx = np.where(
        np.isfinite(lat) & (lat >= domain.south) & (lat <= domain.north)
    )[0]
    x_idx = np.where(
        np.isfinite(lon) & (lon >= domain.west) & (lon <= domain.east)
    )[0]
    if y_idx.size == 0 or x_idx.size == 0:
        raise ValueError(
            "Il bounding box richiesto non interseca la griglia WRF: "
            f"N={domain.north}, W={domain.west}, S={domain.south}, E={domain.east}."
        )
    return {
        y_dim: slice(int(y_idx.min()), int(y_idx.max()) + 1, stride),
        x_dim: slice(int(x_idx.min()), int(x_idx.max()) + 1, stride),
    }


def subset_geographic_domain(
    ds: xr.Dataset,
    domain: Domain,
    *,
    stride: int = 1,
) -> xr.Dataset:
    """Ritaglia una griglia rettilinea o curvilinea e applica lo stride."""
    return ds.isel(geographic_subset_indexers(ds, domain, stride=stride))

def ensure_lat_lon_order(ds: xr.Dataset) -> xr.Dataset:
    """Orienta la griglia da nord a sud e da ovest a est.

    Funziona sia con coordinate rettilinee 1D sia con le coordinate
    curvilinee 2D tipiche dei GRIB WRF. Il ribaltamento viene applicato alle
    dimensioni della griglia con ``isel`` e quindi coinvolge coerentemente
    tutte le variabili del dataset.
    """
    y_dim, x_dim, lat, lon, curvilinear = horizontal_grid_info(ds)

    if curvilinear:
        lat_first = _finite_edge_mean(lat, axis=0, index=0)
        lat_last = _finite_edge_mean(lat, axis=0, index=-1)
    else:
        lat_first = float(np.asarray(lat, dtype=float)[0])
        lat_last = float(np.asarray(lat, dtype=float)[-1])

    if np.isfinite(lat_first) and np.isfinite(lat_last) and lat_first < lat_last:
        logging.debug(
            "Ribalto la dimensione %s: latitudine da sud-nord a nord-sud", y_dim
        )
        ds = ds.isel({y_dim: slice(None, None, -1)})

    # Rilegge la griglia dopo l'eventuale ribaltamento nord-sud.
    y_dim, x_dim, lat, lon, curvilinear = horizontal_grid_info(ds)
    if curvilinear:
        lon_first = _finite_edge_mean(lon, axis=1, index=0)
        lon_last = _finite_edge_mean(lon, axis=1, index=-1)
    else:
        lon_first = float(np.asarray(lon, dtype=float)[0])
        lon_last = float(np.asarray(lon, dtype=float)[-1])

    if np.isfinite(lon_first) and np.isfinite(lon_last) and lon_first > lon_last:
        logging.debug(
            "Ribalto la dimensione %s: longitudine da est-ovest a ovest-est", x_dim
        )
        ds = ds.isel({x_dim: slice(None, None, -1)})

    return ds


def _same_grid(ds_a: xr.Dataset, ds_b: xr.Dataset, atol: float = 1.0e-6) -> bool:
    _, _, lat_a, lon_a, _ = horizontal_grid_info(ds_a)
    _, _, lat_b, lon_b, _ = horizontal_grid_info(ds_b)
    return (
        lat_a.shape == lat_b.shape
        and lon_a.shape == lon_b.shape
        and np.allclose(lat_a, lat_b, rtol=0.0, atol=atol, equal_nan=True)
        and np.allclose(lon_a, lon_b, rtol=0.0, atol=atol, equal_nan=True)
    )


def align_datasets(ds_sl: xr.Dataset, ds_pl: xr.Dataset) -> tuple[xr.Dataset, xr.Dataset]:
    ds_sl = ensure_lat_lon_order(ds_sl)
    ds_pl = ensure_lat_lon_order(ds_pl)

    if "level" not in ds_pl.dims:
        available = ", ".join(ds_pl.dims) or "nessuna"
        raise ValueError(
            "La dimensione dei livelli isobarici 'level' non è stata costruita. "
            f"Dimensioni disponibili: {available}. "
            "Verificare di avere incluso i file isobaricInhPa-* delle variabili "
            "gh, t, q, u e v."
        )

    # Pressione ordinata dal basso verso l'alto: 1000 -> 100 hPa.
    level_values = np.asarray(ds_pl["level"].values).reshape(-1)
    if level_values.size > 1 and level_values[0] < level_values[-1]:
        ds_pl = ds_pl.sortby("level", ascending=False)

    common_time = np.intersect1d(ds_sl.time.values, ds_pl.time.values)
    if len(common_time) == 0:
        raise ValueError("Nessun timestamp comune fra single levels e pressure levels.")

    ds_sl = ds_sl.sel(time=common_time)
    ds_pl = ds_pl.sel(time=common_time)

    if _same_grid(ds_sl, ds_pl):
        return ds_sl, ds_pl

    sl_y, sl_x, sl_lat, sl_lon, sl_curv = horizontal_grid_info(ds_sl)
    pl_y, pl_x, pl_lat, pl_lon, pl_curv = horizontal_grid_info(ds_pl)

    # Interpolazione xarray ammessa solo per una griglia rettilinea con
    # latitude/longitude vere dimensioni indicizzate.
    rectilinear_indexed = (
        not sl_curv
        and not pl_curv
        and sl_y == "latitude"
        and sl_x == "longitude"
        and pl_y == "latitude"
        and pl_x == "longitude"
    )
    if rectilinear_indexed:
        logging.warning(
            "Griglie rettilinee non coincidenti: interpolo i pressure levels "
            "sulla griglia single levels."
        )
        ds_pl = ds_pl.interp(
            latitude=ds_sl["latitude"],
            longitude=ds_sl["longitude"],
        )
        return ds_sl, ds_pl

    details = (
        f"single={sl_lat.shape}/{sl_y},{sl_x}; "
        f"pressure={pl_lat.shape}/{pl_y},{pl_x}"
    )
    raise ValueError(
        "Le griglie WRF curvilinee di superficie e dei livelli isobarici non "
        "coincidono e non possono essere interpolate con xarray.interp. "
        "Usare file appartenenti alla stessa corsa e allo stesso dominio. "
        + details
    )


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


def first_bracketed_threshold_crossing(
    z: np.ndarray,
    y: np.ndarray,
    threshold: float,
    zmin: float,
) -> tuple[float, bool, float]:
    """Diagnostica il primo superamento della soglia sui livelli originali.

    Restituisce ``(crossing, bracketed, first_exceeding_height)``.

    - Se il crossing è delimitato da un livello sotto soglia e uno sopra soglia,
      ``crossing`` contiene la quota interpolata e ``bracketed=True``.
    - Se il primo livello disponibile sopra ``zmin`` supera già la soglia, la
      quota esatta non è risolta. In questo caso ``crossing=nan``,
      ``bracketed=False`` e ``first_exceeding_height`` costituisce soltanto un
      limite superiore: il vero crossing si trova a una quota inferiore o uguale.
    - Se non esiste alcun superamento utilizzabile, tutti i valori diagnostici
      sono non validi.
    """
    mask = np.isfinite(z) & np.isfinite(y) & (z >= zmin)
    zz = np.asarray(z[mask], dtype=float)
    yy = np.asarray(y[mask], dtype=float)
    if len(zz) < 2:
        return np.nan, False, np.nan

    order = np.argsort(zz)
    zz = zz[order]
    yy = yy[order]
    idx = np.where(yy >= threshold)[0]
    if len(idx) == 0:
        return np.nan, False, np.nan

    i = int(idx[0])
    first_exceeding_height = float(zz[i])

    # Il primo livello utile è già sopra soglia: non si attribuisce una quota
    # puntuale all'ABL, ma si conserva l'informazione fisicamente robusta che il
    # crossing è avvenuto non oltre questo livello.
    if i == 0:
        return np.nan, False, first_exceeding_height

    z0, z1 = float(zz[i - 1]), float(zz[i])
    y0, y1 = float(yy[i - 1]), float(yy[i])
    if not np.isfinite([z0, z1, y0, y1]).all() or y0 >= threshold or y1 == y0:
        return np.nan, False, first_exceeding_height

    frac = (threshold - y0) / (y1 - y0)
    crossing = float(z0 + np.clip(frac, 0.0, 1.0) * (z1 - z0))
    return crossing, True, first_exceeding_height


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
) -> tuple[float, int, int, float]:
    """Stima adattiva dell'ABL dal Bulk Richardson number.

    Il crossing ``RiB >= 0.33`` viene accettato come altezza puntuale solo se è
    verticalmente delimitato da due livelli originali. In caso contrario:

    1. viene usato PBLH/HPBL del WRF, se finito e positivo;
    2. se PBLH non è valida ma il primo livello già sopra soglia si trova a una
       quota <= ``min_classifiable_abl_m``, non viene inventata un'altezza ABL:
       si restituisce soltanto un limite superiore sufficiente a concludere che
       la colonna deve restare ``Surface plume``;
    3. negli altri casi il profilo rimane non classificabile.

    Restituisce ``(altezza, fonte, crossing_rib_risolto, limite_superiore)``.

    Codici ``fonte``:
      0 = non disponibile;
      1 = crossing Bulk Richardson delimitato;
      2 = PBLH/HPBL WRF;
      3 = solo limite superiore RiB, certamente sotto la soglia operativa ABL.
    """
    def valid_pblh() -> bool:
        # Mantiene esattamente il comportamento della v1.7.2 per tutti i casi
        # in cui PBLH è disponibile: qualsiasi valore finito e positivo viene
        # accettato e poi limitato all'intervallo operativo.
        return bool(np.isfinite(fallback_blh) and fallback_blh > 0.0)

    def fallback(first_exceeding_height: float = np.nan) -> tuple[float, int, int, float]:
        if valid_pblh():
            return (
                float(np.clip(fallback_blh, thresholds.min_abl_m, thresholds.max_abl_m)),
                2,
                0,
                first_exceeding_height,
            )

        # Inferenza conservativa: non si assegna un valore ABL. Se il primo
        # livello già sopra soglia è entro 600 m, il primo crossing definitorio
        # deve trovarsi a quota <= quel livello e quindi certamente sotto la
        # soglia minima richiesta per le classi 1-4.
        if (
            np.isfinite(first_exceeding_height)
            and (
                first_exceeding_height <= thresholds.min_classifiable_abl_m
                or np.isclose(
                    first_exceeding_height,
                    thresholds.min_classifiable_abl_m,
                    rtol=0.0,
                    atol=0.05,
                )
            )
        ):
            return np.nan, 3, 0, float(first_exceeding_height)
        return np.nan, 0, 0, first_exceeding_height

    if not np.isfinite([theta_v_surface, u10, v10]).all():
        return fallback()

    wind2 = (u - u10) ** 2 + (v - v10) ** 2
    wind2 = np.maximum(wind2, 0.5)
    rib = (G / theta_v_surface) * (theta_v - theta_v_surface) * z_agl / wind2
    abl, bracketed, first_exceeding_height = first_bracketed_threshold_crossing(
        z_agl, rib, RIB_CRITICAL, zmin=200.0
    )

    if (
        bracketed
        and np.isfinite(abl)
        and thresholds.min_abl_m <= abl <= thresholds.max_abl_m
    ):
        return float(abl), 1, 1, first_exceeding_height

    return fallback(first_exceeding_height)


def original_layer_gradient(
    z: np.ndarray,
    x: np.ndarray,
    zmin: float,
    zmax: float,
    *,
    min_levels: int,
    min_span_m: float,
) -> tuple[float, int, bool]:
    """Gradiente lineare sui soli livelli originali e relativo controllo qualità."""
    mask = np.isfinite(z) & np.isfinite(x) & (z >= zmin) & (z <= zmax)
    count = int(mask.sum())
    if count < min_levels or zmax <= zmin:
        return np.nan, count, False
    zz = np.asarray(z[mask], dtype=float)
    xx = np.asarray(x[mask], dtype=float)
    span = float(np.nanmax(zz) - np.nanmin(zz))
    if not np.isfinite(span) or span < min_span_m:
        return np.nan, count, False
    return float(np.polyfit(zz, xx, deg=1)[0]), count, True


def rh_top_is_resolved(z: np.ndarray, abl_m: float, max_gap_m: float) -> bool:
    """Verifica che il top ABL sia racchiuso da livelli originali non troppo distanti."""
    zz = np.asarray(z, dtype=float)
    zz = zz[np.isfinite(zz)]
    below = zz[zz <= abl_m]
    above = zz[zz >= abl_m]
    if len(below) == 0 or len(above) == 0:
        return False
    gap = float(np.nanmin(above) - np.nanmax(below))
    return np.isfinite(gap) and 0.0 <= gap <= max_gap_m


def maximum_shear_height_original(
    z: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    *,
    zmin: float = 200.0,
    zmax: float = 6000.0,
    min_levels: int = 5,
    max_layer_gap_m: float = 1000.0,
) -> tuple[float, int, bool]:
    """
    Localizza il massimo shear usando esclusivamente differenze tra livelli originali.

    Il massimo è considerato risolto solo se ci sono almeno ``min_levels`` nel
    layer di ricerca, il massimo non cade sul primo/ultimo intervallo e gli
    intervalli adiacenti al massimo non superano ``max_layer_gap_m``.
    """
    mask = (
        np.isfinite(z) & np.isfinite(u) & np.isfinite(v)
        & (z >= zmin) & (z <= zmax)
    )
    zz, uu, vv = clean_profile(z[mask], u[mask], v[mask])
    nlev = int(len(zz))
    if nlev < min_levels:
        return np.nan, nlev, False

    dz = np.diff(zz)
    valid = np.isfinite(dz) & (dz > 0.0)
    shear = np.full(len(dz), np.nan, dtype=float)
    shear[valid] = np.hypot(np.diff(uu)[valid], np.diff(vv)[valid]) / dz[valid]
    if not np.isfinite(shear).any():
        return np.nan, nlev, False

    i = int(np.nanargmax(shear))
    z_mid = float(0.5 * (zz[i] + zz[i + 1]))
    interior = 0 < i < len(shear) - 1
    local_gaps = dz[max(0, i - 1): min(len(dz), i + 2)]
    gaps_ok = len(local_gaps) > 0 and np.isfinite(local_gaps).all() and float(np.max(local_gaps)) <= max_layer_gap_m
    return z_mid, nlev, bool(interior and gaps_ok)


# -----------------------------------------------------------------------------
# Classificazione adattiva
# -----------------------------------------------------------------------------

def continuous_score(
    lcl_ratio: float,
    theta_grad_abl: float,
    stability_overshoot_capable: bool,
    stability_deep_capable: bool,
    gamma_theta_ft: float,
    shear_distance_ratio: float,
    rh_top_abl: float,
    th: Thresholds,
) -> float:
    """Score 0-100 ripesato sui diagnostici realmente disponibili."""
    if not np.isfinite(lcl_ratio):
        return np.nan

    weighted_terms: list[tuple[float, float]] = []
    lcl_term = float(np.clip(1.0 - abs(lcl_ratio - 1.0), 0.0, 1.0))
    weighted_terms.append((30.0, lcl_term))

    if np.isfinite(theta_grad_abl):
        stability_term = float(np.clip(
            (th.theta_grad_overshoot_max - theta_grad_abl)
            / (th.theta_grad_overshoot_max + 1.0e-3),
            0.0,
            1.0,
        ))
    else:
        stability_term = 1.0 if stability_deep_capable else (0.5 if stability_overshoot_capable else 0.0)
    weighted_terms.append((25.0, stability_term))

    if np.isfinite(gamma_theta_ft):
        cap_term = float(np.clip(
            (th.gamma_theta_strong_cap - gamma_theta_ft)
            / (th.gamma_theta_strong_cap - 2.5e-3),
            0.0,
            1.0,
        ))
        weighted_terms.append((20.0, cap_term))

    if np.isfinite(shear_distance_ratio):
        weighted_terms.append((15.0, float(np.clip(1.0 - shear_distance_ratio, 0.0, 1.0))))

    if np.isfinite(rh_top_abl):
        weighted_terms.append((10.0, float(np.clip((rh_top_abl - 55.0) / 40.0, 0.0, 1.0))))

    total_weight = sum(weight for weight, _ in weighted_terms)
    if total_weight <= 0:
        return np.nan
    return float(np.clip(100.0 * sum(weight * term for weight, term in weighted_terms) / total_weight, 0.0, 100.0))


def classify_pyroconvection_adaptive(
    *,
    abl_m: float,
    lcl_m: float,
    theta_grad_abl: float,
    stability_overshoot_capable: bool,
    stability_deep_capable: bool,
    gamma_theta_ft: float,
    cap_resolved: bool,
    rh_top_abl: float,
    rh_resolved: bool,
    z_shear_m: float,
    shear_resolved: bool,
    classification_mode: int,
    th: Thresholds,
) -> tuple[int, float, float, float, int]:
    """Ladder adattiva con motivazione esplicita della classe finale."""
    if not np.isfinite([abl_m, lcl_m]).all() or abl_m <= 0:
        return NODATA_CLASS, np.nan, np.nan, np.nan, 0

    lcl_ratio = float(lcl_m / abl_m)
    shear_distance_ratio = np.nan
    if shear_resolved and np.isfinite(z_shear_m):
        shear_distance_ratio = float(
            min(abs(z_shear_m - abl_m), abs(z_shear_m - lcl_m)) / abl_m
        )

    score = continuous_score(
        lcl_ratio=lcl_ratio,
        theta_grad_abl=theta_grad_abl,
        stability_overshoot_capable=stability_overshoot_capable,
        stability_deep_capable=stability_deep_capable,
        gamma_theta_ft=gamma_theta_ft,
        shear_distance_ratio=shear_distance_ratio,
        rh_top_abl=rh_top_abl,
        th=th,
    )

    # Filtro operativo coerente con il prodotto di riferimento: sotto 600 m
    # una ABL notturna/poco sviluppata resta Surface plume.
    if abl_m < th.min_classifiable_abl_m:
        return 0, score, lcl_ratio, shear_distance_ratio, 1

    if not stability_overshoot_capable:
        return 1, score, lcl_ratio, shear_distance_ratio, 2

    if lcl_ratio > th.lcl_ratio_overshoot_max:
        return 1, score, lcl_ratio, shear_distance_ratio, 3

    if 1.0 <= lcl_ratio <= th.lcl_ratio_overshoot_max:
        return 2, score, lcl_ratio, shear_distance_ratio, 4

    if not stability_deep_capable:
        return 1, score, lcl_ratio, shear_distance_ratio, 5

    if not np.isfinite(rh_top_abl) or rh_top_abl < th.rh_top_moist:
        return 1, score, lcl_ratio, shear_distance_ratio, 6

    if not rh_resolved:
        return 3, score, lcl_ratio, shear_distance_ratio, 7

    if not cap_resolved or not np.isfinite(gamma_theta_ft):
        return 3, score, lcl_ratio, shear_distance_ratio, 8

    if gamma_theta_ft > th.gamma_theta_weak_cap:
        return 3, score, lcl_ratio, shear_distance_ratio, 9

    shear_close = bool(
        shear_resolved
        and np.isfinite(shear_distance_ratio)
        and shear_distance_ratio <= th.shear_distance_deep
    )

    if classification_mode == 3:
        if not shear_close:
            return 3, score, lcl_ratio, shear_distance_ratio, 10
        return 4, score, lcl_ratio, shear_distance_ratio, 12

    return 4, score, lcl_ratio, shear_distance_ratio, 11


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

    scalar_values = [surface_geopotential, surface_pressure_pa, t2m_k, td2m_k, u10, v10]
    if not np.isfinite(scalar_values).all() or surface_pressure_pa <= 0:
        return metrics

    p_pa = np.asarray(pressure_hpa, dtype=float) * 100.0
    z_asl = np.asarray(geopotential, dtype=float) / G
    orog_m = surface_geopotential / G
    z_agl = z_asl - orog_m

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

    z_agl, p_pa, temperature_k, q, u, v = clean_profile(z_agl, p_pa, temperature_k, q, u, v)
    metrics.n_levels_valid = int(len(z_agl))
    if len(z_agl) < 4:
        return metrics

    theta = potential_temperature(temperature_k, p_pa)
    theta_v = virtual_potential_temperature(theta, q)

    q_sfc = specific_humidity_from_dewpoint(td2m_k, surface_pressure_pa)
    theta_sfc = float(potential_temperature(t2m_k, surface_pressure_pa))
    theta_v_sfc = theta_sfc * (1.0 + 0.61 * q_sfc)

    abl_m, abl_source, rib_crossing_resolved, abl_upper_bound_m = estimate_abl_rib(
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
    metrics.abl_source = abl_source
    metrics.abl_upper_bound_m = abl_upper_bound_m
    metrics.rib_crossing_resolved = rib_crossing_resolved

    # Caso speciale e rigorosamente limitato: il primo livello originale sopra
    # 200 m supera già RiB=0.33, PBLH non è valida e quel livello è <= 600 m.
    # Non si inventa l'altezza ABL e non si eseguono le diagnostiche che la
    # richiedono; si usa soltanto il limite superiore per assegnare Surface plume.
    if abl_source == 3:
        metrics.pyro_class = 0
        metrics.classification_mode = 4
        metrics.classification_quality = 1
        metrics.classification_reason = 13
        return metrics

    if not np.isfinite(abl_m):
        return metrics

    lcl_m = lcl_height_bolton_m(t2m_k, td2m_k, surface_pressure_pa)
    if not np.isfinite(lcl_m):
        return metrics

    z_top = min(thresholds.max_profile_m, float(np.nanmax(z_agl)))
    z_grid = np.arange(50.0, z_top + 25.0, 50.0)
    if len(z_grid) < 10:
        return metrics

    theta_grid = interp_to_grid(z_agl, theta, z_grid)
    t_grid = interp_to_grid(z_agl, temperature_k, z_grid)
    q_grid = interp_to_grid(z_agl, q, z_grid)
    p_grid = interp_to_grid(z_agl, p_pa, z_grid)
    rh_grid = relative_humidity_from_q_t_p(q_grid, t_grid, p_grid)

    # Stabilità ABL: gradiente diretto solo se il layer è realmente campionato.
    abl_low = max(100.0, 0.20 * abl_m)
    abl_high = max(abl_low + 150.0, 0.80 * abl_m)
    abl_high = min(abl_high, abl_m - 50.0 if abl_m > 250.0 else abl_m)
    theta_grad_abl, n_abl, stability_resolved = original_layer_gradient(
        z_agl,
        theta,
        abl_low,
        abl_high,
        min_levels=thresholds.min_levels_abl_gradient,
        min_span_m=thresholds.min_span_abl_gradient_m,
    )
    metrics.n_levels_abl = n_abl
    if stability_resolved:
        stability_overshoot_capable = bool(theta_grad_abl <= thresholds.theta_grad_overshoot_max)
        stability_deep_capable = bool(theta_grad_abl <= thresholds.theta_grad_stable)
        metrics.stability_method = 2
    else:
        # Proxy validato per colonne scarsamente campionate: profondità di mixing.
        theta_grad_abl = np.nan
        proxy_capable = bool(abl_m >= thresholds.mixing_depth_proxy_m)
        stability_overshoot_capable = proxy_capable
        stability_deep_capable = proxy_capable
        metrics.stability_method = 1

    # Cap sopra ABL: classificazione deep solo se almeno due livelli originali
    # coprono un tratto verticale sufficiente. Il valore interpolato resta una
    # diagnostica, ma non abilita classe 4.
    ft_low = abl_m + 200.0
    ft_high = min(abl_m + 1200.0, z_top)
    gamma_direct, n_cap, cap_resolved = original_layer_gradient(
        z_agl,
        theta,
        ft_low,
        ft_high,
        min_levels=thresholds.min_levels_cap_gradient,
        min_span_m=thresholds.min_span_cap_gradient_m,
    )
    metrics.n_levels_cap = n_cap
    gamma_theta_ft = gamma_direct if cap_resolved else linear_gradient(z_grid, theta_grid, ft_low, ft_high)

    # RH media al top ABL; la risoluzione è valutata sui livelli originali che
    # racchiudono il top, non sul numero di punti della griglia interpolata.
    rh_top_abl = mean_in_layer(z_grid, rh_grid, max(50.0, abl_m - 150.0), abl_m + 150.0)
    rh_resolved = rh_top_is_resolved(z_agl, abl_m, thresholds.max_rh_bracket_gap_m)

    # Shear su livelli originali. Nessun massimo viene inventato a 50 m.
    z_shear_m, n_shear, shear_resolved = maximum_shear_height_original(
        z_agl,
        u,
        v,
        min_levels=thresholds.min_levels_shear,
        max_layer_gap_m=thresholds.max_shear_layer_gap_m,
    )
    metrics.n_levels_shear = n_shear

    if stability_resolved and cap_resolved and rh_resolved and shear_resolved:
        classification_mode = 3
    elif stability_resolved:
        classification_mode = 2
    else:
        classification_mode = 1

    classification_quality = classification_mode
    if abl_source == 2:
        classification_quality = max(1, classification_quality - 1)
    if abl_m < thresholds.low_confidence_abl_m:
        classification_quality = max(1, classification_quality - 1)
    if not cap_resolved or not rh_resolved:
        classification_quality = max(1, classification_quality - 1)

    pyro_class, score, lcl_ratio, shear_distance_ratio, classification_reason = classify_pyroconvection_adaptive(
        abl_m=abl_m,
        lcl_m=lcl_m,
        theta_grad_abl=theta_grad_abl,
        stability_overshoot_capable=stability_overshoot_capable,
        stability_deep_capable=stability_deep_capable,
        gamma_theta_ft=gamma_theta_ft,
        cap_resolved=cap_resolved,
        rh_top_abl=rh_top_abl,
        rh_resolved=rh_resolved,
        z_shear_m=z_shear_m,
        shear_resolved=shear_resolved,
        classification_mode=classification_mode,
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
    metrics.classification_mode = classification_mode
    metrics.classification_quality = classification_quality
    metrics.classification_reason = classification_reason
    metrics.cap_resolved = int(cap_resolved)
    metrics.rh_resolved = int(rh_resolved)
    metrics.shear_resolved = int(shear_resolved)
    return metrics


# -----------------------------------------------------------------------------
# Elaborazione dominio
# -----------------------------------------------------------------------------

def _drop_extra_singleton_dims(
    da: xr.DataArray,
    required_dims: Sequence[str],
    label: str,
) -> xr.DataArray:
    """Elimina dimensioni extra solo se singleton, altrimenti segnala l'errore."""
    out = da
    extras = [dim for dim in out.dims if dim not in required_dims]
    for dim in extras:
        if out.sizes[dim] != 1:
            raise ValueError(
                f"La variabile {label!r} contiene la dimensione extra {dim!r} "
                f"di lunghezza {out.sizes[dim]}; dimensioni: {out.dims}."
            )
        out = out.isel({dim: 0}, drop=True)
    return out


def _dimension_coordinate_values(ds: xr.Dataset, dim: str) -> np.ndarray:
    if dim in ds.coords and ds[dim].ndim == 1 and ds[dim].dims == (dim,):
        return np.asarray(ds[dim].values)
    return np.arange(ds.sizes[dim], dtype=np.int32)



def ensure_provinces_geojson(
    *,
    explicit_path: Path | None = None,
    cache_dir: Path = DEFAULT_BOUNDARIES_CACHE_DIR,
    url: str = OPENPOLIS_PROVINCES_URL,
    required: bool = True,
) -> Path | None:
    """Restituisce il GeoJSON dei confini provinciali italiani.

    Se non viene fornito un file locale, scarica una copia nella cache. Lo stesso
    layer viene usato sia per tracciare le province sia, in assenza di LSM nel
    WRF, per derivare una land-sea mask geometrica sul territorio italiano.
    """
    if explicit_path is not None:
        path = Path(explicit_path)
        if not path.exists():
            raise FileNotFoundError(f"File dei confini provinciali non trovato: {path}")
        return path

    cache_dir = Path(cache_dir)
    path = cache_dir / "limits_IT_provinces.geojson"
    if path.exists() and path.stat().st_size > 1024:
        return path

    try:
        cache_dir.mkdir(parents=True, exist_ok=True)
        logging.info("Scarico i confini provinciali ISTAT/Openpolis: %s", url)
        response = requests.get(url, timeout=120)
        response.raise_for_status()
        tmp = path.with_suffix(path.suffix + ".part")
        tmp.write_bytes(response.content)
        if tmp.stat().st_size <= 1024:
            raise RuntimeError("GeoJSON scaricato vuoto o incompleto")
        tmp.replace(path)
        logging.info("Confini provinciali salvati in: %s", path)
        return path
    except Exception as exc:
        if required:
            raise RuntimeError(
                "Impossibile ottenere i confini provinciali necessari per la maschera del mare. "
                "Specificare --provinces-file con un GeoJSON/Shapefile locale oppure rendere "
                "disponibile la connessione al primo avvio."
            ) from exc
        logging.warning("Confini provinciali non disponibili: %s", exc)
        return None


@lru_cache(maxsize=8)
def _load_geojson_geometries_cached(path_string: str) -> tuple[object, ...]:
    """Carica e memorizza geometrie GeoJSON/Shapefile come oggetti Shapely."""
    path = Path(path_string)
    try:
        from shapely.geometry import shape
    except Exception as exc:
        raise RuntimeError(
            "Per usare la maschera geografica e i confini provinciali è richiesto shapely. "
            "Installare con: pip install shapely"
        ) from exc

    geometries: list[object] = []
    if path.suffix.lower() in {".json", ".geojson"}:
        payload = json.loads(path.read_text(encoding="utf-8"))
        for feature in payload.get("features", []):
            geometry = feature.get("geometry")
            if geometry:
                geom = shape(geometry)
                if not geom.is_empty:
                    geometries.append(geom)
    else:
        try:
            import fiona
            with fiona.open(path) as source:
                for feature in source:
                    geometry = feature.get("geometry")
                    if geometry:
                        geom = shape(geometry)
                        if not geom.is_empty:
                            geometries.append(geom)
        except Exception as exc:
            raise RuntimeError(f"Impossibile leggere il file vettoriale {path}: {exc}") from exc

    if not geometries:
        raise RuntimeError(f"Nessuna geometria valida trovata in {path}")
    return tuple(geometries)


def load_boundary_geometries(path: Path) -> tuple[object, ...]:
    return _load_geojson_geometries_cached(str(Path(path).resolve()))


def derive_land_mask_from_boundaries(
    latitude: np.ndarray,
    longitude: np.ndarray,
    boundaries_path: Path,
) -> np.ndarray:
    """Crea una maschera terra/mare sui centri cella usando le province italiane."""
    try:
        import shapely
    except Exception as exc:
        raise RuntimeError(
            "La land-sea mask geometrica richiede shapely. Installare con: pip install shapely"
        ) from exc

    lat = np.asarray(latitude, dtype=float)
    lon = np.asarray(longitude, dtype=float)
    if lat.ndim == 1 and lon.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon, lat)
    elif lat.ndim == 2 and lon.ndim == 2 and lat.shape == lon.shape:
        lat2d, lon2d = lat, lon
    else:
        raise ValueError(
            "Coordinate non compatibili con la maschera geografica: "
            f"latitude={lat.shape}, longitude={lon.shape}"
        )

    geoms = load_boundary_geometries(boundaries_path)
    merged = shapely.union_all(np.asarray(geoms, dtype=object))
    valid = np.isfinite(lat2d) & np.isfinite(lon2d)
    mask = np.zeros(lat2d.shape, dtype=np.float32)
    if valid.any():
        # intersects_xy include anche i punti esattamente sul bordo costiero.
        mask[valid] = shapely.intersects_xy(
            merged,
            lon2d[valid],
            lat2d[valid],
        ).astype(np.float32)
    return mask


def normalize_land_sea_mask_values(da: xr.DataArray) -> xr.DataArray:
    """Normalizza LSM espresso come 0-1, 0-100 o maschera binaria."""
    finite = da.where(np.isfinite(da))
    try:
        vmax = float(finite.max(skipna=True).values)
    except Exception:
        vmax = np.nan
    if np.isfinite(vmax) and vmax > 1.5:
        da = da / 100.0
    return xr.where(np.isfinite(da), da.clip(min=0.0, max=1.0), 0.0)


def draw_shapely_boundaries_plain(ax, geometries, *, linewidth: float = 0.45, color: str = "0.15", zorder: int = 8) -> None:
    """Disegna Polygon/MultiPolygon su assi matplotlib senza Cartopy."""
    def draw_polygon(poly) -> None:
        x, y = poly.exterior.xy
        ax.plot(x, y, color=color, linewidth=linewidth, zorder=zorder)
        for ring in poly.interiors:
            x, y = ring.xy
            ax.plot(x, y, color=color, linewidth=max(0.2, linewidth * 0.7), zorder=zorder)

    for geom in geometries:
        if geom.geom_type == "Polygon":
            draw_polygon(geom)
        elif geom.geom_type == "MultiPolygon":
            for poly in geom.geoms:
                draw_polygon(poly)


def process_domain(
    ds_sl: xr.Dataset,
    ds_pl: xr.Dataset,
    thresholds: Thresholds,
    mask_sea: bool = True,
    store_diagnostics: bool = True,
    fallback_land_boundaries: Path | None = None,
) -> xr.Dataset:
    # Alias comuni nei prodotti WRF/MeteoHub (NetCDF e GRIB).
    t2m = convert_temperature_to_k(find_var(ds_sl, ["t2m", "T2", "2t", "temperature_2m", "2m_temperature"]))
    d2m = find_optional_var(ds_sl, ["d2m", "TD2", "2d", "dewpoint_2m", "2m_dewpoint_temperature"])
    if d2m is not None:
        d2m = convert_temperature_to_k(d2m)
    else:
        rh2 = find_var(ds_sl, ["rh2", "RH2", "2r", "relative_humidity_2m", "2m_relative_humidity"])
        d2m = dewpoint_from_t_rh(t2m, rh2)
    u10_da = find_var(ds_sl, ["u10", "U10", "10u", "u_component_of_wind_10m", "10m_u_component_of_wind"])
    v10_da = find_var(ds_sl, ["v10", "V10", "10v", "v_component_of_wind_10m", "10m_v_component_of_wind"])
    blh = find_var(ds_sl, ["blh", "PBLH", "hpbl", "boundary_layer_height"])
    z_sfc = height_or_geopotential_to_geopotential(find_var(ds_sl, ["z_sfc", "orog", "HGT", "hgt", "terrain_height", "surface_geopotential", "z"]))
    sp_raw = find_optional_var(ds_sl, ["sp", "PSFC", "surface_pressure", "pres_sfc"])
    if sp_raw is not None:
        sp = convert_pressure_to_pa(sp_raw)
    else:
        msl = convert_pressure_to_pa(find_var(ds_sl, ["prmsl", "msl", "MSLP", "mean_sea_level_pressure"]))
        terrain_m = z_sfc / G
        # Equazione ipsometrica semplificata. È un fallback operativo quando MeteoHub
        # distribuisce PRMSL ma non la pressione al suolo.
        sp = msl * np.exp(-G * terrain_m / (RD * xr.where(t2m > 180.0, t2m, 288.15)))
        sp.attrs.update({"units": "Pa", "long_name": "surface pressure estimated from MSLP, terrain and T2m"})
        logging.warning("Pressione al suolo non presente: stimata da PRMSL, orografia e T2m")
    lsm = find_optional_var(ds_sl, ["lsm", "LANDMASK", "landmask", "land_sea_mask", "land_fraction"])

    z_pl = height_or_geopotential_to_geopotential(find_var(ds_pl, ["gh", "z", "geopotential_height", "geopotential", "HGT"]))
    t_pl = convert_temperature_to_k(find_var(ds_pl, ["t", "T", "temperature", "air_temperature"]))
    q_pl = find_optional_var(ds_pl, ["q", "QVAPOR", "specific_humidity", "sh"])
    if q_pl is None:
        rh_pl = find_var(ds_pl, ["r", "rh", "RH", "relative_humidity"])
        q_pl = specific_humidity_from_rh_da(rh_pl, t_pl, ds_pl["level"])
    u_pl = find_var(ds_pl, ["u", "U", "u_component_of_wind", "ua"])
    v_pl = find_var(ds_pl, ["v", "V", "v_component_of_wind", "va"])

    sl_y, sl_x, sl_lat, sl_lon, sl_curvilinear = horizontal_grid_info(ds_sl)
    pl_y, pl_x, _, _, _ = horizontal_grid_info(ds_pl)

    pl_order = ("time", "level", pl_y, pl_x)
    sl_order = ("time", sl_y, sl_x)

    def normalize_pl(da: xr.DataArray, label: str) -> xr.DataArray:
        if "time" not in da.dims:
            da = da.expand_dims(time=ds_pl.time)
        da = _drop_extra_singleton_dims(da, pl_order, label)
        missing = [dim for dim in pl_order if dim not in da.dims]
        if missing:
            raise ValueError(
                f"La variabile {label!r} non contiene le dimensioni richieste {missing}; "
                f"dimensioni presenti: {da.dims}."
            )
        return da.transpose(*pl_order)

    z_pl = normalize_pl(z_pl, str(z_pl.name or "geopotential"))
    t_pl = normalize_pl(t_pl, str(t_pl.name or "temperature"))
    q_pl = normalize_pl(q_pl, str(q_pl.name or "specific humidity"))
    u_pl = normalize_pl(u_pl, str(u_pl.name or "u wind"))
    v_pl = normalize_pl(v_pl, str(v_pl.name or "v wind"))

    def normalize_sl(
        da: xr.DataArray,
        label: str,
        reference: xr.DataArray | None = None,
    ) -> xr.DataArray:
        if "time" not in da.dims:
            da = da.expand_dims(time=ds_sl.time)
        da = _drop_extra_singleton_dims(da, sl_order, label)
        if reference is not None:
            da = da.broadcast_like(reference)
        missing = [dim for dim in sl_order if dim not in da.dims]
        if missing:
            raise ValueError(
                f"La variabile {label!r} non contiene le dimensioni richieste {missing}; "
                f"dimensioni presenti: {da.dims}."
            )
        return da.transpose(*sl_order)

    t2m = normalize_sl(t2m, str(t2m.name or "2 m temperature"))
    d2m = normalize_sl(d2m, str(d2m.name or "2 m dewpoint"), t2m)
    u10_da = normalize_sl(u10_da, str(u10_da.name or "10 m u wind"), t2m)
    v10_da = normalize_sl(v10_da, str(v10_da.name or "10 m v wind"), t2m)
    sp = normalize_sl(sp, str(sp.name or "surface pressure"), t2m)
    blh = normalize_sl(blh, str(blh.name or "boundary layer height"), t2m)
    z_sfc = normalize_sl(z_sfc, str(z_sfc.name or "surface geopotential"), t2m)
    if lsm is not None:
        lsm = normalize_sl(lsm, str(lsm.name or "land sea mask"), t2m)
        lsm = normalize_land_sea_mask_values(lsm)
        logging.info("Maschera terra/mare: uso il campo %s presente nei dati WRF", lsm.name or "LSM")
    elif mask_sea:
        if fallback_land_boundaries is None:
            raise RuntimeError(
                "Campo LSM/LANDMASK assente e nessun file di confini disponibile per derivare la maschera del mare"
            )
        geometric_mask = derive_land_mask_from_boundaries(
            latitude=np.asarray(sl_lat),
            longitude=np.asarray(sl_lon),
            boundaries_path=fallback_land_boundaries,
        )
        base_dims = (sl_y, sl_x) if geometric_mask.ndim == 2 else t2m.dims[-2:]
        base_mask = xr.DataArray(geometric_mask, dims=base_dims)
        lsm = base_mask.expand_dims(time=t2m.time).transpose(*sl_order).broadcast_like(t2m)
        lsm.name = "land_sea_mask"
        lsm.attrs.update(
            long_name="land sea mask derived from Italian province polygons",
            units="1",
            source=str(fallback_land_boundaries),
        )
        logging.warning(
            "Campo LSM/LANDMASK assente: land-sea mask derivata geometricamente da %s",
            fallback_land_boundaries,
        )
    else:
        lsm = xr.ones_like(t2m, dtype=np.float32)
        lsm.name = "land_sea_mask"
        logging.info("Opzione --keep-sea: nessuna mascheratura marina applicata")

    times = ds_sl.time.values
    levels = np.asarray(ds_pl.level.values, dtype=float)
    nt = len(times)
    ny = ds_sl.sizes[sl_y]
    nx = ds_sl.sizes[sl_x]
    if ds_pl.sizes[pl_y] != ny or ds_pl.sizes[pl_x] != nx:
        raise ValueError(
            "Dimensioni orizzontali diverse dopo l'allineamento: "
            f"single=({ny}, {nx}), pressure=({ds_pl.sizes[pl_y]}, {ds_pl.sizes[pl_x]})."
        )
    shape = (nt, ny, nx)

    # 255 distingue NoData/non classificabile da una vera classe 0 Surface plume.
    pyro_class = np.full(shape, NODATA_CLASS, dtype=np.uint8)
    classification_mode = np.zeros(shape, dtype=np.uint8)
    classification_quality = np.zeros(shape, dtype=np.uint8)
    classification_reason = np.zeros(shape, dtype=np.uint8)
    if store_diagnostics:
        score = np.full(shape, np.nan, dtype=np.float32)
        abl_m = np.full(shape, np.nan, dtype=np.float32)
        lcl_m = np.full(shape, np.nan, dtype=np.float32)
        lcl_ratio = np.full(shape, np.nan, dtype=np.float32)
        z_shear_m = np.full(shape, np.nan, dtype=np.float32)
        shear_distance_ratio = np.full(shape, np.nan, dtype=np.float32)
        theta_grad_abl = np.full(shape, np.nan, dtype=np.float32)
        gamma_theta_ft = np.full(shape, np.nan, dtype=np.float32)
        rh_top_abl = np.full(shape, np.nan, dtype=np.float32)
        abl_source = np.zeros(shape, dtype=np.uint8)
        abl_upper_bound_m = np.full(shape, np.nan, dtype=np.float32)
        rib_crossing_resolved = np.zeros(shape, dtype=np.uint8)
        stability_method = np.zeros(shape, dtype=np.uint8)
        n_levels_valid = np.zeros(shape, dtype=np.uint8)
        n_levels_abl = np.zeros(shape, dtype=np.uint8)
        n_levels_cap = np.zeros(shape, dtype=np.uint8)
        n_levels_shear = np.zeros(shape, dtype=np.uint8)
        cap_resolved = np.zeros(shape, dtype=np.uint8)
        rh_resolved = np.zeros(shape, dtype=np.uint8)
        shear_resolved = np.zeros(shape, dtype=np.uint8)
    else:
        score = abl_m = lcl_m = lcl_ratio = None
        z_shear_m = shear_distance_ratio = None
        theta_grad_abl = gamma_theta_ft = rh_top_abl = None
        abl_source = abl_upper_bound_m = rib_crossing_resolved = stability_method = None
        n_levels_valid = n_levels_abl = n_levels_cap = n_levels_shear = None
        cap_resolved = rh_resolved = shear_resolved = None
        logging.info(
            "Modalità class-only: salvo classe, modalità, qualità, motivazione e land mask"
        )

    # La versione precedente eseguiva ``.values`` sull'intero intervallo
    # temporale, caricando contemporaneamente in RAM cinque campi 4D ai livelli
    # isobarici. Qui viene decodificato un solo timestamp alla volta.
    lsm_out = np.full(shape, np.nan, dtype=np.float32)
    bytes_per_time = 5 * len(levels) * ny * nx * np.dtype(np.float32).itemsize
    logging.info(
        "Caricamento pressure levels per singolo timestamp: circa %.1f MB "
        "per i cinque campi principali", bytes_per_time / 1024**2
    )

    total = nt * ny * nx
    done = 0
    for it, timestamp in enumerate(times):
        logging.info("Elaboro %s (%d/%d)", pd.Timestamp(timestamp), it + 1, nt)

        # I cast float32 limitano ulteriormente il picco di memoria senza perdita
        # significativa rispetto alla precisione originaria dei GRIB.
        z_pl_np = np.asarray(z_pl.isel(time=it).values, dtype=np.float32)
        t_pl_np = np.asarray(t_pl.isel(time=it).values, dtype=np.float32)
        q_pl_np = np.asarray(q_pl.isel(time=it).values, dtype=np.float32)
        u_pl_np = np.asarray(u_pl.isel(time=it).values, dtype=np.float32)
        v_pl_np = np.asarray(v_pl.isel(time=it).values, dtype=np.float32)

        t2m_np = np.asarray(t2m.isel(time=it).values, dtype=np.float32)
        d2m_np = np.asarray(d2m.isel(time=it).values, dtype=np.float32)
        u10_np = np.asarray(u10_da.isel(time=it).values, dtype=np.float32)
        v10_np = np.asarray(v10_da.isel(time=it).values, dtype=np.float32)
        sp_np = np.asarray(sp.isel(time=it).values, dtype=np.float32)
        blh_np = np.asarray(blh.isel(time=it).values, dtype=np.float32)
        z_sfc_np = np.asarray(z_sfc.isel(time=it).values, dtype=np.float32)
        lsm_np = np.asarray(lsm.isel(time=it).values, dtype=np.float32)
        lsm_out[it] = lsm_np

        for iy in range(ny):
            for ix in range(nx):
                done += 1
                if mask_sea and lsm_np[iy, ix] < thresholds.min_land_fraction:
                    continue

                m = process_profile(
                    pressure_hpa=levels,
                    geopotential=z_pl_np[:, iy, ix],
                    temperature_k=t_pl_np[:, iy, ix],
                    q=q_pl_np[:, iy, ix],
                    u=u_pl_np[:, iy, ix],
                    v=v_pl_np[:, iy, ix],
                    surface_geopotential=float(z_sfc_np[iy, ix]),
                    surface_pressure_pa=float(sp_np[iy, ix]),
                    t2m_k=float(t2m_np[iy, ix]),
                    td2m_k=float(d2m_np[iy, ix]),
                    u10=float(u10_np[iy, ix]),
                    v10=float(v10_np[iy, ix]),
                    blh_fallback_m=float(blh_np[iy, ix]),
                    thresholds=thresholds,
                )

                pyro_class[it, iy, ix] = np.uint8(m.pyro_class)
                classification_mode[it, iy, ix] = np.uint8(m.classification_mode)
                classification_quality[it, iy, ix] = np.uint8(m.classification_quality)
                classification_reason[it, iy, ix] = np.uint8(m.classification_reason)
                if store_diagnostics:
                    assert score is not None and abl_m is not None and lcl_m is not None
                    assert lcl_ratio is not None and z_shear_m is not None
                    assert shear_distance_ratio is not None and theta_grad_abl is not None
                    assert gamma_theta_ft is not None and rh_top_abl is not None
                    score[it, iy, ix] = m.score
                    abl_m[it, iy, ix] = m.abl_m
                    lcl_m[it, iy, ix] = m.lcl_m
                    lcl_ratio[it, iy, ix] = m.lcl_ratio
                    z_shear_m[it, iy, ix] = m.z_shear_m
                    shear_distance_ratio[it, iy, ix] = m.shear_distance_ratio
                    theta_grad_abl[it, iy, ix] = m.theta_grad_abl
                    gamma_theta_ft[it, iy, ix] = m.gamma_theta_ft
                    rh_top_abl[it, iy, ix] = m.rh_top_abl
                    assert abl_source is not None and abl_upper_bound_m is not None and rib_crossing_resolved is not None and stability_method is not None
                    assert n_levels_valid is not None and n_levels_abl is not None
                    assert n_levels_cap is not None and n_levels_shear is not None
                    assert cap_resolved is not None and rh_resolved is not None and shear_resolved is not None
                    abl_source[it, iy, ix] = np.uint8(m.abl_source)
                    abl_upper_bound_m[it, iy, ix] = m.abl_upper_bound_m
                    rib_crossing_resolved[it, iy, ix] = np.uint8(m.rib_crossing_resolved)
                    stability_method[it, iy, ix] = np.uint8(m.stability_method)
                    n_levels_valid[it, iy, ix] = np.uint8(min(m.n_levels_valid, 255))
                    n_levels_abl[it, iy, ix] = np.uint8(min(m.n_levels_abl, 255))
                    n_levels_cap[it, iy, ix] = np.uint8(min(m.n_levels_cap, 255))
                    n_levels_shear[it, iy, ix] = np.uint8(min(m.n_levels_shear, 255))
                    cap_resolved[it, iy, ix] = np.uint8(m.cap_resolved)
                    rh_resolved[it, iy, ix] = np.uint8(m.rh_resolved)
                    shear_resolved[it, iy, ix] = np.uint8(m.shear_resolved)

        # Libera esplicitamente i grandi array prima di decodificare il tempo
        # successivo, utile soprattutto su Windows.
        del (
            z_pl_np, t_pl_np, q_pl_np, u_pl_np, v_pl_np,
            t2m_np, d2m_np, u10_np, v10_np, sp_np, blh_np,
            z_sfc_np, lsm_np,
        )
        gc.collect()
        logging.debug("Avanzamento %.1f%%", 100.0 * done / total)

    valid_cells = pyro_class != NODATA_CLASS
    n_valid = int(valid_cells.sum())
    if n_valid:
        mode_counts = {
            mode: int(((classification_mode == mode) & valid_cells).sum())
            for mode in (1, 2, 3, 4)
        }
        quality_counts = {quality: int(((classification_quality == quality) & valid_cells).sum()) for quality in (1, 2, 3)}
        shallow_count = int(((classification_reason == 1) & valid_cells).sum())
        shallow_bound_count = int(((classification_reason == 13) & valid_cells).sum())
        stable_count = int(((classification_reason == 2) & valid_cells).sum())
        logging.info(
            "Classificazione adattiva: validi=%d; proxy-noshear=%d, resolved-noshear=%d, full=%d, shallow-upper-bound=%d; qualità bassa=%d, media=%d, alta=%d",
            n_valid,
            mode_counts[1], mode_counts[2], mode_counts[3], mode_counts[4],
            quality_counts[1], quality_counts[2], quality_counts[3],
        )
        logging.info(
            "Motivi principali: ABL calcolata < %.0f m=%d; ABL certamente < %.0f m da limite RiB=%d; mixed layer stabile/poco rimescolato=%d",
            thresholds.min_classifiable_abl_m, shallow_count,
            thresholds.min_classifiable_abl_m, shallow_bound_count,
            stable_count,
        )
    else:
        logging.warning("Nessuna cella classificabile: tutte le celle sono NoData")

    coords: dict[str, object] = {"time": times}
    if sl_curvilinear:
        coords[sl_y] = _dimension_coordinate_values(ds_sl, sl_y)
        coords[sl_x] = _dimension_coordinate_values(ds_sl, sl_x)
        coords["latitude"] = ((sl_y, sl_x), sl_lat)
        coords["longitude"] = ((sl_y, sl_x), sl_lon)
    else:
        if sl_y == "latitude":
            coords["latitude"] = sl_lat
        else:
            coords[sl_y] = _dimension_coordinate_values(ds_sl, sl_y)
            coords["latitude"] = ((sl_y,), sl_lat)
        if sl_x == "longitude":
            coords["longitude"] = sl_lon
        else:
            coords[sl_x] = _dimension_coordinate_values(ds_sl, sl_x)
            coords["longitude"] = ((sl_x,), sl_lon)

    dims = ("time", sl_y, sl_x)
    data_vars: dict[str, object] = {
        "pyro_class": (dims, pyro_class),
        "classification_mode": (dims, classification_mode),
        "classification_quality": (dims, classification_quality),
        "classification_reason": (dims, classification_reason),
        "land_sea_mask": (dims, lsm_out),
    }
    if store_diagnostics:
        assert score is not None and abl_m is not None and lcl_m is not None
        assert lcl_ratio is not None and z_shear_m is not None
        assert shear_distance_ratio is not None and theta_grad_abl is not None
        assert gamma_theta_ft is not None and rh_top_abl is not None
        data_vars.update({
            "pyro_score": (dims, score),
            "abl_height": (dims, abl_m),
            "lcl_height": (dims, lcl_m),
            "lcl_abl_ratio": (dims, lcl_ratio),
            "maximum_shear_height": (dims, z_shear_m),
            "shear_distance_ratio": (dims, shear_distance_ratio),
            "theta_gradient_abl": (dims, theta_grad_abl),
            "gamma_theta_free_troposphere": (dims, gamma_theta_ft),
            "relative_humidity_abl_top": (dims, rh_top_abl),
            "abl_source": (dims, abl_source),
            "abl_upper_bound": (dims, abl_upper_bound_m),
            "rib_crossing_resolved": (dims, rib_crossing_resolved),
            "stability_method": (dims, stability_method),
            "n_original_levels_valid": (dims, n_levels_valid),
            "n_original_levels_abl": (dims, n_levels_abl),
            "n_original_levels_cap": (dims, n_levels_cap),
            "n_original_levels_shear": (dims, n_levels_shear),
            "cap_resolved": (dims, cap_resolved),
            "rh_top_resolved": (dims, rh_resolved),
            "shear_resolved": (dims, shear_resolved),
        })

    out = xr.Dataset(
        data_vars=data_vars,
        coords=coords,
        attrs={
            "title": "WRF 1.5 km adaptive forecast atmospheric predisposition to pyroconvection",
            "method": "Adaptive ladder using original-level sampling; RiB ABL accepted only when vertically bracketed, otherwise WRF PBLH fallback; if PBLH is invalid and the first RiB exceedance is at or below 600 m, only a conservative shallow-ABL upper bound is used to assign Surface plume",
            "software_version": "pyroconvection_wrf15_v1_7_3",
            "reference": "Castellnou et al. (2022), JGR Atmospheres, doi:10.1029/2022JD036920",
            "warning": "WRF is a deterministic forecast. Classes indicate atmospheric predisposition under an intense-fire forcing, not occurrence probability.",
            "horizontal_grid": "curvilinear" if sl_curvilinear else "rectilinear",
            "sea_mask_applied": "true" if mask_sea else "false",
            "land_sea_mask_source": str(lsm.attrs.get("source", lsm.name or "WRF LSM")),
            "class_0": CLASS_NAMES[0],
            "class_1": CLASS_NAMES[1],
            "class_2": CLASS_NAMES[2],
            "class_3": CLASS_NAMES[3],
            "class_4": CLASS_NAMES[4],
        },
    )

    out["latitude"].attrs.update(ds_sl["latitude"].attrs)
    out["longitude"].attrs.update(ds_sl["longitude"].attrs)
    out["pyro_class"].attrs.update(
        long_name="pyroconvection atmospheric class",
        flag_values=np.array([0, 1, 2, 3, 4], dtype=np.uint8),
        flag_meanings="surface_plume convective_plume overshooting_pyrocu resilient_pyrocu pyrocb",
        comment="255 is NoData/non-classifiable and is distinct from class 0 Surface plume",
    )
    out["classification_mode"].attrs.update(
        long_name="adaptive classification ladder",
        flag_values=np.array([0, 1, 2, 3, 4], dtype=np.uint8),
        flag_meanings="invalid proxy_noshear resolved_noshear full shallow_abl_upper_bound",
    )
    out["classification_quality"].attrs.update(
        long_name="classification quality",
        flag_values=np.array([0, 1, 2, 3], dtype=np.uint8),
        flag_meanings="invalid low medium high",
        comment="Quality is reduced by PBLH fallback or unresolved cap/RH-top",
    )
    out["classification_reason"].attrs.update(
        long_name="dominant reason for the final pyroconvection class",
        flag_values=np.array(sorted(CLASSIFICATION_REASON_NAMES), dtype=np.uint8),
        flag_meanings=" ".join(
            CLASSIFICATION_REASON_NAMES[key]
            for key in sorted(CLASSIFICATION_REASON_NAMES)
        ),
        comment="Audit field for the deterministic classification ladder",
    )
    if store_diagnostics:
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
        out["abl_source"].attrs.update(
            long_name="source of ABL height",
            flag_values=np.array([0, 1, 2, 3], dtype=np.uint8),
            flag_meanings="invalid bulk_richardson wrf_pblh_fallback rib_shallow_upper_bound",
            comment="Code 3 does not provide an exact ABL height; it only proves that the first RiB crossing is at or below abl_upper_bound",
        )
        out["abl_upper_bound"].attrs.update(
            long_name="upper bound on ABL height from the first unresolved RiB exceedance",
            units="m",
            comment="Used only as a conservative constraint; it is not an estimated ABL height",
        )
        out["rib_crossing_resolved"].attrs.update(
            long_name="vertical bracketing quality of the Richardson ABL crossing",
            flag_values=np.array([0, 1], dtype=np.uint8),
            flag_meanings="not_resolved resolved",
            comment="A value of 0 with abl_source=2 or 3 includes the case where the first pressure level above 200 m already exceeds RiB 0.33",
        )
        out["stability_method"].attrs.update(
            long_name="ABL stability diagnostic method",
            flag_values=np.array([0, 1, 2], dtype=np.uint8),
            flag_meanings="invalid mixing_depth_proxy original_level_gradient",
        )
        for name in ["n_original_levels_valid", "n_original_levels_abl", "n_original_levels_cap", "n_original_levels_shear"]:
            out[name].attrs["long_name"] = name.replace("_", " ")
        for name in ["cap_resolved", "rh_top_resolved", "shear_resolved"]:
            out[name].attrs.update(
                flag_values=np.array([0, 1], dtype=np.uint8),
                flag_meanings="not_resolved resolved",
            )

    return out


# -----------------------------------------------------------------------------
# Esportazione
# -----------------------------------------------------------------------------

def save_netcdf(ds: xr.Dataset, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    encoding: dict[str, dict] = {}
    for name, da in ds.data_vars.items():
        if np.issubdtype(da.dtype, np.integer):
            encoding[name] = {"zlib": True, "complevel": 4, "dtype": "u1"}
            if name == "pyro_class":
                encoding[name]["_FillValue"] = NODATA_CLASS
        else:
            encoding[name] = {"zlib": True, "complevel": 4, "dtype": "f4"}
    def write_netcdf3_fallback(source: xr.Dataset) -> None:
        fallback = source.copy(deep=False)
        fallback_encoding: dict[str, dict] = {}
        for name, da in source.data_vars.items():
            if np.issubdtype(da.dtype, np.integer):
                if name == "pyro_class":
                    signed = da.astype(np.int16)
                    fallback[name] = xr.where(da == NODATA_CLASS, np.int16(-32768), signed)
                    fallback_encoding[name] = {"dtype": "i2", "_FillValue": -32768}
                else:
                    fallback[name] = da.astype(np.int16)
                    fallback_encoding[name] = {"dtype": "i2"}
            else:
                fallback_encoding[name] = {"dtype": "f4"}
        fallback.to_netcdf(path, engine="scipy", encoding=fallback_encoding)

    engines = xr.backends.list_engines()
    modern_engine = "netcdf4" if "netcdf4" in engines else ("h5netcdf" if "h5netcdf" in engines else None)
    if modern_engine is None:
        logging.warning(
            "Backend NetCDF4/HDF5 non disponibile: uso NetCDF3 con interi signed a 16 bit."
        )
        write_netcdf3_fallback(ds)
    else:
        try:
            ds.to_netcdf(path, engine=modern_engine, encoding=encoding)
        except (ValueError, TypeError, ModuleNotFoundError) as exc:
            logging.warning(
                "Encoding %s non supportato (%s): uso il fallback NetCDF3.",
                modern_engine,
                exc,
            )
            write_netcdf3_fallback(ds)
    logging.info("Salvato NetCDF: %s", path)


def prepare_raster_array(
    da: xr.DataArray,
    nodata: int = NODATA_CLASS,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, bool]:
    """Prepara matrice e coordinate, supportando griglie 1D e curvilinee 2D."""
    arr = np.asarray(da.values)
    if arr.ndim != 2:
        raise ValueError(f"Il raster deve essere 2D; forma trovata: {arr.shape}")

    _, _, lat, lon, curvilinear = horizontal_grid_info(da)
    lat = np.asarray(lat)
    lon = np.asarray(lon)

    if curvilinear:
        lat_first = _finite_edge_mean(lat, axis=0, index=0)
        lat_last = _finite_edge_mean(lat, axis=0, index=-1)
        if np.isfinite(lat_first) and np.isfinite(lat_last) and lat_first < lat_last:
            lat = lat[::-1, :]
            lon = lon[::-1, :]
            arr = arr[::-1, :]

        lon_first = _finite_edge_mean(lon, axis=1, index=0)
        lon_last = _finite_edge_mean(lon, axis=1, index=-1)
        if np.isfinite(lon_first) and np.isfinite(lon_last) and lon_first > lon_last:
            lat = lat[:, ::-1]
            lon = lon[:, ::-1]
            arr = arr[:, ::-1]
    else:
        if lat[0] < lat[-1]:
            lat = lat[::-1]
            arr = arr[::-1, :]
        if lon[0] > lon[-1]:
            lon = lon[::-1]
            arr = arr[:, ::-1]

    arr = np.where(np.isfinite(arr), arr, nodata).astype(np.uint8)
    return arr, lat, lon, curvilinear


def _sample_indices(size: int, maximum: int = 25) -> list[int]:
    if size <= maximum:
        return list(range(size))
    indices = np.linspace(0, size - 1, maximum, dtype=int).tolist()
    return sorted(set(indices + [0, size - 1]))


def save_geotiffs(ds: xr.Dataset, outdir: Path) -> list[Path]:
    try:
        import rasterio
        from affine import Affine
        from rasterio.control import GroundControlPoint
        from rasterio.crs import CRS
        from rasterio.transform import from_bounds
    except ImportError as exc:
        logging.warning("rasterio non installato: GeoTIFF non prodotti (%s)", exc)
        return []

    tif_dir = outdir / "geotiff"
    tif_dir.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []

    for timestamp in ds.time.values:
        da = ds["pyro_class"].sel(time=timestamp)
        arr, lat, lon, curvilinear = prepare_raster_array(da)

        profile: dict[str, object] = {
            "driver": "GTiff",
            "height": arr.shape[0],
            "width": arr.shape[1],
            "count": 1,
            "dtype": "uint8",
            "nodata": NODATA_CLASS,
            "compress": "DEFLATE",
            "predictor": 2,
        }
        gcps: list[GroundControlPoint] | None = None

        if not curvilinear:
            dx = abs(float(np.nanmedian(np.diff(lon))))
            dy = abs(float(np.nanmedian(np.diff(lat))))
            if not np.isfinite(dx) or not np.isfinite(dy) or dx <= 0 or dy <= 0:
                raise ValueError("Passo geografico non valido per la scrittura GeoTIFF.")
            west = float(np.nanmin(lon) - dx / 2.0)
            east = float(np.nanmax(lon) + dx / 2.0)
            south = float(np.nanmin(lat) - dy / 2.0)
            north = float(np.nanmax(lat) + dy / 2.0)
            profile.update(
                crs="EPSG:4326",
                transform=from_bounds(west, south, east, north, arr.shape[1], arr.shape[0]),
            )
        else:
            # Una griglia WRF in lat/lon è curvilinea e non possiede una singola
            # trasformazione affine esatta. Il GeoTIFF viene quindi georiferito
            # con una rete di Ground Control Points in EPSG:4326.
            profile.update(transform=Affine.identity())
            gcps = []
            for row in _sample_indices(arr.shape[0]):
                for col in _sample_indices(arr.shape[1]):
                    x = float(lon[row, col])
                    y = float(lat[row, col])
                    if np.isfinite(x) and np.isfinite(y):
                        gcps.append(
                            GroundControlPoint(
                                row=float(row) + 0.5,
                                col=float(col) + 0.5,
                                x=x,
                                y=y,
                            )
                        )
            if len(gcps) < 4:
                raise ValueError("Coordinate insufficienti per creare i GCP del GeoTIFF WRF.")

        tag = pd.Timestamp(timestamp).strftime("%Y%m%dT%H%MZ")
        path = tif_dir / f"pyro_class_{tag}.tif"
        with rasterio.open(path, "w", **profile) as dst:
            dst.write(arr, 1)
            if gcps is not None:
                dst.gcps = (gcps, CRS.from_epsg(4326))
            dst.update_tags(
                title="Atmospheric predisposition to pyroconvection",
                class_0=CLASS_NAMES[0],
                class_1=CLASS_NAMES[1],
                class_2=CLASS_NAMES[2],
                class_3=CLASS_NAMES[3],
                class_4=CLASS_NAMES[4],
                valid_time_utc=str(pd.Timestamp(timestamp)),
                georeferencing=(
                    "GCP EPSG:4326 for curvilinear WRF grid"
                    if curvilinear
                    else "Affine EPSG:4326"
                ),
            )
        outputs.append(path)
        logging.info("Salvato GeoTIFF: %s", path)
    return outputs


def matplotlib_objects():
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.patches import Patch

    cmap = ListedColormap([CLASS_COLORS[i] for i in range(5)])
    cmap.set_bad("#D9D9D9")
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], cmap.N)
    handles = [
        Patch(facecolor=CLASS_COLORS[i], edgecolor="0.35", label=CLASS_NAMES[i])
        for i in range(5)
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



def add_province_boundaries(
    ax,
    cartopy,
    *,
    enabled: bool = True,
    provinces_file: Path | None = None,
    linewidth: float = 0.55,
) -> None:
    """Disegna i confini provinciali ISTAT/Openpolis su assi Cartopy o semplici."""
    if not enabled or provinces_file is None:
        return
    try:
        geometries = load_boundary_geometries(provinces_file)
        if cartopy is not None:
            ccrs, _ = cartopy
            ax.add_geometries(
                geometries,
                crs=ccrs.PlateCarree(),
                facecolor="none",
                edgecolor="0.12",
                linewidth=linewidth,
                zorder=9,
            )
        else:
            draw_shapely_boundaries_plain(
                ax,
                geometries,
                linewidth=linewidth,
                color="0.12",
                zorder=9,
            )
    except Exception as exc:
        logging.warning("Impossibile disegnare i confini provinciali: %s", exc)


def plot_one_map(
    da: xr.DataArray,
    path: Path,
    use_cartopy: bool = True,
    draw_regional_boundaries: bool = False,
    regions_shapefile: Path | None = None,
    draw_province_boundaries: bool = True,
    provinces_file: Path | None = None,
    land_sea_mask: xr.DataArray | None = None,
    title_prefix: str = "Tipi di piroconvezione - WRF 1.5 km",
) -> None:
    plt, cmap, norm, handles = matplotlib_objects()
    cartopy = get_cartopy(use_cartopy)

    timestamp = pd.Timestamp(da.time.values)
    lon = da.longitude.values
    lat = da.latitude.values
    values = np.ma.masked_equal(da.values, NODATA_CLASS)
    if land_sea_mask is not None:
        sea = np.asarray(land_sea_mask.values, dtype=float) < DEFAULT_THRESHOLDS.min_land_fraction
        values = np.ma.masked_where(sea, values)

    if cartopy:
        ccrs, cfeature = cartopy
        fig = plt.figure(figsize=(9.2, 7.0))
        ax = plt.axes(projection=ccrs.PlateCarree())
        mesh = ax.pcolormesh(lon, lat, values, cmap=cmap, norm=norm, shading="nearest", transform=ccrs.PlateCarree())
        ax.coastlines(resolution="10m", linewidth=0.8, zorder=7)
        ax.add_feature(cfeature.BORDERS.with_scale("10m"), linewidth=0.45, zorder=7)
        add_regional_boundaries(
            ax,
            cartopy,
            enabled=draw_regional_boundaries,
            shapefile=regions_shapefile,
            linewidth=0.65,
        )
        add_province_boundaries(
            ax,
            cartopy,
            enabled=draw_province_boundaries,
            provinces_file=provinces_file,
            linewidth=0.65,
        )
        # OCEAN sopra il raster: garantisce una mascheratura visiva netta anche
        # lungo la costa, mentre le isole restano visibili come terra.
        ax.add_feature(cfeature.OCEAN.with_scale("10m"), facecolor="#B7DDE8", zorder=4)
        ax.add_feature(cfeature.LAKES.with_scale("10m"), facecolor="#B7DDE8", edgecolor="0.5", linewidth=0.3, zorder=5)
        ax.set_extent([float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max())])
        gl = ax.gridlines(draw_labels=True, linewidth=0.25, alpha=0.45)
        gl.top_labels = False
        gl.right_labels = False
    else:
        fig, ax = plt.subplots(figsize=(9.2, 7.0))
        mesh = ax.pcolormesh(lon, lat, values, cmap=cmap, norm=norm, shading="nearest")
        add_province_boundaries(
            ax,
            None,
            enabled=draw_province_boundaries,
            provinces_file=provinces_file,
            linewidth=0.65,
        )
        ax.set_xlabel("Longitudine")
        ax.set_ylabel("Latitudine")
        ax.grid(True, linewidth=0.25, alpha=0.45)

    ax.set_title(f"{title_prefix}\n{timestamp:%Y-%m-%d %H:%M UTC}", fontsize=13)
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 0.01))
    fig.text(
        0.5,
        0.055,
        "Predisposizione atmosferica sotto forcing di incendio intenso; WRF 1.5 km = previsione deterministica",
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
    draw_province_boundaries: bool = True,
    provinces_file: Path | None = None,
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
            draw_province_boundaries=draw_province_boundaries,
            provinces_file=provinces_file,
            land_sea_mask=ds["land_sea_mask"].sel(time=timestamp),
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
    draw_province_boundaries: bool = True,
    provinces_file: Path | None = None,
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
        values = np.ma.masked_equal(ds["pyro_class"].isel(time=idx).values, NODATA_CLASS)
        if "land_sea_mask" in ds:
            sea = np.asarray(ds["land_sea_mask"].isel(time=idx).values, dtype=float) < DEFAULT_THRESHOLDS.min_land_fraction
            values = np.ma.masked_where(sea, values)
        if cartopy:
            ccrs, cfeature = cartopy
            ax.pcolormesh(lon, lat, values, cmap=cmap, norm=norm, shading="nearest", transform=ccrs.PlateCarree())
            ax.coastlines(resolution="10m", linewidth=0.35, zorder=7)
            ax.add_feature(cfeature.BORDERS.with_scale("10m"), linewidth=0.2, zorder=7)
            add_regional_boundaries(
                ax,
                cartopy,
                enabled=draw_regional_boundaries,
                shapefile=regions_shapefile,
                linewidth=0.35,
            )
            add_province_boundaries(
                ax,
                cartopy,
                enabled=draw_province_boundaries,
                provinces_file=provinces_file,
                linewidth=0.38,
            )
            ax.add_feature(cfeature.OCEAN.with_scale("10m"), facecolor="#B7DDE8", zorder=4)
            ax.set_extent([float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max())])
        else:
            ax.pcolormesh(lon, lat, values, cmap=cmap, norm=norm, shading="nearest")
            add_province_boundaries(
                ax,
                None,
                enabled=draw_province_boundaries,
                provinces_file=provinces_file,
                linewidth=0.38,
            )
            ax.set_xticks([])
            ax.set_yticks([])

        ax.set_title(pd.Timestamp(timestamp).strftime("%Y-%m-%d %HZ"), fontsize=7.2, pad=2)

    fig.suptitle("Tipi di piroconvezione - WRF 1.5 km forecast", fontsize=13, y=0.985)
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 0.015), fontsize=8)
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
        description=("Elabora dati WRF 1.5 km ItaliaMeteo-CIMA-LaMMA con una ladder "
                     "adattiva alla risoluzione verticale reale e produce mappe di "
                     "predisposizione atmosferica alla piroconvezione."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--wrf-files", nargs="+", type=Path, help="File WRF NetCDF/GRIB contenenti surface e pressure levels")
    parser.add_argument("--download-wrf", action="store_true", help="Scarica automaticamente i dati open WRF da MeteoHub")
    parser.add_argument("--wrf-run", default="latest", help="Corsa da scaricare: latest oppure YYYYMMDDHH")
    parser.add_argument("--download-dir", type=Path, default=Path("dati_wrf"), help="Directory dei download WRF")
    parser.add_argument("--download-variables", nargs="+", default=list(DEFAULT_WRF_VARIABLES), help="Shortname MeteoHub da scaricare")
    parser.add_argument("--overwrite-downloads", action="store_true", help="Riscarica i file già presenti")
    parser.add_argument("--meteohub-wrf-base", default=METEOHUB_WRF_BASE, help=argparse.SUPPRESS)
    parser.add_argument("--surface-files", nargs="+", type=Path, help="File separati di superficie; in alternativa a --wrf-files")
    parser.add_argument("--pressure-files", nargs="+", type=Path, help="File separati ai livelli isobarici; in alternativa a --wrf-files")
    parser.add_argument(
        "--pressure-levels", nargs="+", default=[str(x) for x in DEFAULT_PRESSURE_LEVELS],
        help="Livelli isobarici da aprire/utilizzare, in hPa; accetta valori separati da spazio o virgola",
    )
    parser.add_argument(
        "--cfgrib-index-dir", type=Path, default=Path(".cfgrib_index"),
        help="Directory degli indici persistenti cfgrib; velocizza le esecuzioni successive",
    )
    parser.add_argument(
        "--bbox", nargs=4, type=float, metavar=("NORTH", "WEST", "SOUTH", "EAST"),
        help=(
            "Ritaglio geografico applicato a ogni file prima del merge; "
            "ad esempio Toscana: 44.6 9.7 42.2 12.6"
        ),
    )
    parser.add_argument(
        "--stride", type=int, default=1, metavar="N",
        help=(
            "Sottocampionamento spaziale applicato dopo la bbox e prima del merge: "
            "1=tutti i punti, 2=una cella ogni 2, 3=una ogni 3, ecc."
        ),
    )
    parser.add_argument(
        "--start",
        help=("Data valida iniziale YYYY-MM-DD. Deve essere usata insieme a --end; "
              "se entrambe sono omesse viene usato l’intero intervallo disponibile della corsa."),
    )
    parser.add_argument(
        "--end",
        help=("Data valida finale YYYY-MM-DD. Deve essere usata insieme a --start; "
              "se entrambe sono omesse viene usato l’intero intervallo disponibile della corsa."),
    )
    parser.add_argument("--hours", nargs="+", default=["00","03","06","09","12","15","18","21"], help="Ore valide UTC")
    parser.add_argument("--outdir", type=Path, default=Path("output_piroconvezione_wrf"))
    parser.add_argument("--keep-sea", action="store_true", help="Non maschera i punti marini")
    parser.add_argument(
        "--class-only", action="store_true",
        help=("Conserva classe, modalità adattiva, qualità, motivazione e land mask; "
              "non alloca le diagnostiche 3D estese"),
    )
    parser.add_argument("--no-individual-png", action="store_true")
    parser.add_argument("--no-panel-png", action="store_true")
    parser.add_argument("--panel-cols", type=int, default=8)
    parser.add_argument("--no-cartopy", action="store_true")
    parser.add_argument("--regional-boundaries", action="store_true", help="Disegna i confini amministrativi regionali nei PNG")
    parser.add_argument("--regions-shapefile", type=Path, help="Shapefile EPSG:4326 alternativo per i confini regionali")
    parser.add_argument(
        "--no-province-boundaries", action="store_true",
        help="Non disegna i confini provinciali nei PNG e nel pannello",
    )
    parser.add_argument(
        "--provinces-file", type=Path,
        help="GeoJSON/Shapefile EPSG:4326 dei confini provinciali; se omesso viene usata la copia ISTAT/Openpolis in cache",
    )
    parser.add_argument(
        "--boundaries-cache-dir", type=Path, default=DEFAULT_BOUNDARIES_CACHE_DIR,
        help="Directory cache per i confini amministrativi scaricati automaticamente",
    )
    parser.add_argument("--verbose", action="store_true")
    return parser


def validate_config(cfg: RunConfig) -> None:
    if (cfg.start is None) != (cfg.end is None):
        raise ValueError(
            "--start e --end devono essere specificati insieme oppure entrambi omessi"
        )
    if cfg.start is not None and cfg.end is not None:
        start = pd.Timestamp(cfg.start)
        end = pd.Timestamp(cfg.end)
        if end < start:
            raise ValueError("--end precede --start")
    if cfg.panel_cols < 1:
        raise ValueError("--panel-cols deve essere >= 1")
    if cfg.stride < 1:
        raise ValueError("--stride deve essere >= 1")
    if cfg.area.north <= cfg.area.south:
        raise ValueError("Nel --bbox NORTH deve essere maggiore di SOUTH")
    if cfg.area.east <= cfg.area.west:
        raise ValueError("Nel --bbox EAST deve essere maggiore di WEST")
    if cfg.regions_shapefile is not None and not cfg.regions_shapefile.exists():
        raise FileNotFoundError(f"Shapefile non trovato: {cfg.regions_shapefile}")
    if cfg.provinces_file is not None and not cfg.provinces_file.exists():
        raise FileNotFoundError(f"File dei confini provinciali non trovato: {cfg.provinces_file}")

def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    setup_logging(args.verbose)

    hours = normalize_hours(args.hours)
    pressure_levels = parse_pressure_levels(args.pressure_levels)
    area = Domain(*args.bbox) if args.bbox else Domain(90, -180, -90, 180)
    cfg = RunConfig(
        start=args.start, end=args.end, hours=hours, area=area,
        outdir=args.outdir, pressure_levels=pressure_levels, download=False,
        mask_sea=not args.keep_sea, make_individual_png=not args.no_individual_png,
        make_panel_png=not args.no_panel_png, panel_cols=args.panel_cols,
        use_cartopy=not args.no_cartopy,
        draw_regional_boundaries=(args.regional_boundaries or args.regions_shapefile is not None),
        regions_shapefile=args.regions_shapefile,
        draw_province_boundaries=not args.no_province_boundaries,
        provinces_file=args.provinces_file,
        boundaries_cache_dir=args.boundaries_cache_dir,
        stride=args.stride,
    )
    validate_config(cfg)
    cfg.outdir.mkdir(parents=True, exist_ok=True)

    # Il layer provinciale serve sia alla grafica sia come fallback robusto della
    # land-sea mask quando il prodotto WRF non contiene LSM/LANDMASK.
    province_boundaries_file = ensure_provinces_geojson(
        explicit_path=cfg.provinces_file,
        cache_dir=cfg.boundaries_cache_dir,
        required=False,
    ) if (cfg.mask_sea or cfg.draw_province_boundaries) else None

    downloaded_files: list[Path] = []
    if args.download_wrf:
        if args.wrf_files or args.surface_files or args.pressure_files:
            parser.error("--download-wrf non può essere combinato con --wrf-files/--surface-files/--pressure-files")
        downloaded_files, resolved_run = download_wrf_open_data(
            run=args.wrf_run, variables=args.download_variables,
            download_dir=args.download_dir, overwrite=args.overwrite_downloads,
            base_url=args.meteohub_wrf_base, pressure_levels=cfg.pressure_levels,
        )
        logging.info("Uso corsa WRF %s (%d file)", resolved_run, len(downloaded_files))

    if not downloaded_files and not args.wrf_files and not (args.surface_files and args.pressure_files):
        parser.error("Specificare --download-wrf, --wrf-files oppure la coppia --surface-files/--pressure-files")
    if downloaded_files:
        single_files, pressure_files = split_surface_pressure_files(downloaded_files)
    elif args.surface_files or args.pressure_files:
        if not args.surface_files or not args.pressure_files:
            parser.error("--surface-files e --pressure-files devono essere specificati insieme")
        single_files = list(args.surface_files)
        pressure_files = list(args.pressure_files)
    else:
        single_files, pressure_files = split_surface_pressure_files(list(args.wrf_files))

    pressure_files = filter_pressure_files(pressure_files, cfg.pressure_levels)

    ds_sl = open_many_netcdf(
        single_files, start=cfg.start, end=cfg.end, hours=cfg.hours,
        domain=cfg.area if args.bbox else None, stride=cfg.stride,
        cfgrib_index_dir=args.cfgrib_index_dir,
    )
    ds_pl = open_many_netcdf(
        pressure_files, start=cfg.start, end=cfg.end, hours=cfg.hours,
        domain=cfg.area if args.bbox else None, stride=cfg.stride,
        cfgrib_index_dir=args.cfgrib_index_dir,
    )
    ds_pl = select_pressure_levels(ds_pl, cfg.pressure_levels)

    # BBox e stride sono già stati applicati ai singoli file prima del merge.
    ds_sl, ds_pl = align_datasets(ds_sl, ds_pl)

    selected_times = pd.to_datetime(ds_sl["time"].values)
    if len(selected_times) == 0:
        raise ValueError("Nessuna scadenza WRF disponibile dopo il filtro temporale.")
    if cfg.start is None and cfg.end is None:
        logging.info(
            "Intervallo temporale automatico dalla corsa: %s -> %s UTC (%d scadenze)",
            selected_times.min().strftime("%Y-%m-%d %H:%M"),
            selected_times.max().strftime("%Y-%m-%d %H:%M"),
            len(selected_times),
        )
    else:
        logging.info(
            "Intervallo temporale richiesto: %s -> %s; scadenze disponibili selezionate: %s -> %s UTC (%d)",
            cfg.start, cfg.end,
            selected_times.min().strftime("%Y-%m-%d %H:%M"),
            selected_times.max().strftime("%Y-%m-%d %H:%M"),
            len(selected_times),
        )

    logging.info(
        "Livelli isobarici disponibili: %s hPa",
        ", ".join(f"{value:g}" for value in np.asarray(ds_pl["level"].values, dtype=float))
    )
    grid_y, grid_x, _, _, grid_curvilinear = horizontal_grid_info(ds_sl)
    logging.info(
        "Dominio: %d tempi, %d righe (%s), %d colonne (%s), %d livelli; griglia %s",
        ds_sl.sizes["time"],
        ds_sl.sizes[grid_y],
        grid_y,
        ds_sl.sizes[grid_x],
        grid_x,
        ds_pl.sizes["level"],
        "curvilinea" if grid_curvilinear else "rettilinea",
    )

    try:
        result = process_domain(
            ds_sl=ds_sl,
            ds_pl=ds_pl,
            thresholds=DEFAULT_THRESHOLDS,
            mask_sea=cfg.mask_sea,
            store_diagnostics=not args.class_only,
            fallback_land_boundaries=province_boundaries_file,
        )
    finally:
        # Chiude i file GRIB non appena l'output è stato materializzato.
        ds_sl.close()
        ds_pl.close()
        gc.collect()

    nc_path = cfg.outdir / "pyroconvection_wrf15.nc"
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
            draw_province_boundaries=cfg.draw_province_boundaries,
            provinces_file=province_boundaries_file,
        )
    if cfg.make_panel_png:
        save_panel_png(
            result,
            cfg.outdir / "pyroconvection_panel.png",
            cols=cfg.panel_cols,
            use_cartopy=cfg.use_cartopy,
            draw_regional_boundaries=cfg.draw_regional_boundaries,
            regions_shapefile=cfg.regions_shapefile,
            draw_province_boundaries=cfg.draw_province_boundaries,
            provinces_file=province_boundaries_file,
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
