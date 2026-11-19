#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pyroconvection_icon2i_meteohub.py
================================

Produce mappe previsionali della predisposizione atmosferica alla
piroconvezione usando ICON-2I distribuito da MeteoHub/ItaliaMeteo.

Classi:
    0 = nessuna classe / condizioni non sufficienti
    1 = Convective plume
    2 = Overshooting PyroCu potenziale
    3 = Resilient PyroCu potenziale
    4 = PyroCb / deep PyroCu potenziale

Il programma:
  1. individua un run ICON-2I nel repository open-data MeteoHub;
  2. scarica i GRIB delle variabili di superficie e ai livelli isobarici;
  3. seleziona area, scadenze e livelli di pressione richiesti;
  4. ricostruisce i profili verticali in quota AGL;
  5. stima l'ABL con il Bulk Richardson number;
  6. calcola LCL, stabilità ABL, gradiente di theta sopra ABL,
     quota del massimo wind shear e RH al top ABL;
  7. classifica ogni pixel;
  8. esporta NetCDF, GeoTIFF, PNG orari e tavola multi-pannello.
     I PNG possono includere i confini amministrativi regionali.

Fonte dati predefinita:
    https://meteohub.agenziaitaliameteo.it/nwp/
    ICON-2I_SURFACE_PRESSURE_LEVELS/

La modalità open-data non richiede credenziali, ma rende disponibili solo
i run presenti nel repository pubblico. I file sono memorizzati localmente
e riutilizzati alle esecuzioni successive.

IMPORTANTE
----------
Le classi 2-4 esprimono predisposizione atmosferica sotto il forcing di un
incendio sufficientemente intenso: non rappresentano una probabilità
calibrata di occorrenza. Le soglie restano "paper-informed" e devono essere
validate sui casi locali.

Installazione consigliata con conda:
    conda install -c conda-forge python=3.11 numpy pandas xarray netcdf4 \
        requests cfgrib eccodes rasterio matplotlib cartopy

Oppure con pip, purché ecCodes sia già disponibile nel sistema:
    pip install numpy pandas xarray netCDF4 requests cfgrib eccodes \
        rasterio matplotlib

Esempio, ultimo run disponibile e scadenze ogni 3 ore fino a +72 h:
    python pyroconvection_icon2i_meteohub.py \
        --run latest \
        --lead-hours 0 3 6 9 12 15 18 21 24 27 30 33 36 39 42 45 48 51 54 57 60 63 66 69 72 \
        --area 45.2 8.8 41.4 13.5 \
        --outdir output_piroconvezione_icon2i

Area: NORD OVEST SUD EST.
"""

from __future__ import annotations

import argparse
import logging
import math
import re
from dataclasses import dataclass
from html.parser import HTMLParser
from pathlib import Path
from typing import Iterable, Mapping, Sequence
from urllib.parse import urljoin

import numpy as np
import pandas as pd
import requests
import xarray as xr
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


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
SCRIPT_VERSION = "1.2.0"


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

@dataclass(frozen=True)
class RunConfig:
    run: str
    lead_hours: tuple[int, ...]
    area: Domain
    outdir: Path
    pressure_levels: tuple[int, ...]
    base_url: str
    download: bool = True
    overwrite_downloads: bool = False
    mask_sea: bool = True
    make_individual_png: bool = True
    make_panel_png: bool = True
    panel_cols: int = 8
    use_cartopy: bool = True
    draw_regional_boundaries: bool = False
    regions_shapefile: Path | None = None
    stride: int = 1
    timeout_s: int = 120


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
# MeteoHub open-data: catalogo run e download GRIB ICON-2I
# -----------------------------------------------------------------------------

DEFAULT_BASE_URL = (
    "https://meteohub.agenziaitaliameteo.it/"
    "nwp/ICON-2I_SURFACE_PRESSURE_LEVELS/"
)

# Cartella MeteoHub -> nome canonico interno.
SURFACE_VARIABLES: dict[str, str] = {
    "T_2M": "t2m",
    "TD_2M": "d2m",
    "U_10M": "u10",
    "V_10M": "v10",
    "PS": "sp",
    "HSURF": "hsurf",
    "FR_LAND": "lsm",
}

PRESSURE_VARIABLES: dict[str, str] = {
    "FI": "z",
    "T": "t",
    "QV": "q",
    "U": "u",
    "V": "v",
}

REQUIRED_VARIABLES: dict[str, str] = {
    **SURFACE_VARIABLES,
    **PRESSURE_VARIABLES,
}


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


def build_http_session() -> requests.Session:
    retry = Retry(
        total=5,
        connect=5,
        read=5,
        backoff_factor=1.0,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET", "HEAD"}),
        raise_on_status=False,
    )
    session = requests.Session()
    session.headers.update(
        {
            "User-Agent": (
                "pyroconvection-icon2i/1.0 "
                "(operational meteorological research; MeteoHub open data)"
            )
        }
    )
    session.mount("https://", HTTPAdapter(max_retries=retry))
    session.mount("http://", HTTPAdapter(max_retries=retry))
    return session


def list_directory_links(
    session: requests.Session,
    url: str,
    timeout_s: int = 120,
) -> list[str]:
    response = session.get(url, timeout=timeout_s)
    response.raise_for_status()
    parser = _HrefParser()
    parser.feed(response.text)
    return [urljoin(url, href) for href in parser.hrefs]


def list_available_runs(
    session: requests.Session,
    base_url: str,
    timeout_s: int = 120,
) -> list[str]:
    runs: list[str] = []
    for link in list_directory_links(session, base_url, timeout_s):
        match = re.search(r"/(\d{10})/?$", link)
        if match:
            runs.append(match.group(1))
    runs = sorted(set(runs))
    if not runs:
        raise RuntimeError(
            f"Nessun run ICON-2I trovato nell'indice MeteoHub: {base_url}"
        )
    return runs


def list_grib_files_for_variable(
    session: requests.Session,
    base_url: str,
    run: str,
    variable_folder: str,
    timeout_s: int = 120,
) -> list[str]:
    variable_url = urljoin(
        base_url.rstrip("/") + "/",
        f"{run}/{variable_folder}/",
    )
    links = list_directory_links(session, variable_url, timeout_s)
    gribs = [
        link
        for link in links
        if re.search(r"\.(grib2?|grb2?)$", link, flags=re.IGNORECASE)
    ]
    if not gribs:
        raise FileNotFoundError(
            f"Nessun file GRIB trovato per {variable_folder}: {variable_url}"
        )
    return sorted(set(gribs))


def download_file(
    session: requests.Session,
    url: str,
    target: Path,
    overwrite: bool = False,
    timeout_s: int = 120,
) -> Path:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and target.stat().st_size > 0 and not overwrite:
        logging.info("File già presente, salto download: %s", target)
        return target

    part = target.with_suffix(target.suffix + ".part")
    if part.exists():
        part.unlink()

    logging.info("Download %s", url)
    with session.get(url, stream=True, timeout=timeout_s) as response:
        response.raise_for_status()
        expected = int(response.headers.get("Content-Length", "0") or 0)
        downloaded = 0
        with part.open("wb") as stream:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if not chunk:
                    continue
                stream.write(chunk)
                downloaded += len(chunk)

    if expected and downloaded != expected:
        part.unlink(missing_ok=True)
        raise IOError(
            f"Download incompleto per {url}: {downloaded} byte su {expected}"
        )

    part.replace(target)
    if target.stat().st_size == 0:
        raise IOError(f"File scaricato ma vuoto: {target}")
    logging.info("Salvato %s (%.1f MiB)", target, target.stat().st_size / 2**20)
    return target


def download_icon2i_open_data(
    cfg: RunConfig,
) -> tuple[str, dict[str, list[Path]]]:
    """
    Scarica i file GRIB per il run richiesto dal repository open-data.

    Quando ``--run latest`` è usato, viene scelto il run più recente che
    espone tutte le cartelle/variabili necessarie. Questo evita di fermarsi
    su un run ancora in fase di pubblicazione.
    """
    session = build_http_session()
    available_runs = list_available_runs(
        session,
        cfg.base_url,
        cfg.timeout_s,
    )

    if cfg.run.strip().lower() == "latest":
        candidate_runs = list(reversed(available_runs))
    else:
        requested = cfg.run.strip()
        if not re.fullmatch(r"\d{10}", requested):
            raise ValueError("--run deve essere 'latest' oppure YYYYMMDDHH.")
        if requested not in available_runs:
            raise FileNotFoundError(
                f"Run {requested} non presente nel repository open-data. "
                f"Disponibili: {', '.join(available_runs)}"
            )
        candidate_runs = [requested]

    run: str | None = None
    urls_by_folder: dict[str, list[str]] | None = None
    last_error: Exception | None = None

    for candidate in candidate_runs:
        try:
            candidate_urls = {
                folder: list_grib_files_for_variable(
                    session=session,
                    base_url=cfg.base_url,
                    run=candidate,
                    variable_folder=folder,
                    timeout_s=cfg.timeout_s,
                )
                for folder in REQUIRED_VARIABLES
            }
        except Exception as exc:
            last_error = exc
            if cfg.run.strip().lower() == "latest":
                logging.warning(
                    "Run %s incompleto/non disponibile: %s. "
                    "Provo il precedente.",
                    candidate,
                    exc,
                )
                continue
            raise

        run = candidate
        urls_by_folder = candidate_urls
        break

    if run is None or urls_by_folder is None:
        raise RuntimeError(
            "Nessun run ICON-2I completo trovato nel repository MeteoHub."
        ) from last_error

    logging.info("Run ICON-2I selezionato: %s UTC", run)
    raw_dir = cfg.outdir / "raw" / run
    files: dict[str, list[Path]] = {}

    for folder, urls in urls_by_folder.items():
        local: list[Path] = []
        for url in urls:
            filename = url.rstrip("/").rsplit("/", 1)[-1]
            target = raw_dir / folder / filename
            local.append(
                download_file(
                    session,
                    url,
                    target,
                    overwrite=cfg.overwrite_downloads,
                    timeout_s=cfg.timeout_s,
                )
            )
        files[folder] = local

    return run, files


def discover_existing_icon_files(
    outdir: Path,
    requested_run: str,
) -> tuple[str, dict[str, list[Path]]]:
    raw_root = outdir / "raw"
    if requested_run.lower() == "latest":
        runs = sorted(
            p.name
            for p in raw_root.glob("[0-9]" * 10)
            if p.is_dir() and re.fullmatch(r"\d{10}", p.name)
        )
        if not runs:
            raise FileNotFoundError(
                f"Nessun run locale trovato in {raw_root}. "
                "Rimuovere --no-download oppure specificare --run."
            )
        run = runs[-1]
    else:
        run = requested_run

    run_dir = raw_root / run
    files: dict[str, list[Path]] = {}
    for folder in REQUIRED_VARIABLES:
        folder_dir = run_dir / folder
        candidates = sorted(
            [
                *folder_dir.glob("*.grib"),
                *folder_dir.glob("*.grib2"),
                *folder_dir.glob("*.grb"),
                *folder_dir.glob("*.grb2"),
            ]
        )
        if not candidates:
            raise FileNotFoundError(
                f"File GRIB mancanti per {folder} in {folder_dir}"
            )
        files[folder] = candidates
    return run, files


def parse_lead_hours(values: Iterable[str | int]) -> tuple[int, ...]:
    leads: list[int] = []
    for item in values:
        value = int(item)
        if value < 0:
            raise ValueError("Le scadenze --lead-hours devono essere >= 0.")
        if value not in leads:
            leads.append(value)
    if not leads:
        raise ValueError("Lista --lead-hours vuota.")
    return tuple(sorted(leads))


# -----------------------------------------------------------------------------
# Apertura e normalizzazione dataset ICON-2I GRIB
# -----------------------------------------------------------------------------

def _rename_common(ds: xr.Dataset) -> xr.Dataset:
    rename: dict[str, str] = {}
    candidates = {
        "isobaricInhPa": "level",
        "pressure_level": "level",
        "lat": "latitude",
        "lon": "longitude",
    }
    for old, new in candidates.items():
        if (old in ds.dims or old in ds.coords) and not (
            new in ds.dims or new in ds.coords
        ):
            rename[old] = new
    if rename:
        ds = ds.rename(rename)
    return ds


def _normalize_valid_time(ds: xr.Dataset) -> xr.Dataset:
    """
    Converte la rappresentazione cfgrib (reference time + step + valid_time)
    in una dimensione temporale unica chiamata ``time``.
    """
    if "valid_time" not in ds.coords:
        if "time" in ds.coords:
            if ds["time"].ndim == 0:
                ds = ds.expand_dims(time=[ds["time"].values])
            return ds.sortby("time")
        raise ValueError("Il GRIB non contiene né valid_time né time.")

    valid = ds["valid_time"]

    if valid.ndim == 0:
        value = valid.values
        if "time" in ds.coords:
            ds = ds.drop_vars("time")
        ds = ds.drop_vars("valid_time")
        return ds.expand_dims(time=[value])

    if valid.ndim != 1:
        raise ValueError(
            "valid_time multidimensionale non supportato. "
            "Usare un singolo run ICON-2I per file."
        )

    source_dim = valid.dims[0]
    values = np.asarray(valid.values)

    if source_dim == "time":
        ds = ds.drop_vars("valid_time")
        ds = ds.assign_coords(time=("time", values))
        ds = ds.sortby("time")
        _, unique_idx = np.unique(ds["time"].values, return_index=True)
        return ds.isel(time=np.sort(unique_idx))

    if "time" in ds.coords:
        ds = ds.drop_vars("time")
    if "valid_time" in ds.coords:
        ds = ds.drop_vars("valid_time")

    ds = ds.assign_coords(_forecast_time=(source_dim, values))
    ds = ds.swap_dims({source_dim: "_forecast_time"})
    ds = ds.rename({"_forecast_time": "time"})
    if source_dim in ds.coords:
        ds = ds.drop_vars(source_dim)

    ds = ds.sortby("time")
    _, unique_idx = np.unique(ds["time"].values, return_index=True)
    return ds.isel(time=np.sort(unique_idx))


def _choose_data_variable(ds: xr.Dataset, source_folder: str) -> str:
    variables = list(ds.data_vars)
    if len(variables) == 1:
        return variables[0]

    # Alcuni GRIB possono esporre coordinate ausiliarie come data variables.
    candidates = [
        name
        for name in variables
        if ds[name].ndim >= 2
        and name.lower() not in {"valid_time", "step", "number"}
    ]
    if len(candidates) == 1:
        return candidates[0]
    raise ValueError(
        f"Impossibile identificare univocamente la variabile in {source_folder}. "
        f"Data variables trovate: {variables}"
    )


def normalize_icon_units(
    da: xr.DataArray,
    canonical_name: str,
) -> xr.DataArray:
    """Normalizza le poche unità critiche usate dalla diagnostica."""
    attrs = dict(da.attrs)
    units = str(attrs.get("units", "")).strip()
    units_l = units.lower().replace(" ", "")
    standard_name = str(attrs.get("standard_name", "")).lower()

    if canonical_name == "z":
        # process_profile richiede geopotenziale [m2 s-2].
        if (
            "geopotential_height" in standard_name
            or units_l in {"m", "gpm", "geopotentialmetre", "geopotentialmeter"}
        ):
            da = da * G
            attrs["units"] = "m2 s-2"
            attrs["conversion_note"] = (
                "Converted from geopotential height to geopotential using g"
            )

    elif canonical_name == "hsurf":
        # process_profile richiede altezza geometrica [m].
        if (
            "geopotential" in standard_name
            and "height" not in standard_name
        ) or any(token in units_l for token in ("m2s-2", "m**2s**-2")):
            da = da / G
            attrs["units"] = "m"
            attrs["conversion_note"] = (
                "Converted from geopotential to geometric height using g"
            )

    elif canonical_name == "sp":
        # process_profile richiede Pa.
        if units_l in {"hpa", "mbar", "mb", "millibar"}:
            da = da * 100.0
            attrs["units"] = "Pa"
            attrs["conversion_note"] = "Converted from hPa/mbar to Pa"

    elif canonical_name == "lsm":
        # La soglia min_land_fraction è espressa fra 0 e 1.
        if "%" in units or units_l in {"percent", "percentage"}:
            da = da / 100.0
            attrs["units"] = "1"
            attrs["conversion_note"] = "Converted from percent to fraction"

    da.attrs = attrs
    return da


def _infer_pressure_level_hpa(
    ds: xr.Dataset,
    da: xr.DataArray,
    path: Path,
) -> float:
    """
    Ricava il livello isobarico del messaggio GRIB.

    Nei file MeteoHub ogni GRIB di pressione contiene normalmente un solo
    livello. cfgrib lo espone quindi come coordinata *scalare*
    ``isobaricInhPa`` (rinominata ``level``), non come dimensione. Se tale
    coordinata viene eliminata prima del concatenamento, il dataset finale
    perde completamente l'asse verticale.

    Il nome del file viene usato esclusivamente come fallback robusto.
    """
    for container in (da, ds):
        if "level" in container.coords:
            values = np.asarray(container.coords["level"].values).reshape(-1)
            finite = values[np.isfinite(values.astype(float))]
            if finite.size == 1:
                return float(finite[0])

    # Attributi GRIB talvolta esposti da cfgrib/ecCodes.
    for attrs in (da.attrs, ds.attrs):
        for key in ("GRIB_level", "level", "pressureLevel"):
            if key in attrs:
                try:
                    return float(attrs[key])
                except (TypeError, ValueError):
                    pass

    match = re.search(
        r"(?:isobaricInhPa|pressure[_-]?level)[-_](\d+(?:\.\d+)?)",
        path.name,
        flags=re.IGNORECASE,
    )
    if match:
        return float(match.group(1))

    raise ValueError(
        f"Impossibile ricavare il livello isobarico dal GRIB {path}. "
        f"Coordinate disponibili: {list(ds.coords)}"
    )


def _promote_scalar_pressure_level(
    ds: xr.Dataset,
    da: xr.DataArray,
    path: Path,
    canonical_name: str,
) -> xr.DataArray:
    """Trasforma il livello isobarico scalare in una vera dimensione."""
    if canonical_name not in set(PRESSURE_VARIABLES.values()):
        return da

    if "level" in da.dims:
        return da

    level_hpa = _infer_pressure_level_hpa(ds, da, path)

    # Rimuove l'eventuale coordinata scalare prima di creare la dimensione.
    if "level" in da.coords:
        da = da.reset_coords("level", drop=True)

    da = da.expand_dims(level=[level_hpa])
    da["level"].attrs.update(
        {
            "long_name": "pressure level",
            "standard_name": "air_pressure",
            "units": "hPa",
            "positive": "down",
            "axis": "Z",
        }
    )
    logging.debug(
        "Livello isobarico preservato per %s: %.0f hPa",
        path.name,
        level_hpa,
    )
    return da


def open_one_icon_grib(
    path: Path,
    source_folder: str,
    canonical_name: str,
) -> xr.Dataset:
    logging.info("Apro GRIB %s", path)
    try:
        ds = xr.open_dataset(
            path,
            engine="cfgrib",
            backend_kwargs={
                "indexpath": "",
                "errors": "raise",
            },
        )
    except Exception as exc:
        raise RuntimeError(
            f"Errore aprendo {path} con cfgrib. "
            "Verificare installazione di cfgrib/ecCodes."
        ) from exc

    ds = _rename_common(ds)
    ds = _normalize_valid_time(ds)
    source_name = _choose_data_variable(ds, source_folder)
    da = ds[source_name].rename(canonical_name)
    da = normalize_icon_units(da, canonical_name)

    # IMPORTANTE: nei file MeteoHub il livello isobarico è una coordinata
    # scalare. Va promosso a dimensione PRIMA della pulizia delle coordinate.
    da = _promote_scalar_pressure_level(
        ds=ds,
        da=da,
        path=path,
        canonical_name=canonical_name,
    )

    # Elimina coordinate scalari specifiche del GRIB (ad esempio
    # heightAboveGround=2/10 m) che possono entrare in conflitto nel merge.
    # Le coordinate dimensionali time/level/latitude/longitude sono preservate.
    auxiliary_coords = [
        name
        for name in da.coords
        if name not in da.dims
    ]
    if auxiliary_coords:
        da = da.reset_coords(auxiliary_coords, drop=True)

    # Topografia e land fraction sono statiche: evita che un solo timestamp
    # limiti il merge con le variabili previsionali.
    if canonical_name in {"hsurf", "lsm"} and "time" in da.dims:
        da = da.isel(time=0, drop=True)

    return da.to_dataset()


def open_variable_gribs(
    paths: Sequence[Path],
    source_folder: str,
    canonical_name: str,
) -> xr.Dataset:
    if not paths:
        raise FileNotFoundError(f"Nessun GRIB specificato per {source_folder}.")

    datasets = [
        open_one_icon_grib(path, source_folder, canonical_name)
        for path in paths
    ]

    if len(datasets) == 1:
        return datasets[0]

    if canonical_name in {"hsurf", "lsm"}:
        return datasets[0]

    try:
        out = xr.combine_by_coords(
            datasets,
            combine_attrs="drop_conflicts",
        )
    except Exception as exc:
        if canonical_name in set(PRESSURE_VARIABLES.values()):
            logging.warning(
                "combine_by_coords fallito per %s (%s): "
                "concateno esplicitamente lungo level.",
                source_folder,
                exc,
            )
            out = xr.concat(
                datasets,
                dim="level",
                data_vars="minimal",
                coords="minimal",
                compat="override",
                join="outer",
            )
        else:
            # Fallback per file suddivisi esclusivamente lungo il tempo.
            out = xr.concat(
                datasets,
                dim="time",
                data_vars="minimal",
                coords="minimal",
                compat="override",
            )

    if "level" in out.coords:
        out = out.sortby("level", ascending=False)
        _, unique_level_idx = np.unique(out["level"].values, return_index=True)
        out = out.isel(level=np.sort(unique_level_idx))

    if "time" in out.coords:
        out = out.sortby("time")
        _, unique_idx = np.unique(out["time"].values, return_index=True)
        out = out.isel(time=np.sort(unique_idx))
    return out


def load_icon2i_datasets(
    files_by_folder: Mapping[str, Sequence[Path]],
) -> tuple[xr.Dataset, xr.Dataset]:
    surface_parts: list[xr.Dataset] = []
    pressure_parts: list[xr.Dataset] = []

    for folder, canonical in SURFACE_VARIABLES.items():
        surface_parts.append(
            open_variable_gribs(files_by_folder[folder], folder, canonical)
        )

    for folder, canonical in PRESSURE_VARIABLES.items():
        pressure_parts.append(
            open_variable_gribs(files_by_folder[folder], folder, canonical)
        )

    ds_sl = xr.merge(
        surface_parts,
        compat="override",
        join="inner",
        combine_attrs="drop_conflicts",
    )
    ds_pl = xr.merge(
        pressure_parts,
        compat="override",
        join="inner",
        combine_attrs="drop_conflicts",
    )
    return ds_sl, ds_pl


def ensure_lat_lon_order(ds: xr.Dataset) -> xr.Dataset:
    if "latitude" not in ds.coords or "longitude" not in ds.coords:
        raise ValueError(
            "Coordinate latitude/longitude non trovate. "
            f"Coordinate disponibili: {list(ds.coords)}"
        )
    if ds["latitude"].ndim != 1 or ds["longitude"].ndim != 1:
        raise ValueError(
            "Questo script richiede ICON-2I su griglia regolare lat/lon 1D. "
            "È stata trovata una griglia curvilinea/non strutturata."
        )

    if ds.latitude.values[0] < ds.latitude.values[-1]:
        ds = ds.sortby("latitude", ascending=False)
    if ds.longitude.values[0] > ds.longitude.values[-1]:
        ds = ds.sortby("longitude", ascending=True)
    return ds


def subset_domain(ds: xr.Dataset, area: Domain, stride: int = 1) -> xr.Dataset:
    ds = ensure_lat_lon_order(ds)
    lon = ds.longitude

    # Converte la bbox in 0..360 solo se la griglia usa tale convenzione.
    west, east = area.west, area.east
    if float(lon.min()) >= 0.0 and west < 0.0:
        west %= 360.0
        east %= 360.0

    lat_mask = (ds.latitude <= area.north) & (ds.latitude >= area.south)
    lon_mask = (ds.longitude >= west) & (ds.longitude <= east)
    ds = ds.where(lat_mask, drop=True).where(lon_mask, drop=True)

    if ds.sizes.get("latitude", 0) == 0 or ds.sizes.get("longitude", 0) == 0:
        raise ValueError(
            "La bbox non interseca la griglia ICON-2I: "
            f"N={area.north}, W={area.west}, S={area.south}, E={area.east}"
        )

    if stride > 1:
        ds = ds.isel(
            latitude=slice(None, None, stride),
            longitude=slice(None, None, stride),
        )
    return ds


def select_pressure_levels(
    ds: xr.Dataset,
    requested_levels: Sequence[int],
) -> xr.Dataset:
    if "level" not in ds.coords:
        raise ValueError(
            "Coordinata isobarica 'level' non trovata dopo il caricamento. "
            "Questo indica che i livelli scalari dei GRIB non sono stati "
            "promossi a dimensione. "
            f"Coordinate disponibili: {list(ds.coords)}"
        )

    available = np.asarray(ds.level.values, dtype=float)
    selected = [
        int(level)
        for level in requested_levels
        if np.any(np.isclose(available, float(level)))
    ]
    missing = [
        int(level)
        for level in requested_levels
        if not np.any(np.isclose(available, float(level)))
    ]
    if missing:
        logging.warning(
            "Livelli richiesti non presenti in ICON-2I e ignorati: %s",
            ", ".join(map(str, missing)),
        )
    if len(selected) < 4:
        raise ValueError(
            f"Troppi pochi livelli disponibili per i profili: {selected}"
        )

    ds = ds.sel(level=selected)
    if ds.level.values[0] < ds.level.values[-1]:
        ds = ds.sortby("level", ascending=False)
    return ds


def select_lead_times(
    ds: xr.Dataset,
    run: str,
    lead_hours: Sequence[int],
) -> xr.Dataset:
    if "time" not in ds.coords:
        # Dataset statico.
        return ds

    run_time = pd.to_datetime(run, format="%Y%m%d%H")
    wanted = np.array(
        [
            np.datetime64(run_time + pd.Timedelta(hours=int(lead)))
            for lead in lead_hours
        ]
    )
    available = np.asarray(ds.time.values).astype("datetime64[ns]")
    wanted = wanted.astype("datetime64[ns]")
    present = np.isin(available, wanted)
    selected = ds.isel(time=np.where(present)[0])

    actual = {
        int(
            round(
                (pd.Timestamp(value) - run_time)
                / pd.Timedelta(hours=1)
            )
        )
        for value in selected.time.values
    }
    missing = sorted(set(map(int, lead_hours)) - actual)
    if missing:
        logging.warning(
            "Scadenze non disponibili nel run %s e ignorate: %s",
            run,
            ", ".join(f"+{value}h" for value in missing),
        )
    if selected.sizes.get("time", 0) == 0:
        raise ValueError(
            f"Nessuna scadenza richiesta presente nel run ICON-2I {run}."
        )
    return selected


def align_datasets(
    ds_sl: xr.Dataset,
    ds_pl: xr.Dataset,
) -> tuple[xr.Dataset, xr.Dataset]:
    ds_sl = ensure_lat_lon_order(ds_sl)
    ds_pl = ensure_lat_lon_order(ds_pl)

    common_time = np.intersect1d(ds_sl.time.values, ds_pl.time.values)
    if len(common_time) == 0:
        raise ValueError(
            "Nessun tempo di validità comune fra campi di superficie "
            "e campi ai livelli isobarici."
        )

    ds_sl = ds_sl.sel(time=common_time)
    ds_pl = ds_pl.sel(time=common_time)

    same_lat = np.array_equal(ds_sl.latitude.values, ds_pl.latitude.values)
    same_lon = np.array_equal(ds_sl.longitude.values, ds_pl.longitude.values)
    if not (same_lat and same_lon):
        logging.warning(
            "Griglie ICON non coincidenti: interpolo i livelli isobarici "
            "sulla griglia dei campi di superficie."
        )
        ds_pl = ds_pl.interp(
            latitude=ds_sl.latitude,
            longitude=ds_sl.longitude,
        )

    return ds_sl, ds_pl


def find_var(
    ds: xr.Dataset,
    candidates: Sequence[str],
) -> xr.DataArray:
    for name in candidates:
        if name in ds.data_vars:
            return ds[name]
    raise KeyError(
        f"Nessuna variabile tra {list(candidates)}. "
        f"Presenti: {list(ds.data_vars)}"
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
    Per il profilo ambientale ICON-2I si usa il riferimento di superficie; la
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
    # non ricavabile dal solo ICON-2I. La classe indica quindi predisposizione.
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
    surface_height_m: float,
    surface_pressure_pa: float,
    t2m_k: float,
    td2m_k: float,
    u10: float,
    v10: float,
    thresholds: Thresholds,
) -> ProfileMetrics:
    metrics = ProfileMetrics()

    scalar_values = [
        surface_height_m,
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
    orog_m = surface_height_m
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
        fallback_blh=np.nan,
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
    run: str,
    mask_sea: bool = True,
) -> xr.Dataset:
    # Nomi canonici assegnati in fase di lettura dei GRIB MeteoHub.
    t2m = find_var(ds_sl, ["t2m"])
    d2m = find_var(ds_sl, ["d2m"])
    u10_da = find_var(ds_sl, ["u10"])
    v10_da = find_var(ds_sl, ["v10"])
    sp = find_var(ds_sl, ["sp"])
    hsurf = find_var(ds_sl, ["hsurf"])
    lsm = find_var(ds_sl, ["lsm"])

    z_pl = find_var(ds_pl, ["z"])
    t_pl = find_var(ds_pl, ["t"])
    q_pl = find_var(ds_pl, ["q"])
    u_pl = find_var(ds_pl, ["u"])
    v_pl = find_var(ds_pl, ["v"])

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
    hsurf = broadcast_sl(hsurf)
    lsm = broadcast_sl(lsm)

    times = ds_sl.time.values
    lats = ds_sl.latitude.values
    lons = ds_sl.longitude.values
    levels = ds_pl.level.values.astype(float)

    nt, ny, nx = len(times), len(lats), len(lons)
    shape = (nt, ny, nx)

    pyro_class = np.full(shape, NODATA_CLASS, dtype=np.uint8)
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
    hsurf_np = np.asarray(hsurf.values)
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

                # 0 = punto terrestre valido ma nessuna classe identificata.
                pyro_class[it, iy, ix] = np.uint8(0)

                m = process_profile(
                    pressure_hpa=levels,
                    geopotential=z_pl_np[it, :, iy, ix],
                    temperature_k=t_pl_np[it, :, iy, ix],
                    q=q_pl_np[it, :, iy, ix],
                    u=u_pl_np[it, :, iy, ix],
                    v=v_pl_np[it, :, iy, ix],
                    surface_height_m=float(hsurf_np[it, iy, ix]),
                    surface_pressure_pa=float(sp_np[it, iy, ix]),
                    t2m_k=float(t2m_np[it, iy, ix]),
                    td2m_k=float(d2m_np[it, iy, ix]),
                    u10=float(u10_np[it, iy, ix]),
                    v10=float(v10_np[it, iy, ix]),
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

    run_time = np.datetime64(
        pd.to_datetime(run, format="%Y%m%d%H"),
        "ns",
    )
    lead_time_hours = (
        (times.astype("datetime64[ns]") - run_time)
        / np.timedelta64(1, "h")
    ).astype(np.int16)

    coords = {
        "time": times,
        "latitude": lats,
        "longitude": lons,
        "forecast_reference_time": run_time,
    }
    dims = ("time", "latitude", "longitude")

    out = xr.Dataset(
        data_vars={
            "lead_time_hours": (("time",), lead_time_hours),
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
            "title": "ICON-2I forecast atmospheric predisposition to pyroconvection",
            "source": "MeteoHub / ItaliaMeteo - ICON-2I surface and pressure levels",
            "model_run_utc": run,
            "method": "Paper-informed diagnostic based on ABL stability, LCL/ABL, shear proximity, free-tropospheric theta gradient and RH at ABL top",
            "reference": "Castellnou et al. (2022), JGR Atmospheres, doi:10.1029/2022JD036920",
            "warning": "Classes indicate atmospheric predisposition under an intense-fire forcing, not calibrated occurrence probability.",
            "class_0": CLASS_NAMES[0],
            "class_1": CLASS_NAMES[1],
            "class_2": CLASS_NAMES[2],
            "class_3": CLASS_NAMES[3],
            "class_4": CLASS_NAMES[4],
        },
    )

    out["lead_time_hours"].attrs.update(
        long_name="forecast lead time",
        units="h",
    )
    out["forecast_reference_time"].attrs.update(
        long_name="ICON-2I forecast reference time",
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
    if path.exists():
        path.unlink()

    engines = xr.backends.list_engines()
    if "netcdf4" in engines:
        engine = "netcdf4"
    elif "h5netcdf" in engines:
        engine = "h5netcdf"
    else:
        engine = "scipy"

    ds_write = ds.copy(deep=False)
    for name in ds_write.variables:
        ds_write[name].encoding = {}
        ds_write[name].attrs.pop("_FillValue", None)

    encoding: dict[str, dict] = {}
    if engine in {"netcdf4", "h5netcdf"}:
        for name in ds_write.data_vars:
            if name == "pyro_class":
                encoding[name] = {
                    "zlib": True,
                    "complevel": 4,
                    "dtype": "u1",
                    "_FillValue": NODATA_CLASS,
                }
            elif name == "lead_time_hours":
                encoding[name] = {
                    "zlib": True,
                    "complevel": 4,
                    "dtype": "i2",
                }
            else:
                encoding[name] = {
                    "zlib": True,
                    "complevel": 4,
                    "dtype": "f4",
                }
    else:
        # NetCDF3/Scipy non supporta uint8: usa int16 e una codifica
        # temporale esplicita, mantenendo comunque il valore nodata 255.
        ds_write["pyro_class"] = ds_write["pyro_class"].astype(np.int16)
        ds_write["pyro_class"].attrs.update(ds["pyro_class"].attrs)
        ds_write["pyro_class"].attrs["flag_values"] = np.array(
            [0, 1, 2, 3, 4],
            dtype=np.int16,
        )
        for name in ds_write.data_vars:
            if name == "pyro_class":
                encoding[name] = {
                    "dtype": "i2",
                    "_FillValue": np.int16(NODATA_CLASS),
                }
            elif name == "lead_time_hours":
                encoding[name] = {"dtype": "i2"}
            else:
                encoding[name] = {"dtype": "f4"}
        time_encoding = {
            "dtype": "i4",
            "units": "hours since 1970-01-01 00:00:00",
            "calendar": "proleptic_gregorian",
        }
        encoding["time"] = dict(time_encoding)
        encoding["forecast_reference_time"] = dict(time_encoding)

    ds_write.to_netcdf(path, engine=engine, encoding=encoding)
    logging.info("Salvato NetCDF (%s): %s", engine, path)


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
    cmap.set_bad(alpha=0.0)
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
    Shapefile, geometrie in coordinate geografiche). In caso contrario usa il
    layer Natural Earth di primo livello amministrativo alla scala 1:10m.

    Natural Earth viene scaricato e memorizzato nella cache Cartopy al primo uso.
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
    title_prefix: str = "Tipi di piroconvezione - ICON-2I Forecast",
) -> None:
    plt, cmap, norm, handles = matplotlib_objects()
    cartopy = get_cartopy(use_cartopy)

    timestamp = pd.Timestamp(da.time.values)
    lon = da.longitude.values
    lat = da.latitude.values
    values = np.ma.masked_equal(da.values, NODATA_CLASS)

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
        "Predisposizione atmosferica sotto forcing di incendio intenso; ICON-2I = previsione",
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
        values = np.ma.masked_equal(
            ds["pyro_class"].isel(time=idx).values,
            NODATA_CLASS,
        )
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

    fig.suptitle("Tipi di piroconvezione - ICON-2I Forecast", fontsize=13, y=0.985)
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
        description=(
            "Scarica ICON-2I da MeteoHub e produce mappe previsionali "
            "di predisposizione atmosferica alla piroconvezione."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--run",
        default="latest",
        help="Run ICON-2I UTC: latest oppure YYYYMMDDHH",
    )
    parser.add_argument(
        "--lead-hours",
        nargs="+",
        default=[str(value) for value in range(0, 73, 3)],
        help="Scadenze previsionali in ore dal run",
    )
    parser.add_argument(
        "--area",
        nargs=4,
        type=float,
        metavar=("NORTH", "WEST", "SOUTH", "EAST"),
        default=[45.2, 8.8, 41.4, 13.5],
        help="Bounding box: nord ovest sud est",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("output_piroconvezione_icon2i"),
    )
    parser.add_argument(
        "--base-url",
        default=DEFAULT_BASE_URL,
        help="Indice open-data MeteoHub del dataset ICON-2I",
    )
    parser.add_argument(
        "--pressure-levels",
        nargs="+",
        default=[
            "1000", "975", "950", "925", "900", "850", "800", "750",
            "700", "650", "600", "550", "500", "450", "400", "350",
            "300", "250", "200",
        ],
        help="Livelli isobarici hPa, separati da spazi o virgole",
    )
    parser.add_argument(
        "--no-download",
        action="store_true",
        help="Usa i GRIB già presenti in OUTDIR/raw/RUN",
    )
    parser.add_argument(
        "--overwrite-downloads",
        action="store_true",
        help="Riscarica anche i GRIB già presenti",
    )
    parser.add_argument(
        "--list-runs",
        action="store_true",
        help="Elenca i run disponibili nel repository MeteoHub e termina",
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help=(
            "Sottocampionamento spaziale: 1 = risoluzione nativa, "
            "2 = un punto ogni due, ecc."
        ),
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=120,
        help="Timeout HTTP in secondi",
    )
    parser.add_argument(
        "--keep-sea",
        action="store_true",
        help="Non maschera i punti marini",
    )
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


def validate_config(cfg: RunConfig) -> None:
    if not (-90 <= cfg.area.south < cfg.area.north <= 90):
        raise ValueError("Latitudini area non valide.")
    if not (-180 <= cfg.area.west < cfg.area.east <= 360):
        raise ValueError("Longitudini area non valide.")
    if cfg.panel_cols < 1:
        raise ValueError("--panel-cols deve essere >= 1")
    if cfg.stride < 1:
        raise ValueError("--stride deve essere >= 1")
    if cfg.timeout_s < 1:
        raise ValueError("--timeout deve essere >= 1")
    if cfg.regions_shapefile is not None and not cfg.regions_shapefile.exists():
        raise FileNotFoundError(
            f"Shapefile dei confini regionali non trovato: {cfg.regions_shapefile}"
        )


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    setup_logging(args.verbose)

    if args.list_runs:
        session = build_http_session()
        for run in list_available_runs(
            session,
            args.base_url,
            args.timeout,
        ):
            print(run)
        return 0

    lead_hours = parse_lead_hours(args.lead_hours)
    pressure_levels = parse_pressure_levels(args.pressure_levels)
    area = Domain(*args.area)
    cfg = RunConfig(
        run=args.run,
        lead_hours=lead_hours,
        area=area,
        outdir=args.outdir,
        pressure_levels=pressure_levels,
        base_url=args.base_url.rstrip("/") + "/",
        download=not args.no_download,
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
        stride=args.stride,
        timeout_s=args.timeout,
    )
    validate_config(cfg)
    cfg.outdir.mkdir(parents=True, exist_ok=True)

    if cfg.download:
        run, files_by_folder = download_icon2i_open_data(cfg)
    else:
        run, files_by_folder = discover_existing_icon_files(
            cfg.outdir,
            cfg.run,
        )

    logging.info("Caricamento run ICON-2I %s UTC", run)
    ds_sl, ds_pl = load_icon2i_datasets(files_by_folder)

    ds_sl = select_lead_times(ds_sl, run, cfg.lead_hours)
    ds_pl = select_lead_times(ds_pl, run, cfg.lead_hours)
    ds_pl = select_pressure_levels(ds_pl, cfg.pressure_levels)

    ds_sl = subset_domain(ds_sl, cfg.area, cfg.stride)
    ds_pl = subset_domain(ds_pl, cfg.area, cfg.stride)
    ds_sl, ds_pl = align_datasets(ds_sl, ds_pl)

    logging.info(
        "Dominio: %d tempi, %d latitudini, %d longitudini, %d livelli",
        ds_sl.sizes["time"],
        ds_sl.sizes["latitude"],
        ds_sl.sizes["longitude"],
        ds_pl.sizes["level"],
    )
    profile_count = (
        ds_sl.sizes["time"]
        * ds_sl.sizes["latitude"]
        * ds_sl.sizes["longitude"]
    )
    if profile_count > 300_000:
        logging.warning(
            "Sono richiesti %,d profili atmosferici. Alla risoluzione "
            "ICON-2I nativa l'elaborazione può essere lunga; per un primo "
            "test usare --stride 2 o --stride 3.",
            profile_count,
        )

    result = process_domain(
        ds_sl=ds_sl,
        ds_pl=ds_pl,
        thresholds=DEFAULT_THRESHOLDS,
        run=run,
        mask_sea=cfg.mask_sea,
    )

    nc_path = cfg.outdir / f"pyroconvection_icon2i_{run}.nc"
    save_netcdf(result, nc_path)
    save_geotiffs(result, cfg.outdir)
    save_class_csv(
        result,
        cfg.outdir / "classi_piroconvezione.csv",
    )

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
            cfg.outdir / f"pyroconvection_icon2i_panel_{run}.png",
            cols=cfg.panel_cols,
            use_cartopy=cfg.use_cartopy,
            draw_regional_boundaries=cfg.draw_regional_boundaries,
            regions_shapefile=cfg.regions_shapefile,
        )

    logging.info(
        "Elaborazione ICON-2I completata in %s",
        cfg.outdir.resolve(),
    )
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
