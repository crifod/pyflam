# SPDX-License-Identifier: AGPL-3.0-or-later
"""Available fuel load from the 10 m Tuscany FBFM40 map, aggregated to a forecast grid.

The pyroconvection product conditions on fire power, and fire power needs an available fuel
load ``w_a`` -- Tory & Kepert's appendix D firepower is ``alpha * h * w_a * dA/dt``. Until now
``w_a`` was a literature band (1.25-4 kg/m2, Australian eucalypt), three times wider than the
margins it was being used to judge.

This derives it from data instead: the 10 m FBFM40 map for Tuscany (CORINE 2018 + Italian
Forest Inventory, cross-walked against FirEUrisk; ``fuel-tos-10m``), converted per fuel model
to an available load and block-averaged onto the forecast grid.

Two choices worth stating, because both move the answer:

* **Which load.** The full grass/shrub load -- 1h + 10h + live herb + live woody -- not the
  fine dead and herbaceous components alone. In shrub models the woody and 10h classes carry
  most of the mass, and omitting them understates SH5 by a factor of ~2.4.
* **A +30 % Mediterranean-basin adjustment** on the model load (empirical, from operational
  use in the region; not part of Scott & Burgan). Applied as an upper bound on the tabulated
  value. It is a local calibration and should travel with that caveat.

The two are corroborated independently: the resulting area-weighted mean over Tuscany is
~1.4 kg/m2, inside the 1.31-3.64 kg/m2 range obtained by inverting Byram through the fireline
intensities Castellnou et al. (2022) observed in the Catalan campaign.

Outputs a two-band GeoTIFF on the target grid: mean available load over burnable cells
(kg/m2), and burnable fraction. Cells with no burnable fuel carry nan load and zero fraction --
distinct states, since "no fuel" and "unsampled" are different things.

Usage:  PYTHONPATH=src python scripts/fuel_load_10m.py [--res 2200] [--out PATH]
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import rasterio
from rasterio.transform import Affine

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from pyflam import fuel_models

FM40_10M = os.environ.get(
    "PYFLAM_FUEL10M",
    "/Users/cristianofoderi/DATI/FUEL_TOS/fuel-tos-10m/data/outputs/FuelModel_FM40_Toscana_10m.tif")
LB_FT2_TO_KG_M2 = 4.882
MED_ADJUST = float(os.environ.get("PYFLAM_FUEL_MED_ADJUST", 1.30))


def load_table(med_adjust: float = MED_ADJUST) -> dict:
    """FBFM40 code -> available fuel load (kg/m2). Non-burnable models are absent."""
    out = {}
    for code in range(1, 256):
        try:
            f = fuel_models.get(code)
        except KeyError:
            continue
        if not f.is_burnable:
            continue
        w = (f.load_1h + f.load_10h + f.load_live_herb + f.load_live_woody) * LB_FT2_TO_KG_M2
        out[code] = w * med_adjust
    return out


def aggregate(src_path: str, res_m: float):
    """Block-average the 10 m load onto a ``res_m`` grid. -> (load, burn_frac, transform, crs)."""
    tbl = load_table()
    lut = np.full(256, np.nan, np.float32)
    for k, v in tbl.items():
        lut[k] = v

    with rasterio.open(src_path) as d:
        f = int(round(res_m / d.res[0]))
        if f < 1:
            raise ValueError(f"target {res_m} m is finer than the source {d.res[0]} m")
        ny, nx = d.height // f, d.width // f
        load = np.full((ny, nx), np.nan, np.float32)
        frac = np.zeros((ny, nx), np.float32)
        for j in range(ny):                       # row-block streaming: the source is ~535 Mcell
            win = rasterio.windows.Window(0, j * f, nx * f, f)
            blk = lut[d.read(1, window=win)]
            blk = blk.reshape(f, nx, f)
            with np.errstate(invalid="ignore"):
                load[j] = np.nanmean(blk, axis=(0, 2))
            frac[j] = np.isfinite(blk).mean(axis=(0, 2))
        tr = d.transform * Affine.scale(f, f)
        return load, frac, tr, d.crs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", type=float, default=2200.0, help="target cell size in metres")
    ap.add_argument("--out", default="docs/fuel_load_tuscany.tif")
    a = ap.parse_args()

    load, frac, tr, crs = aggregate(FM40_10M, a.res)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with rasterio.open(a.out, "w", driver="GTiff", height=load.shape[0], width=load.shape[1],
                       count=2, dtype="float32", crs=crs, transform=tr, nodata=np.nan,
                       compress="deflate") as o:
        o.write(load, 1); o.set_band_description(1, "available_fuel_load_kg_m2")
        o.write(frac, 2); o.set_band_description(2, "burnable_fraction")

    ok = np.isfinite(load)
    print(f"wrote {a.out}  {load.shape[1]}x{load.shape[0]} @ {a.res:.0f} m")
    print(f"  burnable fraction: mean {frac.mean():.2f}, cells with any fuel {100*ok.mean():.1f}%")
    print(f"  available load over burnable cells: p10 {np.nanpercentile(load,10):.2f}  "
          f"median {np.nanmedian(load):.2f}  p90 {np.nanpercentile(load,90):.2f} kg/m2")
    print(f"  (+{100*(MED_ADJUST-1):.0f}% Mediterranean adjustment applied)")


if __name__ == "__main__":
    main()
