# SPDX-License-Identifier: AGPL-3.0-or-later
"""Per-fire available fuel load from the EFFIS European Fuel Map.

The capability margin needs ``w_a`` per fire (Tory & Kepert appendix D,
``FP = alpha * h * w_a * dA/dt``). Until now every fire carried the same 1.49 kg/m2 -- the
median of the Tuscany 10 m FBFM40 map (:mod:`scripts.fuel_load_10m`) -- applied to Catalan,
Greek, Dutch and Chilean fires alike. That is the single crudest term in the chain, and
docs/graf_vs_pyflam_2026-07-26.md sec. 21.10 charges ~2.7x of the fire-side gap to it.

The EFFIS European Fuel Map fixes it for the European fires. It classifies Europe at 250 m into
the **Anderson (1982) NFFL 13** models -- of which it uses 1-10 -- and pyflam already implements
those as :data:`pyflam.fuel_models.STANDARD_13`, so the crosswalk is a lookup rather than an
invented mapping.

Download (18 MB zipped, 622 MB as GeoTIFF, EPSG:3035):

    https://data.effis.emergency.copernicus.eu/effis/applications/data-and-services/FuelMap_LAEA.zip

Point ``PYFLAM_EFFIS_FUEL`` at the extracted ``FuelMap2000_NFFL_LAEA.tif``. The raster is far too
large to vendor; only the per-fire result is stored, in
``docs/pyroconv_validation/fuel_effis.json``.

Three choices worth stating, because each moves the answer:

* **Sampling radius.** The equivalent radius of the fire's own burnt area (floored at 2 km),
  not a fixed buffer -- the relevant fuel is what the fire actually ran through, and the portal
  publishes each fire's area.
* **Which load.** 1h + 10h + live herb + live woody, matching ``fuel_load_10m`` so the two
  sources are directly comparable. Anderson 13 carries no live-herb load at all, so its totals
  run lower than FBFM40's for equivalent vegetation.
* **The +30 % Mediterranean adjustment** is carried over from the FBFM40 work, where it is an
  operational calibration for the Mediterranean basin. It is applied **only** to fires inside
  that basin -- a Dutch heath fire gets the unadjusted value. Like the original, it travels with
  the caveat that it is a local calibration, not part of Anderson (1982).

Chilean fires fall outside the map and keep the default; they are reported as ``null``.

The **per-fire NFFL class histogram is stored alongside** each result, so any other fuel
assignment -- in particular a different treatment of the unclassified agricultural cells -- can
be re-tested from the stored JSON without the 622 MB raster.

That treatment is not a detail. It decides the headline result of sec. 21: Guissona, the only
fire in the corpus that clears its PFT, sits on ground EFFIS leaves 93 % unclassified. Calling
that cereal NFFL 1 (short grass) would drop its firepower from 342 to 103 GW against a 139 GW
threshold and the detection would vanish. **NFFL 3 (tall grass) is the correct assignment for
Catalan cereal at harvest**, and it keeps it. See sec. 21.13.

Usage:  PYTHONPATH=src:scripts python scripts/fuel_load_effis.py [--out PATH]
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
import valdata
from pyflam import fuel_models

EFFIS_TIF = os.environ.get(
    "PYFLAM_EFFIS_FUEL",
    os.path.expanduser("~/DATI/EFFIS/FuelMap2000_NFFL_LAEA.tif"))
LB_FT2_TO_KG_M2 = 4.882
MED_ADJUST = float(os.environ.get("PYFLAM_FUEL_MED_ADJUST", 1.30))
MIN_RADIUS_M = 2000.0
# NFFL model for EFFIS's unclassified (agricultural) cells. 3 = tall grass, confirmed as the
# correct assignment for Catalan cereal at harvest by C. Foderi (2026-07-29). It is applied to
# every unclassified cell, so it also covers Greek agriculture, where it is plausible but not
# separately confirmed -- those fires run 40-60 % unclassified, against 93-99 % in the Catalan
# cereal belt, so the assignment matters most exactly where it is validated.
AGRI_NFFL = 3
# Mediterranean basin, for the +30 % adjustment: Iberia through Greece, south of the Alps.
MED_BOX = (34.0, 46.0, -10.0, 30.0)          # south, north, west, east


def load_table(med_adjust: float) -> dict:
    """NFFL code -> available fuel load (kg/m2). Non-burnable models absent."""
    out = {}
    for code in range(1, 14):
        try:
            f = fuel_models.get(code)
        except KeyError:
            continue
        if not f.is_burnable:
            continue
        w = (f.load_1h + f.load_10h + f.load_live_herb + f.load_live_woody) * LB_FT2_TO_KG_M2
        out[code] = w * med_adjust
    return out


def in_mediterranean(lat: float, lon: float) -> bool:
    s, n, w, e = MED_BOX
    return s <= lat <= n and w <= lon <= e


def sample_fire(src, lat: float, lon: float, area_ha: float | None, transformer,
                agri_nffl: int = 0):
    """-> dict of w_a and coverage around one fire, or None if off-map.

    ``agri_nffl`` assigns unclassified (class 0) cells an NFFL model. EFFIS does not classify
    agricultural land at all, and several Catalan fires -- Guissona 92 % unclassified, Granyena
    96 %, Santa Coloma 81 % -- burned through the Lleida/Segarra cereal belt. Left at 0 the
    sampled ``w_a`` describes only the few shrub and timber patches, which is not what burned.
    """
    import rasterio

    x, y = transformer.transform(lon, lat)
    r = MIN_RADIUS_M
    if area_ha and area_ha > 0:
        r = max(r, math.sqrt(area_ha * 1.0e4 / math.pi))
    win = rasterio.windows.from_bounds(x - r, y - r, x + r, y + r, src.transform)
    win = win.round_offsets().round_lengths()
    if win.width <= 0 or win.height <= 0:
        return None
    try:
        blk = src.read(1, window=win, boundless=True, fill_value=0)
    except Exception:
        return None
    if blk.size == 0 or not (blk > 0).any():
        return None

    tbl = load_table(MED_ADJUST if in_mediterranean(lat, lon) else 1.0)
    lut = np.full(256, np.nan, np.float32)
    for k, v in tbl.items():
        lut[k] = v
    unclassified = float((blk == 0).mean())
    if agri_nffl:
        lut[0] = tbl.get(agri_nffl, np.nan)
    vals = lut[blk]
    ok = np.isfinite(vals)
    if not ok.any():
        return None
    u, c = np.unique(blk, return_counts=True)
    hist = {int(k): int(v) for k, v in zip(u, c)}
    return dict(w_a=float(np.nanmean(vals)), covered=float(ok.mean()),
                unclassified=unclassified, radius_m=float(r), n_cells=int(ok.sum()),
                nffl_histogram=hist)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(valdata.DATA, "fuel_effis.json"))
    ap.add_argument("--tif", default=EFFIS_TIF)
    ap.add_argument("--agri-nffl", type=int, default=AGRI_NFFL, metavar="N",
                    help=f"assign unclassified (agricultural) cells NFFL model N "
                         f"(default {AGRI_NFFL}, tall grass); 0 leaves them out of the mean")
    a = ap.parse_args()

    if not os.path.exists(a.tif):
        raise SystemExit(f"EFFIS fuel map not found at {a.tif}\n"
                         f"Download FuelMap_LAEA.zip (see this module's docstring) and set "
                         f"PYFLAM_EFFIS_FUEL to the extracted GeoTIFF.")
    import rasterio
    from pyproj import Transformer

    fires = valdata.load("portal_fires")
    src = rasterio.open(a.tif)
    tr = Transformer.from_crs("EPSG:4326", src.crs, always_xy=True)

    out = {}
    print(f"{'fire':26s} {'w_a':>6s} {'unclass':>9s} {'radius km':>9s} {'cells':>7s}  med")
    for f in sorted(fires, key=lambda x: x["slug"]):
        if f.get("lat") is None:
            continue
        res = sample_fire(src, f["lat"], f["lon"], f.get("area"), tr, a.agri_nffl)
        if res is None:
            out[f["slug"]] = None
            print(f"{f['slug'][:26]:26s} {'--':>6s}   (outside the EFFIS map)")
            continue
        med = in_mediterranean(f["lat"], f["lon"])
        out[f["slug"]] = dict(w_a_kg_m2=round(res["w_a"], 3),
                              covered_fraction=round(res["covered"], 3),
                              unclassified_fraction=round(res["unclassified"], 3),
                              radius_m=round(res["radius_m"]), n_cells=res["n_cells"],
                              mediterranean_adjust=med, agri_nffl=a.agri_nffl,
                              nffl_histogram=res["nffl_histogram"])
        flag = "  <- mostly unclassified" if res["unclassified"] > 0.5 else ""
        print(f"{f['slug'][:26]:26s} {res['w_a']:6.2f} {res['unclassified']:9.2f} "
              f"{res['radius_m']/1000:9.1f} {res['n_cells']:7d}  "
              f"{'yes' if med else 'no'}{flag}")

    vals = [v["w_a_kg_m2"] for v in out.values() if v]
    json.dump(out, open(a.out, "w"), indent=1)
    print(f"\nwrote {a.out}: {len(vals)} fires on-map, {len(out)-len(vals)} off-map")
    if vals:
        print(f"w_a over on-map fires: min {min(vals):.2f}  median {np.median(vals):.2f}  "
              f"max {max(vals):.2f} kg/m2   (Tuscany FBFM40 default was 1.49)")


if __name__ == "__main__":
    main()
