#!/usr/bin/env python
"""Add canopy base height + bulk density to the Tuscany landscape from GEDI-derived
canopy height, so the crown-fire path can run on real data.

There is no global canopy-base-height (CBH) / canopy-bulk-density (CBD) product:
GEDI and the GEDI-calibrated canopy-height maps (here Meta/WRI High-Resolution
Canopy Height, Tolan et al. 2024, on the open AWS bucket dataforgood-fb-data) give
canopy *height* only. CBH and CBD are therefore **derived** from canopy height +
canopy cover with transparent fire-science heuristics -- the standard approach when
only height/cover are available (cf. Cruz, Alexander & Wakimoto 2003 on canopy fuel
stratum characteristics). They are estimates, not field measurements.

Run `tests/fetch_canopy_tuscany.sh` (gdalwarp of the Meta tiles to the Tuscany grid)
first, then this script.

Usage:
    PYTHONPATH=src python tests/build_canopy_landscape_tuscany.py \
        --height /tmp/canopy/avg.tif --out /Users/.../canopy_tuscany.lcp
"""

from __future__ import annotations

import argparse
import os

import numpy as np

import pyflam
from pyflam import fuel_models
from pyflam.units import mph_to_ft_per_min

DATA = "/Users/cristianofoderi/DATI/FUEL_TOS/Toscana_Fuel_data"
GEOTIFFS = {
    "fuel_model": f"{DATA}/Fuel_Model_40_tos_epsg3035.tif",
    "slope": f"{DATA}/slope_deg_3035.tif",
    "aspect": f"{DATA}/aspect_deg_azimut_3035.tif",
    "elevation": f"{DATA}/dem_3035_sm.tif",
    "canopy_cover": f"{DATA}/cancov.tif",
}

# Heuristic constants (documented; tune for the forest type).
LOAD_MAX = 2.0          # max available canopy fuel load (kg/m^2) at full cover
CR_OPEN, CR_DENSE = 0.40, 0.80   # crown ratio at 0% and 100% cover (fuller crown ->
#                                  lower crown base -> more crown-prone)
CBD_MIN, CBD_MAX = 0.01, 0.40    # kg/m^3 clamp
MIN_CANOPY_HEIGHT = 2.0          # m; below this a cell is treated as non-forest


def derive_canopy(height_m, cover_pct):
    """CBH (m), CBD (kg/m^3), canopy_height (m) from canopy height + cover.

    Crown ratio rises with cover (denser stands carry crown fuel lower), so the
    crown base CBH = CH*(1-CR) falls with cover. CBD = canopy fuel load / canopy
    depth, with load scaled by cover. Cells with little canopy -> no crown fuel.
    """
    ch = np.maximum(height_m, 0.0)
    cc = np.clip(cover_pct, 0.0, 100.0) / 100.0
    forest = (ch >= MIN_CANOPY_HEIGHT) & (cc > 0.05)

    cr = CR_OPEN + (CR_DENSE - CR_OPEN) * cc          # crown ratio 0.4..0.8
    cbh = np.where(forest, np.maximum(ch * (1.0 - cr), 0.3), 0.0)
    depth = np.maximum(ch - cbh, 1.0)
    load = LOAD_MAX * cc                              # kg/m^2
    cbd = np.where(forest, np.clip(load / depth, CBD_MIN, CBD_MAX), 0.0)
    return cbh, cbd, np.where(forest, ch, 0.0)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--height", default="/tmp/canopy/avg.tif",
                    help="Meta canopy-height raster warped to the Tuscany grid")
    ap.add_argument("--out", default="/Users/cristianofoderi/DATI/FUEL_TOS/"
                    "pyflam_canopy_tuscany/canopy_tuscany.lcp")
    args = ap.parse_args(argv)
    outdir = os.path.dirname(args.out)
    os.makedirs(outdir, exist_ok=True)

    import rasterio
    ls = pyflam.Landscape.from_geotiffs(GEOTIFFS)
    for band in ("slope", "aspect", "elevation", "canopy_cover"):
        a = getattr(ls, band, None)
        if a is not None:
            setattr(ls, band, np.nan_to_num(np.asarray(a, dtype=float), nan=0.0))

    with rasterio.open(args.height) as d:
        h = d.read(1).astype(float)
        nd = d.nodata
    if nd is not None:
        h[h == nd] = 0.0
    h[~np.isfinite(h)] = 0.0
    # Infer units: heights are metres (0-60); the Meta uint16 avg may be scaled.
    hmax = np.nanpercentile(h[h > 0], 99) if np.any(h > 0) else 0.0
    scale = 1.0
    if hmax > 600:
        scale = 0.01     # cm -> m
    elif hmax > 100:
        scale = 0.1      # dm -> m
    height_m = h * scale
    print(f"canopy height: 99th pct raw {hmax:.0f} -> scale {scale} -> "
          f"{np.nanpercentile(height_m[height_m>0], 99):.1f} m, "
          f"forest cover {100*(height_m>=MIN_CANOPY_HEIGHT).mean():.1f}% of grid")

    cover = np.asarray(ls.canopy_cover, dtype=float)
    cbh, cbd, ch = derive_canopy(height_m, cover)

    # Store in LCP integer units: CBH/CH = m*10, CBD = (kg/m^3)*100.
    ls.canopy_height = np.rint(ch * 10).astype(np.int32)
    ls.canopy_base_height = np.rint(cbh * 10).astype(np.int32)
    ls.canopy_bulk_density = np.rint(cbd * 100).astype(np.int32)

    for name, arr in (("canopy_height_m", ch), ("canopy_base_height_m", cbh),
                      ("canopy_bulk_density_kg_m3", cbd)):
        ls.to_geotiff(os.path.join(outdir, f"{name}.tif"), arr)
    ls.to_lcp(args.out)
    print(f"wrote {args.out} (full canopy stack) + derived band rasters in {outdir}")

    # Verify the crown path now runs on the real landscape.
    burnable = np.zeros(ls.shape, bool)
    for num in np.unique(ls.fuel_model):
        try:
            if fuel_models.get(int(num)).is_burnable:
                burnable |= ls.fuel_model == int(num)
        except KeyError:
            pass
    sc = dict(m_1h=0.05, m_10h=0.06, m_100h=0.07, m_live_herb=0.6, m_live_woody=0.9)
    for model in ("cruz2005", "rothermel1991"):
        out = pyflam.crownfire.crown_fire_potential(
            ls, foliar_moisture=100.0, wind_20ft_ft_per_min=mph_to_ft_per_min(30),
            wind_midflame=pyflam.midflame_field(ls, mph_to_ft_per_min(30)),
            crown_spread=model, **sc)
        ft = out["fire_type"][burnable]
        d = {n: int((ft == t).sum()) for t, n in enumerate(("surface", "passive", "active"))}
        print(f"  crown types ({model:13s}) on real Tuscany: {d}")


if __name__ == "__main__":
    main()
