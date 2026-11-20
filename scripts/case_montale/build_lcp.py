"""Montale 2017 case — build a 10 m LCP for the Tobbiana (Montale, PT) fire.

Target grid: 10 m, EPSG:3035 (the fuel CRS), 6x6 km centred on the ignition
(43.96095 N, 11.04026 E -> 3035 X=4404680 Y=2317385). Layers:
  fuel_model  FBFM40 10 m (fuel-tos-10m, EPSG:3035)   window read
  elevation   PT 10 m DTM (EPSG:3003)                 reprojected to 3035
  slope/aspect derived from the reprojected DEM
  canopy      canopy_tuscany.lcp (100 m, EPSG:3035)   resampled to 10 m
"""
import os, sys, numpy as np, warnings
warnings.simplefilter("ignore")
import rasterio
from rasterio.windows import from_bounds
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin, Affine
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from pyflam.landscape import Landscape

FUEL = "/Users/cristianofoderi/DATI/FUEL_TOS/fuel-tos-10m/data/outputs/FuelModel_FM40_Toscana_10m.tif"
DTM = "/Users/cristianofoderi/Dropbox/BKgen19/Desktop/incendio_Montale/PT_dtm10x10_epsg3003_md5_35203e3a75bc0b7701b6a351dc310334/PT.asc"
CANOPY = "/Users/cristianofoderi/DATI/FUEL_TOS/pyflam_canopy_tuscany/canopy_tuscany.lcp"
OUT = os.path.join(os.path.dirname(__file__), "..", "..", "docs", "case_montale_2017", "landscape")

CX, CY = 4404680.0, 2317385.0            # ignition in EPSG:3035
HALF, RES = 3000.0, 10.0                  # 6x6 km, 10 m
N = int(2 * HALF / RES)                    # 600
WEST, NORTH = CX - HALF, CY + HALF
DST = from_origin(WEST, NORTH, RES, RES)   # target transform (3035)
SHAPE = (N, N)


def main():
    os.makedirs(OUT, exist_ok=True)

    # --- fuel: window read (already 3035 10 m) ---
    with rasterio.open(FUEL) as src:
        win = from_bounds(WEST, NORTH - 2*HALF, WEST + 2*HALF, NORTH, src.transform)
        fuel = src.read(1, window=win, out_shape=SHAPE, resampling=Resampling.nearest).astype("int16")

    # --- elevation: reproject DTM 3003 -> 3035 grid (explicit array + src_crs/nodata) ---
    elev = np.full(SHAPE, np.nan, "float32")
    with rasterio.open(DTM) as src:
        a = src.read(1).astype("float32")
        reproject(a, elev, src_transform=src.transform, src_crs="EPSG:3003",
                  src_nodata=src.nodata, dst_transform=DST, dst_crs="EPSG:3035",
                  dst_nodata=np.nan, resampling=Resampling.bilinear)
    med = np.nanmedian(elev[np.isfinite(elev)])
    elev = np.where(np.isfinite(elev), elev, med)

    # --- slope (%) + aspect (deg) from the DEM ---
    gy, gx = np.gradient(elev.astype("float64"), RES, RES)      # gy south-positive (row down)
    slope_pct = np.hypot(gx, gy) * 100.0
    aspect = (np.degrees(np.arctan2(gx, gy)) + 360.0) % 360.0   # downslope azimuth

    # --- canopy: read the 100 m LCP, resample its bands to the 10 m grid ---
    can = Landscape.from_lcp(CANOPY)
    cw, cn, ccs = can.west, can.north, can.cellsize_x
    src_tr = from_origin(cw, cn, ccs, ccs)
    canopy = {}
    for name in ("canopy_cover", "canopy_height", "canopy_base_height", "canopy_bulk_density"):
        arr = getattr(can, name)
        if arr is None:
            canopy[name] = np.zeros(SHAPE, "int16"); continue
        dst = np.zeros(SHAPE, "float32")
        reproject(arr.astype("float32"), dst, src_transform=src_tr, src_crs="EPSG:3035",
                  dst_transform=DST, dst_crs="EPSG:3035", resampling=Resampling.bilinear)
        canopy[name] = np.rint(dst).astype("int16")

    ls = Landscape(
        fuel_model=fuel, slope=slope_pct.astype("float32"), elevation=elev.astype("float32"),
        aspect=aspect.astype("float32"),
        canopy_cover=canopy["canopy_cover"], canopy_height=canopy["canopy_height"],
        canopy_base_height=canopy["canopy_base_height"], canopy_bulk_density=canopy["canopy_bulk_density"],
        cellsize_x=RES, cellsize_y=RES, west=WEST, north=NORTH, crs="EPSG:3035")
    lcp_path = os.path.join(OUT, "montale_10m.lcp")
    ls.to_lcp(lcp_path, latitude=44)

    # --- report ---
    fm = fuel[fuel > 0]
    with open(os.path.join(OUT, "landscape_report.md"), "w") as f:
        f.write(f"""# Montale 2017 — landscape (LCP) build report

**Grid:** 10 m, EPSG:3035, {N}x{N} cells ({2*HALF/1000:.0f}x{2*HALF/1000:.0f} km), centred on the
ignition (43.96095 N, 11.04026 E; 3035 X={CX:.0f} Y={CY:.0f}).
**Extent (3035):** W {WEST:.0f} .. E {WEST+2*HALF:.0f}, S {NORTH-2*HALF:.0f} .. N {NORTH:.0f}.

## Layers
| band | source | note |
|---|---|---|
| fuel_model (FBFM40) | fuel-tos-10m FuelModel_FM40 (3035, 10 m) | window read |
| elevation | PT 10 m DTM (EPSG:3003) | reprojected to 3035 |
| slope / aspect | derived from the DEM | gradient |
| canopy cover/height/CBH/CBD | canopy_tuscany.lcp (100 m) | resampled to 10 m |

## Statistics
- elevation: {np.nanmin(elev):.0f}–{np.nanmax(elev):.0f} m (median {np.nanmedian(elev):.0f})
- slope: median {np.median(slope_pct):.0f} %, max {np.max(slope_pct):.0f} %
- fuel models present (FBFM40 codes): {sorted(np.unique(fm).tolist())}
- burnable cells: {int((fuel>0).sum())}/{fuel.size} ({100*(fuel>0).mean():.0f} %)
- canopy cover: median {np.median(canopy['canopy_cover']):.0f} %, height median {np.median(canopy['canopy_height']):.0f} (LCP units)

Output: `montale_10m.lcp`
""")
    print(f"OK LCP -> {lcp_path}")
    print(f"  fuel models: {sorted(np.unique(fm).tolist())}")
    print(f"  elev {np.nanmin(elev):.0f}-{np.nanmax(elev):.0f} m, slope med {np.median(slope_pct):.0f}%, "
          f"burnable {100*(fuel>0).mean():.0f}%")


if __name__ == "__main__":
    main()
