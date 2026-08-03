# Montale 2017 — landscape (LCP) build report

**Grid:** 10 m, EPSG:3035, 600x600 cells (6x6 km), centred on the
ignition (43.96095 N, 11.04026 E; 3035 X=4404680 Y=2317385).
**Extent (3035):** W 4401680 .. E 4407680, S 2314385 .. N 2320385.

## Layers
| band | source | note |
|---|---|---|
| fuel_model (FBFM40) | fuel-tos-10m FuelModel_FM40 (3035, 10 m) | window read |
| elevation | PT 10 m DTM (EPSG:3003) | reprojected to 3035 |
| slope / aspect | derived from the DEM | gradient |
| canopy cover/height/CBH/CBD | canopy_tuscany.lcp (100 m) | resampled to 10 m |

## Statistics
- elevation: 69–1041 m (median 427)
- slope: median 35 %, max 182 %
- fuel models present (FBFM40 codes): [91, 102, 104, 105, 123, 124, 142, 145, 147, 161, 162, 163, 165, 186, 189, 201]
- burnable cells: 360000/360000 (100 %)
- canopy cover: median 80 %, height median 73 (LCP units)

Output: `montale_10m.lcp`
