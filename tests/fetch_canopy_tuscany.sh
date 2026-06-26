#!/bin/sh
# Fetch GEDI-derived canopy height + cover for Tuscany and warp them onto the
# fuel-data grid (EPSG:3035, 100 m, 2166x2470), as inputs for deriving canopy base
# height / bulk density (see build_canopy_landscape_tuscany.py).
#
# Source: Meta / WRI High-Resolution Canopy Height (Tolan et al. 2024), GEDI-
# calibrated, on the open AWS bucket dataforgood-fb-data (no login). The 10-degree
# EPSG:4326 tiles covering Tuscany (40-50N, 0-20E) are lat=50.0_lon={0,10}.
#
# Usage:  sh tests/fetch_canopy_tuscany.sh [outdir]
set -e
OUT="${1:-/tmp/canopy}"
mkdir -p "$OUT"
BASE="https://dataforgood-fb-data.s3.amazonaws.com/forests/v1/alsgedi_global_v6_float_epsg4326_v3_10deg"

export GDAL_DISABLE_READDIR_ON_OPEN=EMPTY_DIR
export CPL_VSIL_CURL_ALLOWED_EXTENSIONS=.tif
export VSI_CACHE=TRUE GDAL_HTTP_MULTIRANGE=YES

# Tuscany extent in EPSG:3035 (bounds of Fuel_Model_40_tos_epsg3035.tif).
for metric in avg cover; do
    echo "warping ${metric} -> ${OUT}/${metric}.tif ..."
    gdalwarp -t_srs EPSG:3035 \
        -te 4295974.6022546 2126578.20606219 4512557.85719512 2373550.69415995 \
        -ts 2166 2470 -r average -overwrite -q \
        "/vsicurl/${BASE}/meta_chm_lat=50.0_lon=0.0_${metric}.tif" \
        "/vsicurl/${BASE}/meta_chm_lat=50.0_lon=10.0_${metric}.tif" \
        "${OUT}/${metric}.tif"
done
echo "done -> ${OUT}/{avg,cover}.tif  (canopy height is cm; build_canopy_landscape_tuscany.py rescales)"
