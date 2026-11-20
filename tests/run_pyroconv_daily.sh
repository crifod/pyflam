#!/bin/bash
# Daily Tuscany pyroconvection-type maps from ICON-2I (run by launchd/cron).
# Generates today's potential + gated panels, GeoTIFFs and PDF into docs/daily/,
# then prunes ICON-2I GRIB caches older than 3 days.
set -euo pipefail

REPO="/Users/cristianofoderi/-softEST/firelab-flammap6_install_0828_2025/pyflam"
PYTHON="${PYROCONV_PYTHON:-/opt/anaconda3/bin/python3}"
RUN="${PYROCONV_RUN:-0}"                       # ICON-2I 00Z run (ready by ~05Z UTC)
export PYTHONPATH="$REPO/src"
# launchd/cron give a minimal PATH; add the dirs holding python deps, pandoc, tectonic.
export PATH="/opt/anaconda3/bin:/opt/homebrew/bin:/usr/local/bin:/usr/bin:/bin:$PATH"
export PYROCONV_CACHE="${PYROCONV_CACHE:-/tmp/pyflam_icon2i/$(date -u +%Y%m%d)${RUN}}"

cd "$REPO"
mkdir -p "$REPO/docs/daily"
DATE="$(date -u +%Y-%m-%d)"
echo "[$(date)] pyroconv_daily $DATE ${RUN}Z -> docs/daily"
"$PYTHON" tests/pyroconv_daily.py "$DATE" "$RUN"

# prune GRIB caches older than 3 days to bound disk use (~600 MB/run)
find /tmp/pyflam_icon2i -mindepth 1 -maxdepth 1 -type d -mtime +3 -exec rm -rf {} + 2>/dev/null || true
echo "[$(date)] done"
