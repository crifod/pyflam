# pyflam GUI suite

A Streamlit multipage app over the `pyflam` fire-behavior pipelines. It is an
orchestration + presentation layer — no new science — so every page maps onto an
existing `pyflam` entry point.

## Install & run

```bash
# from the repo root
pip install -e ".[gui,geo,atmos]"     # gui = Streamlit; geo = rasterio/pyproj; atmos = weather I/O
streamlit run pyflam_gui/Home.py
```

`pyflam_gui/core/__init__.py` adds the sibling `src/` to `sys.path`, so the app
also runs from a source checkout without an editable install.

## Pages

| Page | Pipeline | Weather sources |
|------|----------|-----------------|
| **Home** | shared setup: area of interest (draw or preset), landscape (.lcp / GeoTIFF / synthetic), output folder | — |
| **1 · Fire-Weather Preview** | `meteo_report` + `pyroconvection_type` (potential & fuel-gated) | any (ICON-2I for the maps) |
| **2 · Realtime Propagation** | `run_realtime` + crown/pyroconvection/plume toggles, energy-flux dynamics, operative analysis, interactive perimeter map | GFS / ICON-2I (present) |
| **3 · Event Reanalysis** | same engine as page 2 at a historical date | ERA5 (past) |
| **4 · Burn Probability** | `burn_probability` with weather **ensembles** + ember **spotting** (FLP, conditional flame length/intensity, fire-size distribution) | ERA5 (past) or GFS / ICON (present) |

### Pages 2 & 3 (shared panel)

Ignitions are placed by clicking the map (or entered manually / at the AOI centre,
multiple points supported). The **run-duration scenario** (start time + total time
+ `dt`) is declared before the run. Toggles: crown fire (with foliar moisture),
pyroconvection feedback, CFD plume, and per-cell (spatial) weather forcing. Results
are tabbed: arrival + interactive perimeter maps, energy-flux time series, operative
driving-force analysis, and a download bundle (arrival/burned GeoTIFF, operative
GeoJSON, meteo CSV, PDF report). *Note: crown fire and per-cell spatial forcing are
mutually exclusive in this pyflam version — the panel falls back to single-column
when both are selected; the energy-flux series is sampled at the AOI centre either way.*

### Page 4 (burn probability)

Ignitions are random (count + seed) or drawn on the map. The context weather can be
a single state or a **weather ensemble** — add rows (date / hour / weight) to sample
the chosen source at several times so fire-to-fire variation drives a meaningful
flame-length distribution. Optional ember spotting (`SpottingModel`) lets fires
cross fuel barriers. Outputs: burn-probability + conditional flame-length/intensity
GeoTIFFs and a per-fire size CSV.

The run-duration scenario (start time + total time) is declared **before** each
propagation run. Every run writes a timestamped folder under the chosen output
directory containing GeoTIFF rasters, GeoJSON vectors, PNG figures and (if pandoc
is available) a PDF report, downloadable as a single `.zip`.

## Offline smoke test (no network/credentials)

On the Home page pick the **Tuscany** preset, build a **Synthetic** landscape,
then on **Fire-Weather Preview** choose the **Constant / idealized** source and
build the meteo report — charts render with no data files. The ICON-2I
pyroconvection maps, GFS/ERA5 propagation and burn probability need the
corresponding optional dependencies and (for ERA5) CDS credentials.

## Shared code

`pyflam_gui/core/` holds the reusable pieces: `aoi` (map + bbox helpers),
`landscape`, `atmosphere` (provider picker), `pyroconv` (classification + fuel
gate, shared with `tests/pyroconv_daily.py`), `outputs` (run-dir + writers) and
`plotting`.
