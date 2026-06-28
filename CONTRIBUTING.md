# Contributing to pyflam

Thanks for your interest in pyflam — an open, multiplatform reimplementation of
FlamMap's fire-behavior science, extended into a weather-driven near-real-time
wildfire tool. Contributions of all kinds are welcome: bug reports, fixes, new
fuel/atmosphere models, validation cases, documentation and the GUI.

By participating you agree to abide by our [Code of Conduct](CODE_OF_CONDUCT.md).

## Getting set up

pyflam uses a `src/` layout and optional-dependency groups. Install in editable
mode with the extras you need:

```bash
git clone https://github.com/crifod/pyflam.git
cd pyflam
python -m pip install -e ".[dev]"                 # core + pytest
# add extras as needed:
#   geo   -> rasterio/pyproj/scikit-image (landscape I/O, GeoJSON, perimeters)
#   atmos -> xarray/cfgrib/netcdf4/cdsapi (GFS/ICON-2I/ERA5 forcing)
#   accel -> numba (faster anisotropic-Eikonal sweep)
#   gui   -> streamlit + folium (the pyflam_gui app)
python -m pip install -e ".[dev,geo,atmos,gui]"   # everything
```

## Running the tests

```bash
python -m pytest                 # full suite (pythonpath=src is preconfigured)
python -m pytest tests/test_rothermel.py -q     # a single file
python -m pytest --cov=pyflam    # with coverage
```

Live/network tests (GFS via Herbie, ERA5 via the Copernicus CDS) are skipped
unless their dependencies and credentials are present; ERA5 also needs
`PYFLAM_LIVE_ERA5=1`. Please keep new unit tests offline and deterministic.

## Running the GUI

```bash
python -m pip install -e ".[gui,geo,atmos]"
streamlit run pyflam_gui/Home.py
```

## Coding guidelines

- **Match the surrounding code.** Mirror the existing naming, comment density and
  numpy-style docstrings; keep the science vectorized (compute per unique fuel
  model and apply across the grid, as in `rothermel`/`mtt`).
- **Keep `src/pyflam` framework-free** — no Streamlit/UI imports in the science
  library. GUI-only code lives in `pyflam_gui/`.
- **Units matter.** pyflam's native units are ft, ft/min, Btu/ft/s; SI is used at
  the atmosphere boundary. Document units on every public function.
- **Add a test** for any behavior change, and a short validation note when a
  result is checked against FlamMap or the literature.
- Make heavy/optional imports (rasterio, xarray, matplotlib, herbie, cdsapi)
  **lazy**, so the base install stays light.

## Pull requests

1. Branch off `main` (e.g. `fix/...`, `feat/...`).
2. Keep each PR focused; write a clear description of *what* and *why*.
3. Ensure `python -m pytest` passes and new code is covered.
4. Fill in the pull-request template. CI (lint + tests + coverage) must be green.

## Reporting bugs / requesting features

Open an issue using the templates under
[`.github/ISSUE_TEMPLATE`](.github/ISSUE_TEMPLATE). For security-sensitive
reports, follow [SECURITY.md](SECURITY.md) instead of opening a public issue.
