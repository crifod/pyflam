# pyflam

An open-source Python reimplementation of the fire-behavior science behind
[FlamMap](https://www.firelab.org/project/flammap) — built up from the published,
peer-reviewed models rather than by decompiling the Windows binary.

**Status: roadmap steps 1–5 implemented** — the Rothermel surface fire spread
model (the scientific core of FlamMap's "Basic Fire Behavior": rate of spread,
reaction intensity, fireline intensity, flame length), landscape I/O, two terrain
wind solvers, crown fire, and a directional Minimum Travel Time spread/perimeter
engine. Validation against real FlamMap output (step 6) is the remaining work.

## Why reimplement instead of decompile

FlamMap is a closed Windows/C++ GUI app, but the science it runs on is open and
documented in USDA Forest Service publications. Reimplementing from the equations
gives clean, testable, cross-platform code; decompiling would give unreadable
machine output and reproduce platform quirks. Validation is done by diffing our
outputs against real FlamMap/BehavePlus runs, not by reading their code.

## Install

```bash
cd pyflam
python -m pip install -e ".[dev]"
```

(Only dependency is NumPy. `pytest` comes with the `dev` extra.)

## Quick start

```python
import pyflam
from pyflam.units import mph_to_ft_per_min, ft_per_min_to_m_per_min

fm = pyflam.get_fuel_model(1)        # FM1, short grass (Anderson 1982)
# or by Scott & Burgan code / number:
fm = pyflam.get_fuel_model("GR2")    # == pyflam.get_fuel_model(102)

fb = pyflam.spread(
    fm,
    m_1h=0.06, m_10h=0.07, m_100h=0.08,   # dead fuel moistures (fractions)
    m_live_herb=0.60, m_live_woody=0.90,  # live fuel moistures
    wind_midflame=mph_to_ft_per_min(5),   # midflame wind speed -> ft/min
    slope=0.30,                           # rise/run (30%)
)

print(fb.rate_of_spread, "ft/min")
print(ft_per_min_to_m_per_min(fb.rate_of_spread), "m/min")
print(fb.flame_length, "ft")
```

All computation is in Rothermel's native English units; `pyflam.units` has the
converters to/from SI.

### Whole-landscape behavior (step 2)

Run the surface model over an entire landscape — the per-cell equivalent of
FlamMap's "Basic Fire Behavior" rasters:

```python
import pyflam
from pyflam.units import mph_to_ft_per_min

# From a FlamMap/FARSITE landscape file (pure Python, no GDAL needed):
ls = pyflam.Landscape.from_lcp("mylandscape.lcp")

# ...or from LANDFIRE GeoTIFFs (needs rasterio):
# ls = pyflam.Landscape.from_geotiffs({"fuel_model": "fbfm40.tif", "slope": "slpd.tif"})

out = pyflam.basic_fire_behavior(
    ls,
    m_1h=0.06, m_10h=0.07, m_100h=0.08,
    m_live_herb=0.60, m_live_woody=0.90,
    wind_midflame=mph_to_ft_per_min(5),
)
out["rate_of_spread"]      # 2D array, ft/min  (nonburnable cells = 0)
out["flame_length"]        # ft
out["fireline_intensity"]  # Btu/ft/s
out["reaction_intensity"]  # Btu/ft^2/min

ls.to_geotiff("ros.tif", out["rate_of_spread"])  # write a result raster
```

Fuel moisture is uniform; fuel model and slope vary per cell. `wind_midflame`
may be a single value **or a 2D field** (see the wind solvers below). Wind and
slope are combined as Rothermel's scalar `1 + phi_w + phi_s` (wind upslope /
maximum-spread case) — directional spread using aspect + wind bearing is planned.

## Terrain wind: two native solvers

Real terrain bends the wind (ridge speed-up, valley channeling, lee wakes), so a
single uniform value is a poor approximation. pyflam computes a gridded
`WindField` two ways, both feeding the same `wind_midflame` interface. (The
science follows [WindNinja](https://github.com/firelab/windninja), USFS — but
pyflam reimplements it natively rather than bridging to that tool.)

### Fast: mass-consistent (`pyflam.windsolver`) — no external deps

The Sasaki **mass-consistent** (diagnostic) model: conserves mass by a
finite-volume variational solve over a terrain-masked 3D grid (just scipy sparse
linear algebra). Seconds to run; divergence-free; self-checked.

```python
from pyflam import windsolver

wf = windsolver.wind_field_from_landscape(ls, speed=5.0, direction=270.0)
# roughness z0 from the fuel grid; elevation from the LCP.
midflame = wf.to_midflame(ls, wind_reduction_factor=0.4)
out = pyflam.basic_fire_behavior(ls, wind_midflame=midflame, m_1h=0.06,
                                 m_10h=0.07, m_100h=0.08,
                                 m_live_herb=0.60, m_live_woody=0.90)
```

### High-fidelity: momentum/RANS CFD (`pyflam.cfd`) — via OpenFOAM

True momentum conservation: a k-ε RANS atmospheric-boundary-layer simulation on a
terrain-following mesh, built/run/read through OpenFOAM (ESI). It captures what
the diagnostic model can't — flow separation and lee wakes — and **non-neutral
stability with diurnal slope flows** (Boussinesq buoyancy + a ground heat flux).

```python
from pyflam import cfd

# neutral:
wf = cfd.wind_field_from_landscape(ls, speed=6.0, direction=270.0)

# diurnal: daytime heating -> enhanced upslope; night cooling -> katabatic.
wf = cfd.solve_rans(dem, cellsize=30, speed=4.0, direction=270.0,
                    buoyant=True, surface_heat_flux=150.0)   # W/m^2 (+day/-night)

midflame = wf.to_midflame(ls, wind_reduction_factor=0.4)
out = pyflam.basic_fire_behavior(ls, wind_midflame=midflame, m_1h=0.06,
                                 m_10h=0.07, m_100h=0.08,
                                 m_live_herb=0.60, m_live_woody=0.90)
```

## Atmospheric forcing: real weather / reanalysis

`pyflam.atmosphere` drives a run from real atmospheric data instead of fixed
scenario inputs — a **forecast** (GFS / HRRR / WRF output) for near-real-time
runs, or a **reanalysis** (ERA5 from Copernicus) for re-analysis runs. It carries
the surface state *and* the convective / energy-flux fields that govern
fire–atmosphere coupling (surface heat flux, CAPE, CIN, PBL height, stability),
and derives pyflam inputs from them:

```python
import pyflam
from pyflam import atmosphere as atm

# 1. get an atmospheric state (here from a downloaded ERA5/GFS/WRF file)
prov = atm.open_atmosphere("era5_2024-08-01.nc", source="era5")   # needs xarray
state = prov.state_at(latitude=43.0, longitude=11.0, time=when)

# 2. derive the fire inputs (NFDRS equilibrium moisture from T/RH, midflame wind):
si = pyflam.spread_inputs_from_state(state)        # m_1h/10h/100h, wind, direction
field = pyflam.spread_field(ls, m_live_herb=0.6, m_live_woody=0.9, **si)

# 3. convection feeds the coupling: background heat flux + a plume enhancement
atm.ambient_surface_heat_flux(state)               # W/m^2 -> buoyant CFD background
atm.convective_plume_factor(state)                 # >1 unstable (CAPE) -> more spotting
fp = atm.atmospheric_firebrand_physics(state)      # firebrands reach farther
```

Run a **time-marching** simulation through evolving weather (the near-real-time /
reanalysis scenario) by handing the provider to the coupled march:

```python
res = pyflam.fire_atmosphere_march(
    ls, ignitions=[(r, c)], total_time=180, dt=15,
    atmosphere=prov, location=(43.0, 11.0), start_time=ignition_datetime,
    m_live_herb=0.6, m_live_woody=0.9,
)   # each 15-min step re-reads the wind, fuel moisture and convective heat flux
```

Three capabilities make this operational:

* **Spatial weather** — `spatial=True` (or `provider.field_on(ls, time)`) samples
  the atmosphere **per cell**, so wind and fuel moisture vary across the domain
  (`spread_field` and `basic_fire_behavior` accept per-cell moisture fields). With
  `spatial=True, plume=True` (`merge_plume_wind`) the fire's CFD plume perturbation
  is superposed onto the per-cell ambient wind — spatial weather *and* the fire's
  own plume together.
* **Live data fetch with caching** — `fetch_gfs(run, fxx)` (NOAA GFS via Herbie;
  open data, no auth — exercised by a live test) and `fetch_era5(cache_path,
  date, time, area)` (Copernicus CDS via `cdsapi`) download once and reuse; both
  return a `GriddedAtmosphere`. Longitudes are wrapped to each source's
  convention (GFS 0–360, ERA5 −180–180), and ERA5's accumulated J/m² surface
  fluxes are converted to instantaneous W/m² upward (`era5_flux_to_watts`).
* **Fuel-moisture memory** — `DeadFuelMoistureModel` is the operational time-lag
  (Nelson-type) model: the 1/10/100-h classes track the equilibrium moisture at
  their own response times, so fuels remember recent humidity instead of snapping
  to the latest value — important for diurnal drying and reanalysis runs.

The convective emphasis is deliberate: an unstable, high-CAPE atmosphere lets the
fire's plume rise higher (stronger indrafts, farther spotting), while a stable
layer caps it — `convective_plume_factor` and `ambient_surface_heat_flux` carry
that into the spotting and the buoyant-RANS coupling. Reading/fetching
forecast/reanalysis data needs the `atmos` extra (`pip install 'pyflam[atmos]'`);
the state physics and the synthetic/`ConstantAtmosphere` provider work with no
extra deps.

This is the step from "a robust FlamMap-equivalent core" toward a weather-driven
operational/research tool — real winds, real moisture, real convection.

### Operational products

```python
import pyflam
from datetime import datetime

# 1. fire-weather variation report over the operational window
report = pyflam.build_meteo_report(
    atmosphere=prov, location=(43.0, 11.0),
    start_time=datetime(2026, 8, 1, 11, 0), total_minutes=360, step_minutes=60)
print(report.summary())          # RH falling, fuels drying, atmosphere destabilising...

# 2. driving-force analysis of the first-run perimeter, by operational sector
arrival = pyflam.minimum_travel_time(field, ignitions=[(r, c)], max_time=90)
op = pyflam.analyze_perimeter(ls, arrival, time=90,
        wind_midflame=mid, wind_direction=270.0, subsectors=2, **sc)
print(op.summary())              # per sector: wind / slope / fuel vectors + resultant + driver
# op.quiver()                    # arrows along the perimeter (needs matplotlib)
```

For each sector (head / right flank / tail / left flank, or finer sub-sectors)
the report gives slope, fuel and wind as separate vectors (bearing + magnitude)
*and* their vector sum, so an analyst sees both the individual drivers and the net
push — e.g. a head running upslope, a flank carried by wind into heavier fuel.

Export it for a GIS / mapping front-end as **GeoJSON** — a contour-traced ordered
**Polygon** perimeter, sector-centroid points (force magnitudes/bearings, dominant
driver, ROS) and one LineString "arrow" per force, in the landscape's CRS or
reprojected to WGS84:

```python
op.write_geojson(ls, "operative.geojson", arrow_length=2000.0, to_wgs84=True)
```

### One-call near-real-time product

`pyflam.run_realtime` does the whole operational pipeline in a single call —
fetch/sample weather → spread through the evolving conditions → perimeter → both
reports:

```python
prod = pyflam.run_realtime(
    ls, ignitions=[(r, c)], atmosphere=prov,          # e.g. fetch_gfs(...)
    location=(43.0, 11.0), start_time=ignition_time,
    total_time=180, dt=30, m_live_herb=0.6, m_live_woody=0.9)
print(prod.summary())                                  # meteo + operative
prod.write_geojson(ls, "operative.geojson", to_wgs84=True)
prod.meteo.variation()                                 # the weather trend
```

Run the meteo report at several horizons to see the trend sharpen — e.g. 3 / 6 /
9 h ahead:

```python
for hours in (3, 6, 9):
    print(pyflam.build_meteo_report(prov, location=(43.0, 11.0),
          start_time=t0, total_minutes=hours * 60, step_minutes=60).summary())
```

### Pyroconvection: the fire reshapes its own wind

An intense fire is not a passive tracer — its convective heat builds a buoyant
plume that pulls air inward and lifts it over the front, locally bending and
strengthening the wind. `pyflam.pyroconvection` closes that loop through the
buoyant RANS solver:

```python
import pyflam

# spread the fire in the ambient wind to get per-cell fireline intensity
field = pyflam.spread_field(ls, wind_midflame=mid, wind_direction=270.0, **sc)

# fire energy -> ground heat-flux field -> buoyant CFD -> plume-modified wind
wf = pyflam.couple_fire_wind(
    ls, field.fireline_intensity,        # Btu/ft/s, per cell
    speed=3.0, direction=270.0,          # ambient inflow (m/s, deg FROM)
    active_mask=pyflam.perimeter_mask(arrival, time=30.0),  # only the burning area
)
# feed the plume-modified wind back into the next increment of spread:
field2 = pyflam.spread_field(ls, wind_midflame=wf.to_midflame(ls),
                             wind_direction=270.0, **sc)
```

`fire_heat_flux` is the physical bridge: Byram intensity `I` (energy per unit
length of front) becomes a cell-averaged convective flux `q = χ_c·I/cellsize`
(W/m²). The solve is quasi-steady (a steady plume for a frozen fire state) and
needs OpenFOAM.

OpenFOAM is an external engine, discovered at runtime via the `openfoam` wrapper
(`PYFLAM_OPENFOAM` or on `PATH`); install e.g. `brew install
gerlero/openfoam/openfoam` on macOS. pyflam writes the mesh and case itself (no
PyFoam dependency) and raises a clear error if OpenFOAM is absent. Caveats:
staircase terrain, a single representative `z0`, and *quasi-steady* diurnal
snapshots; buoyant RANS over steep terrain can be hard to converge. See
`tests/CFD_VALIDATION.md` for the benchmark harness (checkMesh, Richards & Hoxey,
Askervein Hill, diurnal slope flows).

## Crown fire (step 4)

Once the surface fire is intense enough it climbs into the canopy. `pyflam.crownfire`
links the surface model to the canonical crown-fire models — Van Wagner (1977)
initiation, Rothermel (1991) active spread, and Scott & Reinhardt (2001) for the
surface/passive/active classification, crown fraction burned, and the torching
and crowning indices. This module works in **SI** (canopy data is metric): canopy
base height in m, bulk density in kg/m³, foliar moisture %, spread rate m/min,
intensity kW/m.

```python
import pyflam
from pyflam.units import mph_to_ft_per_min

sc = dict(m_1h=0.06, m_10h=0.07, m_100h=0.08, m_live_herb=0.60, m_live_woody=0.90)
surface = pyflam.spread(pyflam.get_fuel_model(10),
                        wind_midflame=mph_to_ft_per_min(0.4 * 25), **sc)

cf = pyflam.crown_fire_behavior(
    surface,
    canopy_base_height=1.5,      # m
    canopy_bulk_density=0.15,    # kg/m^3
    foliar_moisture=100,         # %
    wind_20ft_ft_per_min=mph_to_ft_per_min(25),
    canopy_fuel_load=1.0,        # kg/m^2 (for crown intensity)
    **sc,
)
cf.fire_type                     # "surface" | "passive" | "active"
cf.crown_fraction_burned         # 0..1
cf.rate_of_spread                # m/min, final (Scott & Reinhardt blend)

# Wind thresholds (20-ft wind, mph) for torching and active crowning:
fm = pyflam.get_fuel_model(10)
pyflam.torching_index(fm, canopy_base_height=1.5, foliar_moisture=100, **sc)
pyflam.crowning_index(0.15, **sc)

# ...or per cell over a landscape (needs CBH + CBD bands):
out = pyflam.crown_fire_potential(ls, foliar_moisture=100,
        wind_20ft_ft_per_min=mph_to_ft_per_min(25),
        wind_midflame=mph_to_ft_per_min(0.4 * 25), **sc)
out["fire_type"]                 # 0 surface, 1 passive, 2 active
```

## Fire growth: directional spread + Minimum Travel Time (steps 3 & 5)

`pyflam.mtt` is the spread/perimeter engine. It first turns the scalar
`1 + phi_w + phi_s` into a **directional elliptical spread** template (step 3):
the wind and slope factors are combined as *vectors* (wind downwind, slope
upslope from the landscape aspect) to give each cell a maximum spread rate, a
heading, and an ellipse eccentricity (from the Anderson 1983 length-to-breadth
ratio). It then runs Finney's (2002) **Minimum Travel Time** algorithm — fire
arrival time as the shortest-time path from the ignition over a lattice of travel
directions, the same idea as FlamMap's MTT.

```python
import pyflam
from pyflam.units import mph_to_ft_per_min

sc = dict(m_1h=0.06, m_10h=0.07, m_100h=0.08, m_live_herb=0.60, m_live_woody=0.90)

res = pyflam.spread_perimeter(
    ls,
    ignitions=[pyflam.ignition_from_xy(ls, x0, y0)],   # or [(row, col), ...]
    wind_midflame=mph_to_ft_per_min(6),
    wind_direction=270.0,             # degrees the wind blows FROM (met)
    **sc,
)
arrival = res["arrival_time"]         # 2D array, minutes (inf = unburned)

burned = pyflam.perimeter_mask(arrival, time=60.0)     # burned area at 60 min
ls.to_geotiff("arrival.tif", arrival)

# The intermediate directional template is reusable on its own:
field = pyflam.spread_field(ls, wind_midflame=mph_to_ft_per_min(6),
                            wind_direction=270.0, **sc)
field.ros_max, field.heading, field.eccentricity       # per-cell ellipse
field.directional_ros(azimuth=45.0)                    # ROS toward NE, per cell
```

Per-cell wind from the terrain solvers feeds straight in: pass the
`windsolver`/`cfd` midflame field as `wind_midflame` and the matching direction
field as `wind_direction`. A `ring` argument trades angular resolution for speed.

**Scale.** The travel-time graph is assembled with vectorized NumPy and the
shortest path is solved by SciPy's C-level multi-source Dijkstra, so MTT runs on
large / high-resolution (meter-scale) landscapes: ~1M cells in under a second,
~4M cells in a few seconds on a laptop. `max_time` bounds the search and makes
fixed-duration runs cheaper. The result is identical to a reference heap Dijkstra
(`pyflam.mtt._mtt_python`, kept for cross-checking); `build_traveltime_graph`
exposes the sparse graph itself.

For huge landscapes the graph build is **chunked** to keep memory near the size
of the final graph. It is a two-pass, direct-CSR assembly: edges are counted
first (an edge exists iff both endpoints burn — needs only the burnable mask),
the CSR arrays are allocated once, then filled band-by-band with no COO
intermediate and no concatenation doubling. `chunk_rows` sets the band height
(`None` auto-chunks above a few million cells; `0` forces single-shot). On the
included 5.35M-cell Tuscany landscape (2470×2166 at 100 m, 24.6M edges) the build
drops from ~690 MB single-shot to ~470 MB chunked, and a full MTT solve from one
ignition takes ~2.6 s. *Note:* the graph itself still scales with cell count, so
a whole region at 1 m would need spatial tiling — chunking bounds the build, not
the final graph size.

## What's implemented

- `pyflam.rothermel` — full Rothermel (1972) surface model with Albini (1976)
  size-class net-load weighting, live fuel moisture of extinction, dynamic
  herbaceous load curing (Scott & Burgan 2005), wind and slope factors, Byram
  fireline intensity and flame length.
- `pyflam.fuel_models` — both standard sets, looked up by number or code:
  - the original 13 (Anderson 1982), numbers `1`–`13`;
  - the Scott & Burgan (2005) set, numbers `91`–`204` / codes like `GR1`, `SH5`,
    `TL3`, including the 40 burnable models and the 5 nonburnable (`NB`) models.
- `pyflam.landscape` — the `Landscape` raster stack and `basic_fire_behavior`,
  which runs the surface model over every cell (vectorized; the expensive terms
  are computed once per unique fuel model).
- `pyflam.io_lcp` — pure-Python reader/writer for FlamMap/FARSITE `.lcp`
  landscape files (validated against GDAL's LCP driver in the test suite).
- GeoTIFF I/O (`Landscape.from_geotiffs` / `Landscape.to_geotiff`) via rasterio,
  an optional dependency: `pip install 'pyflam[geo]'`.
- `pyflam.wind` — the shared `WindField` type plus a generic ESRI ASCII reader.
- `pyflam.windsolver` — native mass-consistent (Sasaki) wind model: a
  finite-volume, divergence-free terrain-wind solver in pure scipy (fast option).
- `pyflam.cfd` — momentum/RANS terrain-wind solver via OpenFOAM (k-ε ABL,
  terrain-following mesh, Boussinesq buoyancy + diurnal slope flows). Builds,
  runs and reads the OpenFOAM case from Python. Both the ground heat flux and the
  surface roughness `z0` may be a scalar **or a per-cell field** — `z0` defaults
  to per-cell roughness from the fuel grid (the ground wall functions vary cell by
  cell; the ABL inlet uses the median).
- `pyflam.pyroconvection` — fire→atmosphere coupling: turns the fire's fireline
  intensity into a convective ground heat-flux field, runs the buoyant RANS, and
  returns the wind including the fire's own plume (indrafts/updraft) — the
  feedback where the fire reshapes the wind that drives it. `couple_fire_wind` is
  the one-shot coupling; `fire_atmosphere_march` time-marches it, re-solving the
  plume wind every `dt` minutes of MTT growth (wind solver injectable, so the loop
  is usable/testable with or without OpenFOAM). Quasi-steady.
- `pyflam.rothermel.SurfaceKernel` — the wind/slope-independent Rothermel terms,
  computed once per fuel + moisture and applied to scalar or array (wind, slope).
- `pyflam.crownfire` — crown fire: Van Wagner (1977) initiation, Rothermel (1991)
  active spread, Scott & Reinhardt (2001) surface/passive/active classification,
  crown fraction burned, and torching/crowning indices (point + landscape).
- `pyflam.mtt` — fire growth: directional elliptical spread (Finney 1998 +
  Anderson 1983 length-to-breadth), the Minimum Travel Time arrival-time /
  perimeter solver (Finney 2002), `burn_probability` (many fixed-duration fires on
  one prebuilt graph), and `spread_with_spotting` (MTT coupled to ember spotting).
- `pyflam.spotting` — ember (firebrand) spotting, two models that both couple
  into MTT growth and burn probability:
  - `SpottingModel` — fast **parameterized** model (flame-length loft + lognormal
    landing), easy to calibrate.
  - `FirebrandPhysics` — **stochastic, physics-based** model where spotting
    *emerges* from the energy system: buoyant-plume loft from fireline intensity
    (Morton-Taylor-Turner), firebrand size → terminal velocity (drag) → loft and
    combustion burnout, downwind transport, and a landing-ignition probability
    that falls with the receiving fuel's moisture. Randomness (Poisson brand
    counts ∝ intensity, lognormal sizes, turbulent bearing, Bernoulli ignition)
    makes the landing pattern a Monte-Carlo outcome; constants are physical, not
    reach calibrations. The particle sub-models are tied to **measured firebrand
    data**: terminal velocity reproduces measured settling speeds (Tohidi & Kaye;
    Manzello), burnout follows Tarifa's d²-law constant (`calibrate_burnout`,
    `from_burning_constant`), and the size distribution uses Manzello's firebrand
    sizes. Its one weakly-constrained length (`front_length`) is calibrated against
    literature spot-distance anchors (`calibrate_front_length`,
    `spot_distance_report`); pass your own measured `(intensity, wind, distance)`
    rows to recalibrate. `spot_distance_distribution` exposes the per-source
    landing distribution for direct comparison with data.
- **Fuel-load factor** — `load_factor` on `spread` / `surface_kernel` /
  `basic_fire_behavior` / `spread_field` scales the standard models' loads (they
  under-represent real loading by ~20-30%); flows through the full Rothermel
  computation, raising fireline intensity (hence flame length and spotting reach).
  Accepts a scalar, a per-fuel `{number: factor}` map, or a per-cell raster (e.g.
  a LANDFIRE-derived load-correction layer). ROS is non-monotonic in load —
  compact litter can slow past the optimum packing ratio.
- `pyflam.wind_reduction` — wind adjustment factor (20-ft → midflame): the
  unsheltered, fuel-depth-based and sheltered, canopy-based forms (Albini &
  Baughman 1979 / Andrews 2012), per cell over a landscape.
- `pyflam.atmosphere` — atmospheric forcing from forecast/reanalysis data
  (GFS/HRRR/WRF/ERA5 via xarray, or a constant/synthetic provider): the
  fire-relevant + convective state (wind, T, RH, surface heat flux, CAPE/CIN, PBL,
  stability) and the physics to derive pyflam inputs (NFDRS equilibrium fuel
  moisture, midflame wind, Monin-Obukhov stability, ambient buoyant heat flux,
  convective plume/spotting enhancement). Drives `fire_atmosphere_march` for
  near-real-time or reanalysis runs.
- `pyflam.meteo_report` — near-real-time fire-weather variation report: samples an
  atmosphere provider across the run window and tracks how temperature, RH, wind,
  dead fuel moisture (per time lag), convective state (CAPE/CIN/PBL/stability) and
  energy fluxes *change* (min/max/mean/range/net-change/trend + the full series).
- `pyflam.operative` — operational analysis of a run perimeter: splits it into
  head / flanks / tail (or finer sub-sectors) and decomposes the spread drive into
  three vectors — **slope**, **fuel** (gradient of intrinsic fireline intensity →
  toward heavier/faster fuel), **wind** — plus their **resultant** and the dominant
  driver per sector. The vectors are the "arrows" a map front-end draws along the
  perimeter (with an optional matplotlib quiver).
- `pyflam.validate` — cell-by-cell comparison of pyflam fields against real
  FlamMap output rasters (bias/RMSE/ratios/correlation/OLS, log-space and
  "within X%" stats, burn/no-burn classification), a parameter scan to recover
  unknown run settings, **perimeter overlap** (`compare_perimeters`: Jaccard /
  Dice / Hausdorff) and **arrival-time agreement** (`compare_arrival_times`) for
  the MTT growth engine vs a FlamMap single-fire (spotting-off) export. See
  `tests/REFERENCE.md` and `tests/validate_flammap_perimeter.py`.
- `pyflam.units` — English/SI unit conversions.

## Testing

```bash
pytest
```

The suite has two layers: physics/property tests (model must obey known
relationships) and a golden-master regression locking current numeric outputs.
See `tests/REFERENCE.md` for how to validate against real FlamMap output — that
diff is the true acceptance criterion before trusting the numbers.

## Roadmap

| Step | Scope | Status |
|------|-------|--------|
| **1** | Rothermel surface fire model + standard 13 fuel models | ✅ |
| **1b** | Scott & Burgan (2005) 40 fuel models + dynamic herbaceous curing | ✅ |
| **2** | Landscape I/O (`.lcp` + GeoTIFF) and vectorized whole-landscape surface behavior | ✅ |
| **2b** | Native mass-consistent wind solver (`pyflam.windsolver`) — fast terrain winds, no binary | ✅ |
| **2c** | Momentum/RANS CFD wind solver (`pyflam.cfd`) — OpenFOAM ABL, stability + diurnal slope flows | ✅ |
| **3** | Directional (vector) spread (`pyflam.mtt.spread_field`) + per-fuel/canopy wind adjustment factor (`pyflam.wind_reduction`) | ✅ |
| **4** | Crown fire: Van Wagner (1977), Rothermel (1991), Scott & Reinhardt (2001) | ✅ |
| **5** | Spread engine: Minimum Travel Time (Finney 2002), elliptical wavelets (Finney 1998) | ✅ |
| **6** | Validation harness vs. real FlamMap rasters (`pyflam.validate`); ROS ~3%, max-spread direction ~1° | 🚧 |

Step 5 implements MTT with ember spotting (`pyflam.spotting`) and a fuel-load
factor; FARSITE-style explicit perimeter looping is an alternative not provided.
Against a real FlamMap run (`tests/REFERENCE.md`): surface ROS matches to ~3% and
max-spread direction to ~1° over 1.6M cells. Burn probability is only partly
recoverable — adding spotting raises it ~36× (to within ~3× of FlamMap), but its
parameters are calibrated and the cell-level pattern is Monte-Carlo-noise-limited;
a crown-fire diff and a spotting-off single-fire perimeter diff remain.

## References

- Rothermel, R.C. 1972. *A mathematical model for predicting fire spread in
  wildland fuels.* USDA FS Research Paper INT-115.
- Albini, F.A. 1976. *Estimating wildfire behavior and effects.* GTR INT-30.
- Anderson, H.E. 1982. *Aids to determining fuel models for estimating fire
  behavior.* GTR INT-122.
- Scott, J.H.; Burgan, R.E. 2005. *Standard fire behavior fuel models: a
  comprehensive set for use with Rothermel's surface fire spread model.*
  RMRS-GTR-153.
- Andrews, P.L. 2018. *The Rothermel surface fire spread model and associated
  developments: A comprehensive explanation.* RMRS-GTR-371.
- Van Wagner, C.E. 1977. *Conditions for the start and spread of crown fire.*
  Canadian Journal of Forest Research 7: 23-34. (Crown fire initiation.)
- Rothermel, R.C. 1991. *Predicting behavior and size of crown fires in the
  Northern Rocky Mountains.* USDA FS Research Paper INT-438. (Active crown ROS.)
- Scott, J.H.; Reinhardt, E.D. 2001. *Assessing crown fire potential by linking
  models of surface and crown fire behavior.* RMRS-RP-29. (Crown fire type, CFB,
  torching/crowning indices.)
- Finney, M.A. 1998. *FARSITE: Fire Area Simulator — model development and
  evaluation.* USDA FS Research Paper RMRS-RP-4. (Elliptical wavelet spread.)
- Finney, M.A. 2002. *Fire growth using minimum travel time methods.* Canadian
  Journal of Forest Research 32: 1420-1424. (MTT spread engine.)
- Anderson, H.E. 1983. *Predicting wind-driven wildland fire size and shape.*
  USDA FS Research Paper INT-305. (Length-to-breadth ratio.)
- Albini, F.A.; Baughman, R.G. 1979. *Estimating windspeeds for predicting
  wildland fire behavior.* USDA FS Research Paper INT-221. (Wind adjustment factor.)
- Andrews, P.L. 2012. *Modeling wind adjustment factor and midflame wind speed
  for Rothermel's surface fire spread model.* USDA FS RMRS-GTR-266.
- Albini, F.A. 1979. *Spot fire distance from burning trees — a predictive
  model.* USDA FS GTR INT-56. (Firebrand spotting.)
- Albini, F.A. 1983. *Potential spotting distance from wind-driven surface
  fires.* USDA FS Research Paper INT-309.
- Sardoy, N. et al. 2008. *Numerical study of ground-level distribution of
  firebrands generated by line fires.* Combustion and Flame 154. (Lognormal
  firebrand landing distribution.)
- Morton, B.R.; Taylor, G.; Turner, J.S. 1956. *Turbulent gravitational
  convection from maintained and instantaneous sources.* Proc. R. Soc. A 234.
  (Buoyant-plume theory for firebrand loft.)
- Clark, T.L.; Coen, J.; Latham, D. 2004. *Description of a coupled
  atmosphere-fire model.* Int. J. Wildland Fire 13. (Fire-atmosphere coupling.)
- Tarifa, C.S. et al. 1965/1967. *On the flight paths and lifetimes of burning
  particles of wood.* Proc. Combustion Institute. (Firebrand d²-law burnout.)
- Manzello, S.L. et al. 2007+. Firebrand-generation experiments (NIST);
  Tohidi, A.; Kaye, N.B. 2017. *Stochastic modeling of firebrand…* Fire Safety J.
  (Measured firebrand sizes and terminal velocities.)
- Forthofer, J.; Butler, B.; Wagenbrenner, N. 2014. *A comparison of three
  approaches for simulating fine-scale surface winds in support of wildland fire
  management.* (WindNinja.) Int. J. Wildland Fire. https://github.com/firelab/windninja
- Sasaki, Y. 1970. *Some basic formalisms in numerical variational analysis.*
  Mon. Weather Rev. (basis of the mass-consistent wind solver.)
- Sherman, C.A. 1978. *A mass-consistent model for wind fields over complex
  terrain.* J. Appl. Meteorology.
- Richards, P.J.; Hoxey, R.P. 1993. *Appropriate boundary conditions for
  computational wind engineering models using the k-ε turbulence model.* (ABL
  boundary conditions for the RANS solver.)
- Wagenbrenner, N.S. et al. 2016. *Downscaling surface wind predictions from
  numerical weather prediction models in complex terrain with WindNinja.*
  (Momentum/CFD wind over terrain.)

## License

MIT (this reimplementation). FlamMap itself is a separate USDA Forest Service
product; nothing here is derived from its binaries.
