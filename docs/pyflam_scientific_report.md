# pyflam — A Scientific, Technical and Operational Report

*An open-source Python reimplementation of wildland-fire behavior science, from
the Rothermel surface model through coupled fire–atmosphere pyroconvection.*

**Status:** roadmap steps 1–5 implemented; step 6 (validation) ongoing. ~8,600 lines
of source across 25 modules; ~568 automated tests. Pure-Python core (NumPy + SciPy);
optional geospatial, atmospheric, and JIT extras; OpenFOAM and Herbie discovered at
runtime. MIT-licensed; CI on Python 3.11/3.12/3.13.

---

## 0. Executive summary

pyflam began as an open, testable reimplementation of the fire-behavior science
behind **FlamMap** — the USDA Forest Service / Missoula Fire Sciences Laboratory
desktop fire-behavior system — built up from the published, peer-reviewed equations
rather than by decompiling the Windows binary. FlamMap is the operational reference
point: it defines the canonical product set (surface rate of spread, fireline
intensity, flame length; crown-fire potential; minimum-travel-time fire growth;
random-ignition burn probability) and the file conventions (`.lcp` landscapes,
`.fms` moisture) that pyflam reads and reproduces.

From that operational baseline pyflam advances in directions FlamMap does not take:
(1) a **weather-driven** pipeline that derives fuel moisture and wind from live
forecast/reanalysis data and conditions them per cell for terrain and canopy;
(2) a **coupled fire–atmosphere** capability — native terrain-wind solvers, a
buoyant-CFD plume that the fire reshapes, and ember spotting that emerges from the
energy budget; (3) a **literature-current crown-fire** model (Cruz et al. 2005/2004)
that corrects the well-documented under-prediction bias of the FlamMap operational
stack; (4) an **anisotropic-Eikonal** alternative to the lattice-Dijkstra spread
solver that removes most of its grid bias; and (5) **vertical-profile pyroconvection
diagnostics** (LCL, Continuous Haines, inverted-V, a Briggs plume and a pyroCb
firepower threshold) that drive a spatially varying plume feedback. The scientific
value is twofold: pyflam is a *transparent, reproducible* implementation of the
established models, and a *research platform* for the fire–atmosphere coupling that
the operational tools approximate or omit.

---

## 1. Rationale: why reimplement, and FlamMap as the baseline

FlamMap is closed Windows/C++ software, but the science it runs on is open and
documented in USDA Forest Service publications. Reimplementing from the equations
yields clean, cross-platform, unit-tested code; decompilation would yield unreadable
machine output and reproduce platform quirks. Validation is therefore done by
**diffing pyflam's outputs against real FlamMap rasters**, not by reading FlamMap's
code.

What pyflam inherits from the FlamMap lineage as its operational specification:

| Operational product | FlamMap | pyflam |
|---|---|---|
| Surface rate of spread / intensity / flame length | ✅ | `rothermel`, `landscape.basic_fire_behavior` |
| Crown fire potential (surface/passive/active) | ✅ | `crownfire` |
| Minimum Travel Time fire growth | ✅ | `mtt.minimum_travel_time` |
| Random-ignition burn probability | ✅ | `mtt.burn_probability` |
| Landscape `.lcp` / moisture `.fms` I/O | ✅ | `io_lcp`, `landscape` |
| Dead-fuel-moisture conditioning | ✅ | `fuel_conditioning` |

Everything beyond that table is where pyflam departs from FlamMap (Sections 5–9).

---

## 2. Scientific core: the Rothermel surface model

The scientific heart of "Basic Fire Behavior" is the **Rothermel (1972)** surface
fire-spread model with Albini's (1976) refinements: a quasi-steady energy balance in
which the rate of spread is the ratio of the heat flux received by unburned fuel to
the heat required for ignition,

> R = I_R · ξ · (1 + φ_w + φ_s) / (ρ_b · ε · Q_ig),

with reaction intensity I_R, propagating flux ratio ξ, wind and slope factors
φ_w, φ_s, bulk density ρ_b, effective heating number ε, and heat of pre-ignition
Q_ig (Rothermel 1972). pyflam implements the full multi-size-class, multi-category
(dead/live) formulation, exposing it both as a one-shot `spread(...)` call and as a
reusable `SurfaceKernel` — the wind/slope-independent terms computed once per fuel +
moisture and applied to scalar or array (per-cell) wind and slope. Fireline intensity
follows Byram (1959): I = H · w · R; flame length follows L = 0.45·I^0.46.

Fuel inputs use the two operational standard sets: **Anderson's (1982) 13** fuel
models and **Scott & Burgan's (2005) 40**, including the dynamic live-herbaceous
curing transfer. A **fuel-load factor** (scalar, per-fuel, or per-cell raster) scales
the standard loads, which under-represent real loading by ~20–30%; because the
Rothermel response is non-monotonic in load (compact litter can slow past the optimum
packing ratio), this flows through the full kernel rather than scaling the output.

*Validation.* Against a real FlamMap run on the Tuscany landscape (1.6M cells),
surface ROS matches to ~3% (regression slope 0.98, r 0.9998).

---

## 3. Landscape and the directional spread field

A `Landscape` is the in-memory band stack (fuel model, slope, aspect, elevation,
canopy cover/height/base-height/bulk-density). pyflam reads and writes FlamMap/FARSITE
`.lcp` files in pure Python (`io_lcp`, no GDAL needed) and LANDFIRE/arbitrary GeoTIFFs
via rasterio, preserving the coordinate reference system for geolocation.

The transition from a *scalar* `1 + φ_w + φ_s` to a *directional* fire is roadmap
step 3 and the basis of fire growth. pyflam combines the Rothermel wind factor
(blowing downwind) and slope factor (pushing upslope, from the landscape aspect) as
**vectors**, giving each cell a maximum spread rate, a heading azimuth, and an ellipse
**eccentricity** from the Anderson (1983) length-to-breadth ratio (capped, as in
FARSITE). The directional rate of spread off the heading ψ is

> R(ψ) = R_max · (1 − e) / (1 − e·cos ψ),

the elliptical wavelet form (Finney 1998). This `SpreadField` is the object every
downstream solver consumes. *Validation:* the per-cell maximum-spread direction
matches FlamMap's `MAX_SPRE_DIR` to ~1° (mean 0.96°) over 1.6M cells.

---

## 4. Fire growth: Minimum Travel Time, and a novel Eikonal alternative

### 4.1 Minimum Travel Time (the FlamMap engine)

Fire arrival time is the **minimum-time path** from the ignition over a lattice of
travel directions, where the time to cross a segment is its length divided by the
harmonic mean of the elliptical spread rate at its endpoints (Finney 2002). pyflam
assembles the travel-time graph with vectorized NumPy (chunked, direct-CSR for huge
grids) and solves the shortest path with SciPy's C-level multi-source Dijkstra, so it
scales to multi-million-cell landscapes (the 5.35M-cell Tuscany landscape solves from
one ignition in ~2.6 s). `max_time` bounds the search for fixed-duration runs.

### 4.2 The scientific limitation, and the novel fix

MTT is *Dijkstra on a lattice*, and a graph offers only as many travel directions as
its neighbor template — so arrival times carry an **angular (lattice) discretization
bias**: a calm fire becomes a faceted polygon and off-lattice bearings are biased.
This is the "different calculation-point juxtapositions" Finney (2002) noted, and
exactly what the numerical-analysis literature on the **anisotropic Eikonal equation**
was built to remove. Fire arrival time satisfies a static Hamilton–Jacobi equation;
the elliptical-with-wind spread law is precisely a **Randers–Finsler metric**
(an ellipse plus a drift), as formalized for fire by Gahtan et al. (2026) and the
Finsler-geodesic-spray literature.

pyflam therefore offers a **selectable propagation engine**: `method="mtt"` (default)
or `method="fast_marching"`, a semi-Lagrangian anisotropic-Eikonal front solver that
consumes the same `SpreadField`. Benchmarked against the analytic
`distance / R(bearing)` on a uniform field, the Eikonal backend roughly **halves**
MTT's bias on a wind-driven ellipse (mean ~4% vs ~12%) and beats even a denser MTT
template. Three interchangeable backends (Numba-JIT, with a NumPy fallback) give
identical fields; the default **heap** backend is a narrow-band Fast-Marching pass
that **prunes with `max_time` like Dijkstra** followed by a bounding-box Gauss–Seidel
correction — for a bounded single fire it is ~2× faster than MTT *and* more accurate
(MTT's lattice paths are longer, so it under-predicts extent). This is, to our
knowledge, a contribution not present in the operational tools: a Finsler-Eikonal
fire-spread solver offered alongside, and validated against, the classical MTT.

*References:* Finney (2002); Sethian & Vladimirsky (2003, Ordered Upwind Methods);
Mirebeau (2014, Finsler fast-marching); Gahtan, Shpund & Bronstein (2026,
differentiable Randers–Finsler Eikonal solvers).

### 4.3 Burn probability and connected metrics

`burn_probability` reproduces a FlamMap MTT random-ignition run's **full output set**:
burn probability plus the *connected* metrics — conditional flame length, conditional
fireline intensity, the per-class flame-length probabilities (FlamMap's `FLP_METRIC`),
and the per-fire size distribution. Beyond FlamMap it accepts a **weather ensemble**
(per-fire weather variation, which is what makes the flame-length distribution
non-degenerate) and solves fires in **batched multi-source Dijkstra** calls. The
conditional fireline intensity matches FlamMap's `FIRE_LINE_INT` mean to ~2% — the
strongest connected-metric agreement obtained so far.

---

## 5. Weather-driven fuel moisture and wind (beyond FlamMap)

FlamMap conditions dead fuel moisture from a few user-typed values; pyflam derives
the moisture **from weather** and conditions it **per cell**.

### 5.1 Atmospheric forcing

`atmosphere` ingests forecast/reanalysis data — **live GFS** via Herbie (open NOAA
data, no auth), **ERA5** via the Copernicus CDS, or any xarray NetCDF/GRIB column —
behind a provider interface, with a `ConstantAtmosphere` for testing. It carries the
fire-relevant surface and convective state (wind, T, RH, surface heat flux, CAPE,
CIN, PBL height) and derives pyflam's inputs: NFDRS equilibrium dead fuel moisture
(Simard 1968), midflame wind, Monin–Obukhov stability, and the ambient buoyant heat
flux. A time-lag (Nelson-type) `DeadFuelMoistureModel` lets the 1/10/100-h classes
*remember* recent humidity rather than snapping to the latest value.

### 5.2 Per-cell dead-fuel-moisture conditioning

`fuel_conditioning` is the analog of FlamMap's "dead fuel moisture conditioning,"
grounded in a dedicated literature review (Holden & Jolly 2011; Rothermel 1983;
Resco de Dios 2015; Nolan 2016). Sun-exposed fuels (south-facing, open, steep toward
the sun) absorb more shortwave, run warmer than the air, and equilibrate to a *lower*
moisture; shaded fuels (north-facing, under canopy, at night) stay near ambient. The
module computes per-cell solar geometry (`solar_position`, with an equation-of-time /
clock-time correction), a terrain-insolation factor, a canopy-transmission shading
term, and conditions the moisture with either the **NFDRS EMC** or a semi-mechanistic
**VPD** submodel (Resco de Dios 2015 / Nolan 2016). `condition_from_weather` is the
run-setup front door: it derives the initial moisture for a **date/time/location**
from a meteo provider (live GFS/ERA5, sampled per cell on a geolocated landscape) or,
when no data exist, from manually entered T/RH. On the real Tuscany landscape this
produces a ~2× fine-fuel-moisture spread across one weather (e.g. 1.4–14% under a hot
dry GFS column), instead of a single landscape-wide value.

### 5.3 Native terrain-wind solvers

Real terrain bends the wind; a single uniform value is a poor approximation. pyflam
computes a gridded `WindField` two ways, both feeding the same midflame interface and
following the WindNinja science (Forthofer et al.) reimplemented natively:
`windsolver` — a **mass-consistent** solver (fast, no external dependency); and `cfd`
— a **momentum/RANS** solver via OpenFOAM (atmospheric boundary-layer inlet, stability,
diurnal slope flows, per-cell roughness from the fuel grid).

---

## 6. Crown fire: correcting the operational under-prediction

### 6.1 The problem with the FlamMap stack

FlamMap classifies crown fire with **Van Wagner (1977)** initiation + **Rothermel
(1991)** active spread (R_active = 3.34·R₁₀) + **Scott & Reinhardt (2001)**
surface/passive/active classification. Cruz & Alexander (2010) showed that this
operational stack — naming FlamMap, FARSITE, NEXUS, BehavePlus and FFE-FVS explicitly
— carries a **significant crown-fire-spread under-prediction bias**, from three
sources: incompatible surface→crown model linkages, the inherent under-prediction of
the Rothermel ROS models, and an *unsubstantiated* reduction of crown spread by crown
fraction burned. The Van Wagner *initiation* physics is sound; the *spread* side is
outdated.

### 6.2 The pyflam approach

pyflam keeps the FlamMap stack as the default but makes the active-spread model
**selectable** (`crown_spread="rothermel1991" | "cruz2005"`). The **Cruz, Alexander &
Wakimoto (2005)** model — `CROS_active = 11.02·U₁₀^0.90·CBD^0.19·exp(−0.17·M)`
(verified against the primary paper) — is the validated successor; on the Cruz path
an active crown fire spreads at the *full* Cruz rate, dropping the unsubstantiated
crown-fraction reduction. A **CFIS logistic crown-fire-initiation** model (Cruz et al.
2004, coefficients verified at source) provides a *probability* of crowning as an
alternative to the deterministic Van Wagner threshold. On a real GEDI-derived canopy
landscape (Section 10) the Cruz model classifies markedly more cells as active crown
fire than Rothermel — the under-prediction made visible.

Because FlamMap and the whole Rothermel+Van Wagner stack are *known to under-predict*
crown fire, pyflam deliberately **does not validate crown fire against FlamMap**
(that would validate toward a biased reference). Instead its crown fire is *faithful
to the literature-validated Cruz models* (unit-tested against their published
coefficients), so it inherits their validation against observed wildfires (Cruz 2005:
~61% of variance over 57 wildfire observations).

---

## 7. Ember spotting from the energy budget

`spotting` provides two firebrand models that both couple into MTT growth and burn
probability. `SpottingModel` is a fast parameterized loft-and-drift model. The novel
one is `FirebrandPhysics`: a **stochastic, physics-based** model in which spotting
*emerges* from the energy system — buoyant-plume loft from fireline intensity
(Morton–Taylor–Turner), firebrand size → terminal velocity (drag) → loft and
combustion burnout (Tarifa's d²-law), downwind transport, and a landing-ignition
probability that falls with the receiving fuel's moisture. Randomness (Poisson brand
counts ∝ intensity, lognormal sizes, turbulent bearing, Bernoulli ignition) makes the
landing pattern a Monte-Carlo outcome; the constants are *physical*, tied to measured
firebrand data (Tohidi & Kaye and Manzello terminal velocities; Manzello size
distribution), not reach calibrations. The one weakly-constrained length is calibrated
against literature spot-distance anchors.

---

## 8. Coupled fire–atmosphere: pyroconvection and crown–plume–spotting feedback

This is pyflam's most distinctive scientific contribution and is absent from the
operational tools.

### 8.1 The plume the fire makes

`pyroconvection` turns the fire's fireline intensity into a convective ground
heat-flux field, runs a **buoyant RANS** (OpenFOAM `buoyantBoussinesqSimpleFoam`),
and returns the wind *including the fire's own plume* (indrafts/updraft) — the
feedback in which the fire reshapes the wind that drives it (verified end-to-end: a
hot patch shifts the mean near-surface wind 2.43 → 3.14 m/s). `fire_atmosphere_march`
time-marches this coupling, re-solving the plume wind every `dt` minutes of MTT growth,
with the wind solver injectable so the loop is usable and **testable without
OpenFOAM**.

### 8.2 Crown → plume → wind → crown feedback

The regime where the plume matters most is crown fire. pyflam closes a genuinely novel
loop (`docs/crown_plume_coupling.md`): with `crown=True` each march step rebuilds a
**crown-aware spread field** from the current plume-modified wind — crowning cells
spread at the Cruz rate and carry the much-higher crown intensity, which feeds straight
into the plume solver *and* the ember spotting, closing **crowning → stronger plume →
higher wind → faster crown spread**. The positive feedback is bounded by
under-relaxation of the wind and a physical cap (a synthetic stress test confirms no
runaway). On real data the crown intensity ran ~40× the surface value, giving far
greater spotting reach.

### 8.3 Vertical-profile pyroconvection diagnostics

A dedicated review of vertical/drift wind and the boundary-layer/condensation-level
relationship established the key, counter-intuitive result: deep convective (pyroCb)
fire is a **vertical** problem — a deep, dry, well-mixed boundary layer (high LCL)
capped by moisture aloft (the "inverted-V" sounding) with elevated lower-tropospheric
instability — **not** high surface CAPE (pyroCb routinely form with near-zero surface
CAPE; mid-level humidity is the discriminator). pyflam implements diagnostics keyed on
that geometry: `lcl_height_m`, the **Continuous Haines** index (Mills & McCaw 2010),
an **inverted-V** detector, and `pyroconvection_potential` (Castellnou et al. 2022;
Peterson et al. 2017). It adds a **Briggs (1969) bent-over plume** — the validated
wildfire-plume form (Lareau & Clements 2017) — and a **PyroCb Firepower Threshold**
(Tory & Kepert 2021): the minimum firepower for the plume to reach condensation
against the capping stability.

These are **wired into the coupled march**: with `pyroconvection=True` the plume
intensity is scaled by a profile-aware `convective_plume_factor`, so a dry/unstable
inverted-V atmosphere drives a stronger plume and a stable one damps it — and in
`spatial=True` mode the factor is **per cell**, so under one moist-aloft column the
locally dry (high-LCL) cells get the pyroconvective boost while moist cells do not.
This re-weights the convective coupling toward the predictors the pyroCb literature
favors over surface CAPE.

---

## 9. Operational products and near-real-time use

pyflam is engineered for analyst-facing operation, not only research:

- **`meteo_report`** — a near-real-time fire-weather variation report sampling the
  atmosphere across the run window, tracking how T, RH, wind, dead fuel moisture (per
  time lag), convective state (CAPE/CIN/PBL/stability) and energy fluxes *change*.
- **`operative`** — operational analysis of a run perimeter: it splits the perimeter
  into head / flanks / tail (or finer sub-sectors) and decomposes the spread drive
  into **slope**, **fuel** (gradient of intrinsic fireline intensity) and **wind**
  vectors plus their resultant and dominant driver — the "arrows" a map front-end
  draws — with **GeoJSON export** (sector forces, contour-traced perimeter polygon,
  reprojected to WGS84).
- **`nrt.run_realtime`** — a one-call near-real-time product tying weather →
  time-marched spread → perimeter → both reports into a single `RunProduct`.

An end-to-end operational scenario is reproducible via `tests/final_run_tuscany.py`:
300 random ignitions on the real 5.35M-cell Tuscany landscape, 28 June 2026, a 36-hour
fire driven by **live GFS** weather with per-cell conditioned moisture, writing each
pipeline phase to its own structured output folder (landscape → weather → moisture →
surface behavior → wind → spread field → ignitions → burn probability + metrics).

---

## 10. Real-data canopy fuels (GEDI) and the crown pipeline

No global canopy-base-height / bulk-density product exists — GEDI and the
GEDI-calibrated canopy-height maps give *height* only; CBH/CBD are either US-LANDFIRE
or must be derived. pyflam includes a reproducible path (`tests/fetch_canopy_tuscany.sh`,
`tests/build_canopy_landscape_tuscany.py`) that downloads **GEDI-calibrated canopy
height** (Meta/WRI High-Resolution Canopy Height, Tolan et al. 2024, open AWS bucket),
warps it onto the Tuscany grid, and **derives** CBH/CBD from height + canopy cover
with transparent, documented fire-science heuristics (crown ratio rising with cover →
lower crown base; CBD = cover-scaled load ÷ canopy depth). The result is a full canopy
`.lcp` on which the entire crown pipeline runs — read → Cruz vs Rothermel
classification → crown-aware spread field → plume-coupled crown march. (CBH/CBD here
are *derived estimates*, clearly flagged as such, not field measurements.)

---

## 11. Validation strategy

The acceptance test is a cell-by-cell diff against real FlamMap output, not internal
self-consistency. `validate` provides the generic machinery — robust field comparison
(bias/RMSE/ratios/correlation/OLS, log-space and "within-X%" stats, burn/no-burn
classification), perimeter overlap (Jaccard/Dice/Hausdorff), arrival-time agreement,
and categorical agreement (confusion matrix + per-class recall, for crown-fire type) —
with data-wiring scripts per landscape (`tests/validate_flammap_*.py`).

Validated to date (Tuscany, 1.6M cells): **surface ROS ~3%**, **max-spread direction
~1°**, **conditional fireline intensity ~2%** vs FlamMap. Partially recoverable: burn
probability (Monte-Carlo and spotting-parameter limited). Open: a crown-fire diff
against *observed* (not FlamMap) ROS, and a spotting-off single-fire perimeter /
time-of-arrival diff (the bundled dataset ships no arrival-time raster). The two test
layers — physics/property tests asserting known relationships, and golden-master
regressions locking current numerics — total ~568 automated tests run in CI across
Python 3.11/3.12/3.13.

---

## 12. Software engineering

Pure-Python core on **NumPy + SciPy**; optional extras gate heavier capabilities —
`geo` (rasterio/pyproj/scikit-image), `atmos` (xarray/cfgrib/netcdf4/cdsapi), `accel`
(Numba JIT for the Eikonal solver). External engines (OpenFOAM, Herbie/ERA5) are
discovered at runtime and their tests self-skip when absent, so the core installs and
runs anywhere. The vectorized kernel design (a `SurfaceKernel` computed once and
applied to per-cell array inputs) is what makes whole-landscape moisture conditioning,
weather forcing and crown classification tractable at multi-million-cell scale. CI is
SHA-pinned to Node-24 actions; coverage is reported to Codecov.

---

## 13. Summary of novel contributions relative to FlamMap

1. **Weather-driven, per-cell fuel moisture** — live GFS/ERA5 → terrain-insolation +
   canopy-shading conditioning (EMC or VPD), vs FlamMap's typed scalars.
2. **Selectable anisotropic-Eikonal spread engine** — a Finsler-Eikonal solver
   alongside MTT, removing most of the lattice bias and, in heap form, beating
   Dijkstra on bounded fires in both speed and accuracy.
3. **Burn probability with a weather ensemble and the full connected-metric set.**
4. **Literature-current crown fire** — Cruz 2005 active spread + Cruz 2004 logistic
   initiation, correcting the documented under-prediction of the operational stack.
5. **Coupled fire–atmosphere** — native terrain winds, a buoyant-CFD plume the fire
   reshapes, physics-based spotting, and a **crown → plume → wind → crown** feedback.
6. **Vertical-profile pyroconvection** — LCL/C-Haines/inverted-V diagnostics, a Briggs
   plume and a pyroCb firepower threshold, driving a per-cell plume feedback keyed on
   the ABL/LCL geometry rather than surface CAPE.
7. **Operational NRT layer** — weather-variation and perimeter-driver reports with
   GeoJSON export, plus a one-call real-time product.

The scientific value is that each of these is grounded in the peer-reviewed
literature and implemented transparently and reproducibly, with the established
FlamMap models retained as defaults and the novel methods offered as validated,
selectable alternatives.

---

## References

- Albini, F.A. (1976). *Estimating wildfire behavior and effects.* USDA FS GTR INT-30.
- Anderson, H.E. (1982). *Aids to determining fuel models for estimating fire behavior.* USDA FS GTR INT-122.
- Anderson, H.E. (1983). *Predicting wind-driven wildland fire size and shape.* USDA FS RP INT-305.
- Briggs, G.A. (1969). *Plume Rise.* USAEC TID-25075.
- Byram, G.M. (1959). *Combustion of forest fuels.* In *Forest Fire: Control and Use.*
- Castellnou, M.; Stoof, C.R.; Vilà-Guerau de Arellano, J.; et al. (2022). *Pyroconvection classification based on atmospheric vertical profiling.* J. Geophys. Res. Atmos. 127, e2022JD036920.
- Cruz, M.G.; Alexander, M.E.; Wakimoto, R.H. (2004). *Modeling the likelihood of crown fire occurrence in conifer forest stands.* Forest Science 50(5), 640–658.
- Cruz, M.G.; Alexander, M.E.; Wakimoto, R.H. (2005). *Development and testing of models for predicting crown fire rate of spread.* Can. J. For. Res. 35(7), 1626–1639.
- Cruz, M.G.; Alexander, M.E. (2010). *Assessing crown fire potential in coniferous forests of western North America: a critique of current approaches.* Int. J. Wildland Fire 19, 377–398.
- Finney, M.A. (1998). *FARSITE: Fire Area Simulator — model development and evaluation.* USDA FS RP RMRS-RP-4.
- Finney, M.A. (2002). *Fire growth using minimum travel time methods.* Can. J. For. Res. 32(8), 1420–1424.
- Forthofer, J.M.; et al. (WindNinja). *Mass-consistent and momentum diagnostic wind models for wildland fire.*
- Gahtan, B.; Shpund, J.; Bronstein, A.M. (2026). *Wildfire Simulation with Differentiable Randers–Finsler Eikonal Solvers.* arXiv:2603.00035.
- Holden, Z.A.; Jolly, W.M. (2011). *Modeling topographic influences on fuel moisture and fire danger in complex terrain.* Forest Ecology and Management 262, 2033–2041.
- Lareau, N.P.; Clements, C.B. (2017). *The Mean and Turbulent Properties of a Wildfire Convective Plume.* J. Appl. Meteorol. Climatol. 56(8).
- Manzello, S.L.; et al. *Firebrand (ember) size and generation measurements.*
- Mills, G.A.; McCaw, W.L. (2010). *Atmospheric stability environments and fire weather in Australia — the Continuous Haines index.* CAWCR Tech. Rep. 20.
- Mirebeau, J.-M. (2014). *Efficient fast marching with Finsler metrics.* Numerische Mathematik.
- Morton, B.R.; Taylor, G.; Turner, J.S. (1956). *Turbulent gravitational convection from maintained and instantaneous sources.* Proc. R. Soc. Lond. A.
- Morvan, D.; Frangieh, N. (2018). *Wildland fires behaviour: wind effect versus Byram's convective number.* Int. J. Wildland Fire 27(10).
- Nelson, R.M. (2000). *Prediction of diurnal change in 10-h fuel stick moisture content.* Can. J. For. Res. 30, 1071–1087.
- Nolan, R.H.; Resco de Dios, V.; Boer, M.M.; et al. (2016). *Predicting dead fine fuel moisture at regional scales using vapour pressure deficit.* Remote Sensing of Environment 174, 100–108.
- Peterson, D.A.; et al. (2017). *Pyrocumulonimbus climatology — mid-troposphere humidity as a pyroCb discriminator.*
- Resco de Dios, V.; et al. (2015). *A semi-mechanistic model for predicting the moisture content of fine litter.* Agricultural and Forest Meteorology 203, 64–73.
- Rothermel, R.C. (1972). *A mathematical model for predicting fire spread in wildland fuels.* USDA FS RP INT-115.
- Rothermel, R.C. (1983). *How to predict the spread and intensity of forest and range fires.* USDA FS GTR INT-143.
- Rothermel, R.C. (1991). *Predicting behavior and size of crown fires in the Northern Rocky Mountains.* USDA FS RP INT-438.
- Scott, J.H.; Reinhardt, E.D. (2001). *Assessing crown fire potential by linking models of surface and crown fire behavior.* USDA FS RMRS-RP-29.
- Scott, J.H.; Burgan, R.E. (2005). *Standard fire behavior fuel models.* USDA FS GTR RMRS-GTR-153.
- Sethian, J.A.; Vladimirsky, A. (2003). *Ordered Upwind Methods for Static Hamilton–Jacobi Equations.* SIAM J. Numer. Anal. 41(1), 325–363.
- Simard, A.J. (1968). *The moisture content of forest fuels.* Canadian Dept. of Forestry.
- Tarifa, C.S.; et al. (1965). *On the flight paths and lifetimes of burning particles of wood.* Proc. Combustion Institute.
- Tohidi, A.; Kaye, N.B. *Aerodynamic characterization of firebrands / terminal velocity.*
- Tolan, J.; et al. (2024). *Very high resolution canopy height maps from RGB imagery (Meta/WRI HRCH).* Remote Sensing of Environment.
- Tory, K.J.; Kepert, J.D. (2021). *Pyrocumulonimbus Firepower Threshold: Assessing the Atmospheric Potential for pyroCb.* Weather and Forecasting 36(2).
- Van Wagner, C.E. (1977). *Conditions for the start and spread of crown fire.* Can. J. For. Res. 7(1), 23–34.

---

## AI-assistance disclosure

This report, and substantial parts of the pyflam implementation it describes, were
produced with the assistance of an AI coding/research tool (Claude). The underlying
fire-science models are reimplementations of the cited peer-reviewed publications; the
literature reviews that motivated the novel components were conducted with AI-assisted
search, and every cited source was located in live searches and (where central)
verified at its primary source. Coefficients of the empirical models were checked
against the original papers before implementation.
