# Dry-pyrocloud fire-induced boundary layer (fireABL) — research report

**Branch:** `feat/dry-pyroconvection-fireabl` · **Status:** diagnostic layer, magnitude recalibrated, multi-fire confirmed
**Grounding data:** GRAF (Bombers de la Generalitat de Catalunya) / Wageningen University & Research pyroconvection pipeline and datasets (`zenodo_6433389`, `zenodo_10823337`, `zenodo_12946249`); Eghdami et al. (2023) WRF-Fire LES.

## 1. Research questions

pyflam's pyroconvection stack was, before this work, **moist-only**: LCL, ABL depth, cap γ-θ, shear, moisture aloft — the Castellnou et al. (2022) ladder for whether a fire forms pyroCu/pyroCb by **condensation**. The literature describes a second, mechanistically distinct route (Castellnou et al. 2022; Castellnou Ribau et al. 2024): a fire's **sensible heat alone** grows a fire-induced boundary layer (the "fireABL") that decouples from the surface with **no condensation** — the *dry pyrocloud* mechanism. This report addresses:

1. **Q1 — Port.** Can the GRAF/WUR stage-4 fireABL encroachment be reproduced faithfully in pyflam, cleanly and vectorised?
2. **Q2 — Physics soundness.** Does the ported mechanism agree with independent, higher-fidelity physics (a coupled WRF-Fire LES of the same event)?
3. **Q3 — Discriminating value.** Does the fireABL decoupling signal add anything over pyflam's existing moist ladder?
4. **Q4 — Magnitude.** Are the absolute fireABL heights trustworthy against sonde observations, and if not, why?
5. **Q5 — Integration.** How does this become a usable layer in the pyflam gridded forecast pipeline without over-claiming?

## 2. Solution

### 2.1 The mechanism (ported)
The fireABL is the fire-forced analogue of daytime convective-boundary-layer growth by **encroachment** (Stull 1988): a heated parcel rises until it reaches neutral buoyancy against the ambient sounding. Two functions in `pyflam.atmosphere`:

- `fire_parcel_theta_excess(F, θ_mean)` — the sensible-heat forcing: convective velocity `w* = (3gFH/2ρθ_v)^(1/3)`, temperature excess `θ' = F/(ρ w*)`.
- `fire_induced_abl_grid(z, θ, θ_mean_below, heat_flux, blh)` — the fireABL top: the first level **above the ABL** where ambient `θ(z)` reaches `θ_mean + θ'`, linearly interpolated, gridded over `(nlev, ny, nx)`.

### 2.2 Two physics fixes over the GRAF original
- **Extrapolation → real profile (Q1/Q4).** The GRAF stage-4 code finds the fireABL by intersecting the parcel with a **linear extrapolation of a 2 km free-troposphere fit**, pushed to 13 km — producing 11 km fireABLs from a 2 km fit. `fire_induced_abl_grid` intersects the **actual `θ(z)` stack** instead (as `parcel_mixing_depth_grid` does), so the result is bounded by the sounding.
- **Front flux → mixed-layer flux (Q4).** The θ' *form* is the standard convective temperature scale and is correct; the **flux** was wrong. GRAF forces the encroachment with the fire-**front** flux `(I/2)/front_depth` (front_depth ≈ 10–20 m, a locally intense value), driving the fireABL ~3× too high. `mixed_layer_fire_flux(I, abl)` spreads the fire's convective power over the **ABL depth** — the scale over which the heat actually mixes.

### 2.3 Deliberately *not* included
- **No dry/moist class label.** The fireABL and a `decoupling = fireABL/ABL` ratio are exposed as **diagnostics**, not a class, because a single-predictor gate is too weak (Q3) and the causal "dry → longer duration" claim is single-lineage/unpublished.
- **No in-plume LCL offset baked in.** The literature +1 km cloud-base offset is **not supported by GRAF's own labels** (a fit prefers ~−0.5 km, §4). The offset is a configurable parameter defaulting to 0 (ambient LCL).

## 3. Tests and validation

| # | Test | Result |
|---|---|---|
| Q1 | **Faithful port** vs reference `run_fireABL.py`, 19 SCQ hours | machine precision (≤ 2.4 m for the grid function; 8.7e-11 m for the closed-form intersection); every intermediate to ~1e-14 |
| Q1 | Extrapolation fix | linear method: 6 hours above the 8 km sounding (runaway); profile method: bounded, 0 |
| Q2 | **LES referee** (Eghdami et al. 2023 WRF-Fire, same SCQ event) | peak fireABL 7999 m ≈ LES cloud top ~8 km; first moist hour 18:00 Z ≈ LES pyroCu onset 18:28 Z |
| Q3 | **Stage A** — decoupling vs GRAF prototype labels, 795 h / 87 fires | ratio separates classes (Kruskal–Wallis p = 9e-8, monotonic); adding sonde decoupling lifts CV-AUC 0.556 → **0.725** over LCL/ABL alone |
| Q4 | **Caveat 1** — multi-fire fireABL vs sonde, 87 fires | r = 0.71, biased low, worsening with severity (pyroCb −1464 m) |
| Q4 | **Caveat 2** — fit the in-plume LCL offset | −498 m (95% CI [−825, −392]); literature +1 km **not** supported |
| Q4 | **Stage C (linear)** — scale-height fit, 835 h | H ≈ 56 m halves the pyroCb bias, CV-validated (MAE 513 → 451) — but not transferable to the profile method |
| Q4 | **Profile-method magnitude** (SCQ, full sounding + sonde) | over-predicts +3620 m, r 0.83; **no scale height fixes it** (best H > 6 km, still biased) |
| Q4 | **Mixed-layer flux** (SCQ) | bias +3584 → **+58 m**, r 0.84 → **0.94**, parameter-free |
| Q4 | **Multi-fire ERA5 confirmation** (SCQ, Martorell, PoblaMassaluca, Torroella; 36 h) | ML flux: bias −15 m, MAE 400 m, r 0.91; decisive on the shallow-ABL case (Torroella +27 m vs front flux −444 m) |
| Q5 | **Live hybrid cycle** (ICON-EU + ICON-2I, Tuscany 2026-07-20, 09/12/15 Z) | decoupling panel rendered; sensible diurnal signal (median 4.3 → 3.0 → 3.7 as the ABL deepens); fireABL bounded |
| — | **Automated suite** | 669 passed, 1 skipped; 6 new fireABL/flux physics tests |

## 4. Limits found — GRAF approach vs pyflam

| Aspect | GRAF/WUR stage 4 | pyflam | Limit found |
|---|---|---|---|
| Root-find | `shapely.LineString` intersection as an ad-hoc solver | closed-form / profile interpolation | GRAF's is an engineering hack; also injects ~0.5 m integer-stepping noise |
| Intersection target | linear **extrapolation** of a 2 km fit to 13 km | the **real** `θ(z)` sounding | GRAF runs the fireABL to 11 km (unphysical) on the hours that matter most |
| Forcing flux | fire-**front** flux `(I/2)/front_depth` | **mixed-layer** flux `(I/2)/ABL` | GRAF's front flux drives the fireABL ~3× too high; front depth is a local, not a mixing, scale |
| Magnitude vs sonde | the ported stage-4 over-predicts ~3× | recalibrated: bias −15 m across 4 fires | GRAF's own **refined** `expected_fireABL` matches sondes (−370 m) — i.e. the public **demo** code is the uncalibrated version |
| dry/moist threshold | ambient LCL | ambient LCL (offset configurable, default 0) | Tory et al. (2018) say cloud base runs above ambient LCL, but GRAF's own labels do **not** support the +1 km offset (fit ≈ −0.5 km) — so it is not baked in |
| "dry → longer duration" | asserted (2024 preprint, one group) | **not** encoded | single-lineage, not peer-reviewed — exposed as a diagnostic, not a claim |
| Portability (broader pipeline) | Windows-only `bindpython.pyd` (Nelson/Rothermel), `pandas.append`, hardcoded paths | pure-Python, vectorised, tested | context from the earlier comparison; not re-solved here |

**Net:** the GRAF mechanism **ranks** columns well (r ≈ 0.8–0.9 vs sondes) and is physically motivated, but its published demo-stage-4 **magnitude** is uncalibrated (~3× high) for two separable reasons — extrapolation and front-flux forcing — both fixed here.

## 5. Integration and implementation in the pyflam pipeline

### 5.1 Library (`src/pyflam/atmosphere.py`)
- `fire_parcel_theta_excess`, `fire_induced_abl_grid`, `mixed_layer_fire_flux` — exported from `pyflam`.
- `_FIRE_PLUME_SCALE_M` carries the full calibration note (front-flux bias, the mixed-layer-flux fix, and the multi-fire numbers).

### 5.2 Gridded forecast engine (`pyflam_gui/core/pyroconv.py`)
- `profile_diagnostics` (ICON-2I pressure-level path) **and** `iconeu_diagnostics` (ICON-EU model-level path) now return `fireabl` and `decoupling` alongside the moist ladder's fields, computed on the real `θ(z)` stack with a documented reference flux `_REFERENCE_FIRE_FLUX_W_M2`.
- `regrid_diagnostics` carries both new fields onto the target grid (so the hybrid product keeps them).

### 5.3 Daily product (`tests/pyroconv_daily.py`)
- `export_diagnostics` writes per-hour `diag_fireabl_*.tif` / `diag_decoupling_*.tif` GeoTIFFs.
- `render_decoupling` produces an 8-hour continuous heatmap; `build_pdf` embeds it beside the moist class maps as a **diagnostic** panel with a caption stating the caveats.

### 5.4 How to read it
The decoupling panel is a **diagnostic, not a class**: its spatial and diurnal **pattern** is trustworthy (high-decoupling cells that the moist ladder scores low are the dry-pyrocloud columns); its **absolute metres** should be read cautiously on the gridded product (the reference-fire flux is a convention). For a *real* modeled fire, force the encroachment with `mixed_layer_fire_flux(I, abl)`, not the front flux.

## 6. Open items
- Multi-fire confirmation rests on 4 fires / 36 h; the decisive shallow-ABL evidence is one fire (Torroella, n=5).
- `observed fire_ABL` is GRAF/sonde-derived (not fully independent); a truly independent test (a non-GRAF case with resolved LES, e.g. Vaz et al. 2025 Portugal) is the next external check.
- The gridded diagnostic's reference-fire flux framing (a fixed flux assuming a fire everywhere) is a convention; a per-cell reference intensity ÷ ABL is the principled successor once the intensity-unit conventions are pinned.
