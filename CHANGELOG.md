# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Mari's pyroconvection reference scripts** in `pyroconv_script_MARI/`: ERA5 v1.1, ICON-2I MeteoHub v1.2 and WRF 1.5 km v1.7.4, with the WRF technical report. These are the originals that pyflam's classifier was ported from. They are credited to Mari as developer in the folder README, `NOTICE`, `README.md`, `CITATION.cff` (references) and `.zenodo.json` (contributor). They are not under the AGPL, and all rights remain with Mari.
- **PyroCb Firepower Threshold, gridded** (`atmosphere.pyrocb_firepower_threshold_grid`) —
  Tory & Kepert (2021) eq 31 with their full section 4b–4c procedure: height-weighted mixed
  layer iterated until the ML-LCL sits inside it, the Tory et al. (2018) saturation-point curve
  at 15 K per g kg⁻¹, a new pseudoadiabat solver, the β-walk to the free-convection height, the
  0.5 K buoyancy buffer, and the vector-mean wind to z_fc. Reproduces five of six published
  manual analyses exactly (Black Saturday 1244 vs 1240 GW, Sedgerly Rd 12 vs 12 GW). Exported
  as a raster with `z_fc`, `delta_theta_fc`, `u_ml`; **never gated on** — the authors state the
  absolute values are unverified and recommend the field for relative threat. The earlier
  scalar `pyrocb_firepower_threshold` targeted the LCL with a bulk stability, which are not the
  paper's quantities and understate PFT by orders of magnitude.
- **ABL by maximum relative humidity** (`atmosphere.max_rh_abl_grid`) — the criterion
  Castellnou Ribau et al. (2025) sec. 2.6 actually use, a moisture rather than a dynamic one.
  Across 59 ambient campaign sondes it agrees with the bulk-Richardson depth to 9 % on
  well-mixed columns but runs 3.8× deeper on collapsed ones, and on the Tuscany grid it moves
  the 18Z median LCL/ABL from 4.37 to 0.55 — convection plume to resilient pyroCu. Provided for
  comparison; the classifier still uses the Rib depth pending a hand check against plotted
  profiles, since the source applies the criterion visually.
- `ICON_EU_MODEL_LEVELS` extended 74→51 to **74→36** (~4 km → ~8.9 km). The PFT needs the
  profile to the −20 °C electrification level; on the short stack every column returned nan,
  and an intermediate 74→39 still left 29 % untestable. ~148 → ~235 MB per step.

- **Entrainment jump and FireCAPE** (`atmosphere.entrainment_jump_grid`,
  `atmosphere.fire_cape_grid`) — the two variables Castellnou et al. (2022) compute and
  pyflam did not. Δθ across the entrainment zone is what sec. 2.1.2 names as the control on
  a parcel's *ability to penetrate above the ABL*; FireCAPE is Potter (2005), the CAPE of a
  fire-heated parcel (not surface CAPE, which this product still omits). Both are
  **diagnostic only** — exported as rasters, deliberately not wired into the ladder. The
  obvious penetration gate θ′/Δθ ≥ 1 is a no-op over Tuscany (p10 = 2.7 even after the θ′
  correction below), and the blocker is not the forcing: forced with each campaign fire's own
  *measured* excess, the encroachment reproduces observed fireABL at ~151 % relative error
  (n = 4, mean |error| 1156 m on observations of 371–2850 m). The gate cannot rest on a chain
  with no demonstrated skill, so the next question is the encroachment model, not its input.
- **Residual-layer depth** (`atmosphere.residual_layer_grid`) — references the parcel to the
  coldest level in the lower column rather than the surface, so the mixed-layer gradient stays
  well-posed after the convective layer decays. By day it reduces to the parcel mixing depth.
  Evening classifiable land over Tuscany rises 48 → 81 %. The depth it returns still looks too
  shallow to be a true residual layer; see `docs/graf_vs_pyflam_2026-07-26.md`.
- **Dry-pyrocloud decoupling map** in the 3-day Tuscany product, alongside the fuel-gated and
  potential class maps: a continuous fireABL/ABL panel, 3 days × 8 hours from 00Z of the run
  day, with a daytime summary table carrying the classifiable-area share.
- **Per-province decoupling metric** (`decoup_max`) in the provincial report and CSV.

### Changed

- **One reference gauge for dry thermodynamics.** `atmosphere.py` now defines p₀ and the
  Poisson exponent once (`_P0_HPA = 1000`, `_P0_PA`, `_KAPPA_DRY = _RD / _CP_DRY = 0.28591`).
  Previously θ was computed forward with a literal `0.286` (`theta_kelvin`,
  `brunt_vaisala_squared`) but inverted with `_KAPPA_DRY` in the PFT's mixed-layer and
  saturation-point steps, so the module carried two gauges. Dry-air c_p is merged at 1004
  J kg⁻¹ K⁻¹ (`_CP_AIR` is now an alias, previously 1005), and the two `287.0` literals are
  now `_RD`. Output drift: θ −0.005 / −0.011 / −0.023 K at 850 / 700 / 500 hPa; Briggs
  buoyancy flux +0.1 % and Obukhov length −0.1 % (both ∝ c_p^∓1); PFT on the synthetic test columns +0.26 % /
  −0.05 %, below the β-grid noise. A new test pins the gauge and fails if a re-typed
  literal comes back. `spotting.py` keeps its own `_CP_AIR = 1005` (out of scope).
  Recommendation 1 of `docs/deep_research_gauge_thermodynamics_2026-11-19.md`.
- **`pyrocb_firepower_threshold_grid(moisture_reference=)`**, an opt-in research variant.
  `"theta"` (the default) is Tory & Kepert's eq 30. `"surface_t"` applies the 15 K per g kg⁻¹
  fire ratio per kelvin of the mixed layer's actual surface temperature, which removes the
  reference-pressure dependence over high ground. On the GRAF capability-margin set it moves
  the critical growth rate by −0.04 to +0.5 % on 24 of 27 columns and by +5.5 to +15.5 % on
  three (Rojals, Guissona, Orista). Guissona's margin stays positive (+0.39 → +0.34), and no
  detection count, Fisher p or rank correlation changes (`docs/gauge_audit/margin_*.txt`).
  The default is unchanged. The cost ladder still uses the paper form.
- **Reference-pressure covariance tests** (Recommendation 2 of the same report). The tests
  move p₀ and assert that θ rescales as a pure dilation, and that `critical_rh_for_cloud_persistence`,
  `fire_cape_grid` and `bulk_richardson_abl_grid` do not move. A further test pins PFT to
  Tory & Kepert's (2021) own gauge (P₀ = 10⁵ Pa, Table 2; eq 30 applies the 15 K per g kg⁻¹ ratio
  to βθ_ML). That gauge dependence is part of the published method and is reproduced, not
  "fixed".
- **Removed a threshold fitted to a third party's forecast, and annotated the provenance of
  every remaining one.** `PyroconvThresholds.rh_top_moist` was lowered from 80 % to 60 %
  because that "brings the deep-pyroCb fraction into agreement" with the Catalan Bombers
  ICON-EU product on a single day (2026-07-15). That is a fit to someone else's output, on the
  gate governing classes 3 and 4, which means any subsequent agreement with that product was
  partly manufactured. Restored to 80 % — not because 80 % is better supported (it is equally
  uncited) but because a fitted value cannot serve as independent corroboration.

  Consequence, visible rather than tuned away: at 27/07 15Z over Tuscany, class 3 falls
  30.1 % → 14.4 % and class 4 falls **5.0 % → 0.5 %**; total ≥ pyroCu is unchanged at 85.7 %,
  since demoted cells land in class 2. If that is wrong, the fix is a moisture criterion with a
  physical basis, not a number chosen to match a map.

  Every field in `PyroconvThresholds` now carries a provenance grade: `[PRIMARY]` (traceable to
  a peer-reviewed source, with citation), `[DERIVED]` (interpolated from a small number of
  published cases), or `[UNSOURCED]` (inherited from the MARI reference port with no citation
  found). The audit found that only `ml_stable` and `lcl_ratio_resilient_max` are primary; the
  two cap brackets rest on two cases each; and `ml_overshoot_max`, `lcl_ratio_overshoot_max`,
  `lcl_ratio_deep_max` and `shear_distance_deep` have no source in Castellnou et al. (2022) at
  all. They are retained because removing them would change the ladder's structure, but they
  are marked as not evidence.

- **`max_rh_abl_grid` added and available as `iconeu_diagnostics(abl_method="maxrh_floored")`,
  but the bulk-Richardson depth remains the default.** The maximum-RH height is the criterion
  the source method uses (Castellnou Ribau et al. 2025 sec. 2.6, "supplemented with numerical
  calculations using the bulk Richardson number"), here with the supplement applied as a floor
  — a moisture maximum below the dynamically diagnosed mixing top is not a capping inversion.

  It was briefly made the default and then reverted, because measured on the ladder this
  product actually runs it is worse. Against 26 ambient campaign sondes at GRAF-labelled fires,
  on the **adaptive** ladder (which resolves to the 5-diagnostic `shear` path in 51 of 52
  cases): `rib` gives 11/26 exact and mean bias **+0.12**, `maxrh_floored` 7/26 and **+0.58**.
  The ordering reverses on the 3-diagnostic `castellnou` ladder (+1.38 against +1.81), which is
  what an earlier analysis measured — and got wrong by validating a ladder the product does not
  use. A deeper ABL shifts both the cap layer and the RH-top window upward, and on the adaptive
  ladder those gates carry the discrimination, so a better LCL/ABL ratio bought at their expense
  is a net loss.

  Two things worth carrying forward from that episode: on the adaptive ladder the classifier is
  **nearly unbiased (+0.12) and exact 42 % of the time** against observed fire classes — not the
  "1.4 classes hot" reported earlier — and `abl_rib` is now exported alongside `abl` so the two
  depths are comparable in any product.
- **The dry-pyrocloud reference fire is prescribed as a temperature excess, not a heat flux**
  (`_REFERENCE_THETA_EXCESS_K = 10.0`, replacing `_REFERENCE_FIRE_FLUX_W_M2 = 200`), and
  `fire_induced_abl_grid` accepts `theta_excess=` directly. Even dimensionally corrected, a
  flux-derived theta* is an order of magnitude below the 0.1–13.1 K anomalies the GRAF
  campaign measured inside real plumes, and no cell-averaged flux closes that gap — it spreads
  the fire's heat over ground that is mostly not burning. 10 K is the intense end of those
  measurements. This halves the decoupling diagnostic (27/07 15Z median 2.9 -> 2.08).

- **`ABL_MIN_M` 600 m → 150 m.** The 600 m gate ("mixing too shallow to classify") had no
  basis in the source method — Castellnou et al. (2022) sec. 2.4.2 specifies the Rib procedure
  and imposes no minimum depth — and is contradicted by the campaign's own sondes: five of the
  eight ambient profiles released beside real pyroconvective wildfires sit at 202–391 m and
  were being discarded unclassified, though their fires grew observed fire-induced boundary
  layers of 451–2850 m. Now a plausibility floor matching the Rib solver's own clip, not a
  discriminating gate.
- **Mixed-layer level count** is now the land-only median over the 12–15Z mature window, with
  a companion `ml_fit_support` (share of land columns clearing `ML_FIT_MIN_PTS`). It was read
  from hour 0, which sampled the collapsed nocturnal layer and reported ~1 level in the very
  sentence meant to establish the hybrid path's resolution; the honest figure is ~13.
- The daily report's ladder paragraph and shear threshold row are derived from the ladder
  actually used, instead of asserting the pressure-level answer. On the hybrid path the full
  5-diagnostic ladder does run, and the report said otherwise.

### Fixed

- **The Fire Weather Preview treated every cell as sea level.** `pyroconv.classify` computed
  the surface potential temperature as `theta_kelvin(T2m, 1000.0)`, so over terrain θ_sfc was
  too cold by the Exner factor: 4.4 / 9.1 / 14.2 K on 950 / 900 / 850 hPa surfaces. That adds
  3–10× the `_ML_STABLE` gate to the mixed-layer gradient, so hills read as stable and the
  preview was biased toward the quiet classes. `ICON2I_FIRE_FIELDS` now downloads `PS`,
  `read_icon2i` returns it, and `classify(..., ps_hpa=)` uses θ(T2m, ps), places the surface
  at the standard height of ps, and measures the ABL above it. Cells less than 250 m below
  850 hPa return `PYROCONV_NODATA` instead of a class. Without `ps_hpa` the old behaviour is
  kept. The terrain logic lives in one helper, `pyroconv.surface_layer_terms`. It skips
  pressure levels at or below the ground when searching for the ABL, and returns the ABL
  height above the surface. The four stand-alone map scripts (`tests/pyroconv_map_tuscany*.py`,
  `tests/pyroconv_icon2i_*.py`) now read surface pressure (Open-Meteo `surface_pressure`, GFS
  `PRES:surface`, ICON-2I `PS`/`sp`) and call the helper. They render nodata cells in grey. The scripts' data paths, run date, output dir and `.lcp` can be overridden with `ICON2I_DIR`/`ICON2I_GRIB`, `RUN_DATE`, `OUTDIR` and `LCP`.
  The pressure-aware path also measures the ABL **above the ground**. The legacy path counted it from sea level, so every cell got the 110 m height of the 1000 hPa level for free, and shallow stable layers cleared the 150 m `ABL_MIN_M` floor. On the real ICON-2I run of 2026-11-19 12Z (full Tuscany, 13,920 cells) the legacy path called 13,919 cells "convection plume". With surface pressure: 7,984 surface plume (mostly fog-capped lowlands with 50–100 m noon inversions), 4,362 convection plume, 161 overshooting, 496 resilient pyroCu, 588 deep pyroCu/pyroCb on the hills, and 329 high-ground nodata. A deep-class ridge cell checked by hand (900 hPa, 15 °C, θ_sfc 297.3 K vs θ₈₅₀ 296.5 K, LCL ≈ ABL ≈ 880 m) is genuinely well mixed. The legacy path had it stable because it put the 2 m state at 1000 hPa.
- **`critical_rh_for_cloud_persistence` equated a θ excess with a temperature cooling.** The
  latent cooling (L_v/c_p)Δq is a ΔT, but the plume excess is a Δθ. The excess is now
  converted at the evaluation level, ΔT = Δθ (p/p₀)^κ. RH_crit rises by 0.5–4.4 points, e.g.
  54.7 → 59.1 % for a 10 K excess at 700 hPa.
- **`fire_parcel_theta_excess` was dimensionally wrong.** It computed `F / (rho * w0)`, which
  carries units of J/kg rather than kelvin — the specific heat capacity was absent from both
  the velocity and the temperature scale — and the result was added straight to a potential
  temperature in `fire_induced_abl_grid`. For a 200 W/m² flux it returned 20.1 K against a
  correct 0.23 K, an overstatement of ~88×. Now the standard Deardorff scaling, with a
  docstring warning that theta* is a mixed-layer *turbulence* scale of a few tenths of a
  kelvin and must not be used as a plume's temperature excess.

- **Unclassifiable columns are no longer drawn as class 0.** `classify_profile` allocated with
  `np.zeros` and left rejected cells at 0, which the legend renders as "Surface plume", so the
  map could not distinguish *no significant convection* from *could not classify*. At 18Z on
  the 2026-07-26 run this painted 99.2 % of land as quiet when 0.0 % of it was a genuine
  surface-plume diagnosis. Rejected cells now carry `PYROCONV_NODATA` (−1) and render grey; a
  fuel-gated cell below 10 MW/m stays class 0, which *is* a diagnosis.
- **`bulk_richardson_abl_height` used Ri_c = 0.25**, matching neither Castellnou et al. (2022)
  (0.33, after Zhang et al. 2014) nor pyflam's own gridded path. The 200 vs 400 m search start
  is retained as a deliberate split — 400 m is a fire-sonde correction for plume indraft, 200 m
  (Zhang) is right for ambient forecast columns — and both constants now carry the reasoning.
- **The 3-day PDF built without its figures.** `build_report` passed pandoc an absolute `.md`
  path while the image links are basenames, so they resolved against the repo root and the
  build silently produced a text-only PDF.

- **Rothermel effective wind-speed limit** (opt-in). `SurfaceKernel` gains
  `effective_wind_limit()` (0.9·I_R, Rothermel 1972 eq. 87) and
  `limit_combined_factor()`, and `rate_of_spread`/`behavior`/`basic_fire_behavior`/
  `spread_field` take `effective_wind_limit=` (default **False**). With it on, the
  combined wind+slope factor is capped at the effective-wind limit as BehavePlus/
  FlamMap do — taming the Rothermel wind-factor singularity that drives spread to
  unphysical values on steep-slope / high-wind cells in high-SAV fuels (surfaced by
  the Montale 2017 case). Left off by default because the classic 0.9·I_R limit
  over-restricts low-intensity fuels at modest wind (Andrews, Cruz & Rothermel 2013),
  which would change verified spread and the golden-master ROS values.

## [0.2.0] - 2026-07-21

Pre-release for field testing against real recorded wildfires. Integrates the
fire-weather / FWI pipeline and the new dry-pyrocloud fire-induced boundary layer.

> **Correction notice (2026-07-27).** The fireABL physics described below was
> dimensionally wrong at release. `fire_parcel_theta_excess` computed `F / (rho * w0)`,
> which carries units of J/kg rather than kelvin, and the result was added directly to a
> potential temperature — overstating the plume excess by ~88×. The "recalibrated forcing"
> credited below to `mixed_layer_fire_flux` was, in effect, a ~1/88 global coefficient
> standing in for the missing specific heat capacity; that is why the validation quoted
> "one global coefficient for the table's unit ambiguity" and why `r = 0.91` held while the
> absolute scale did not. The `≤2.5 m` agreement with the reference pipeline reproduced the
> reference's own error rather than validating the physics.
>
> Corrected in [Unreleased]. With `cp` restored the chain needs **no fitted coefficient**:
> a 10 MW/m fire over a 1500 m ABL now yields θ′ = 1.31 K against 131 K before, inside the
> 0.1–13.1 K band measured in real plumes by the GRAF campaign, and the resulting fireABL
> tracks observed plume tops to a mean error of ~675 m across the campaign sondes. Anything
> derived from the 0.2.0 fireABL magnitudes — including published decoupling rasters — should
> be regenerated.

### Added

- **Dry-pyrocloud fire-induced boundary layer.** `atmosphere.fire_induced_abl_grid`
  and `atmosphere.fire_parcel_theta_excess` add the sensible-heat-only,
  condensation-free boundary-layer decoupling mechanism (the "fireABL") from the
  GRAF/WUR pipeline (Castellnou et al. 2022; Castellnou Ribau et al. 2024) — the
  dry counterpart to the existing moist LCL/cap/shear classifier. The fireABL top
  is found by intersecting the fire-heated parcel with the actual `theta(z)`
  profile stack (bounded by the sounding), fixing the linear-extrapolation
  runaway in the original. Validated against the reference pipeline on the Santa
  Coloma de Queralt case (≤2.5 m) and against the Eghdami et al. (2023) WRF-Fire
  LES anchors. The in-plume LCL offset for the dry/moist split is left as a
  configurable parameter (default 0 = ambient LCL); the literature +1 km value is
  not baked in, as it is not supported by the GRAF prototype labels.
  **Magnitude — recalibrated forcing.** The reference stage-4 forcing (fire-front
  flux) drives the fireABL too high, and no scale height fixes it. `mixed_layer_fire_flux`
  spreads the fire's convective power over the ABL depth instead of the flaming front.
  Confirmed across the 4 fires with ERA5 soundings + sonde fireABLs (36 h): bias −15 m,
  MAE 400 m, r 0.91, and it correctly captures the shallow-ABL strong-decoupling case
  the front flux misses by ~4×. Recommended forcing for a real fire is
  `mixed_layer_fire_flux(I, abl)`. The gridded `decoupling` diagnostic still uses a
  fixed reference flux (a convention, since it assumes a fire everywhere) — read its
  spatial/diurnal pattern qualitatively.

## [0.1.3] - 2026-07-02

Submission-ready release — no library/behavior changes.

### Changed

- **Paper title sharpened** to lead with the distinctive contribution (coupled
  fire–atmosphere pyroconvection and fire-danger assessment), and a suggested
  JOSS subject/editor-track note added to `paper/SUBMISSION.md`.
  ([#17](https://github.com/crifod/pyflam/pull/17))

### Fixed

- **Citation now uses the Zenodo *concept* DOI** `10.5281/zenodo.21129258`
  (all-versions) instead of the v0.1.0 version DOI, in `CITATION.cff` and the
  README badge/citation. ([#16](https://github.com/crifod/pyflam/pull/16))

## [0.1.2] - 2026-07-02

Publication release — no library/behavior changes.

### Changed

- **JOSS paper restructured to the current JOSS mandatory sections** — adds
  *State of the field* (with an explicit build-vs-contribute justification),
  *Software design*, *Research impact statement*, and an *AI usage disclosure*
  section; trims *Statement of need* to the need itself.
  ([#13](https://github.com/crifod/pyflam/pull/13))

## [0.1.1] - 2026-07-02

Packaging and publication release — no library/behavior changes.

### Added

- **JOSS paper** (`paper/paper.md`, `paper/paper.bib`) describing pyflam, with a
  SHA-pinned GitHub Action (`draft-pdf.yml`) that typesets it to a PDF artifact.
  ([#8](https://github.com/crifod/pyflam/pull/8))
- **Zenodo archiving + citation metadata** — `.zenodo.json` and `CITATION.cff`,
  and the concept DOI [10.5281/zenodo.21129259](https://doi.org/10.5281/zenodo.21129259)
  wired into the README and citation file.
  ([#6](https://github.com/crifod/pyflam/pull/6), [#7](https://github.com/crifod/pyflam/pull/7))
- **JOSS submission aids** — cover note and pre-submission checklist
  (`paper/SUBMISSION.md`). ([#9](https://github.com/crifod/pyflam/pull/9))

## [0.1.0] - 2026-07-02

First tagged release.

### Changed

- **License: relicensed from MIT to AGPL-3.0-or-later.** pyflam is now
  distributed under the GNU Affero General Public License v3.0 or later. This is
  a strong network-copyleft license: anyone who distributes a modified version,
  **or offers it to users over a network** (e.g. hosting the Streamlit GUI in
  `pyflam_gui/` as a service), must make the complete corresponding source of
  their modified version available under the same terms.
  ([#2](https://github.com/crifod/pyflam/pull/2))
  - `LICENSE` now contains the verbatim GNU AGPL v3 text.
  - New `NOTICE` file records the copyright holder, the AGPL warranty notice, and
    the FlamMap independence disclaimer.
  - SPDX headers in all `src/` and `pyflam_gui/` source files updated to
    `AGPL-3.0-or-later`.
  - `pyproject.toml`, `README.md`, and the scientific-report docs updated to
    reflect the new license.

[Unreleased]: https://github.com/crifod/pyflam/compare/v0.1.3...HEAD
[0.1.3]: https://github.com/crifod/pyflam/compare/v0.1.2...v0.1.3
[0.1.2]: https://github.com/crifod/pyflam/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/crifod/pyflam/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/crifod/pyflam/releases/tag/v0.1.0
