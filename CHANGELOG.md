# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

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
