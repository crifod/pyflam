# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

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
