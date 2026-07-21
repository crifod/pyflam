# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
