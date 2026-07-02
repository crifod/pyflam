# JOSS submission — cover note & pre-submission checklist

Working notes for submitting pyflam to the **Journal of Open Source Software (JOSS)**.
Not part of the paper; a private aid for the submitting author.

---

## Suggested subject / editor track

- **Primary subject (domain / editor track):** computational wildland-fire-behavior
  modelling — **Earth & Environmental Sciences** (fire science / geoscientific
  modelling); scientific-computing / geospatial as secondary.
- **Lead theme to foreground:** an open, cross-platform operational wildfire-behavior
  toolkit whose distinctive contribution is **coupled fire–atmosphere pyroconvection and
  fire-danger assessment** — the defensible gap no operational open tool fills (routine,
  per-cell, forecast-driven danger mapping, incl. ICON-2I over Italy). The FlamMap-style
  product set is the validated credibility base, not the novelty.
- **Frontmatter tags already signal this:** `wildfire`, `fire behavior`,
  `fire–atmosphere coupling`, `pyroconvection`.
- **Suggested-reviewer expertise:** wildland-fire behavior; coupled fire–atmosphere /
  pyroconvection; scientific-Python / geospatial software.

---

## How JOSS submission works (mechanics)

1. Go to <https://joss.theoj.org> → **Submit a paper**.
2. Provide:
   - **Repository URL:** `https://github.com/crifod/pyflam`
   - **Branch (with the paper):** `main` (paper lives at `paper/paper.md`)
   - **Software version:** the released version under review (see checklist item on archiving)
   - **Archive DOI:** the Zenodo DOI (finalised at *acceptance* — see note below)
3. Paste the **cover note** (below) into the "comments to the editor" box.
4. An editor runs a pre-review, assigns reviewers, and review happens **openly in a
   GitHub issue** against a reviewer checklist. Reviewers install the software, run the
   tests, read the docs, and check the paper.
5. **At acceptance:** you make a final tagged release and a fresh Zenodo archive; the
   editor records that archive DOI as the version of record.

> **Archive-DOI note.** The existing archive `10.5281/zenodo.21129259` (v0.1.0) predates
> `paper/paper.md`. That's fine for *starting* the submission, but the **version of
> record is set at the end of review** — you will cut a new release (e.g. `v0.1.1` or
> `v0.2.0`) that includes the paper and any review changes, archive it, and give JOSS
> that DOI. Don't burn effort perfecting the archive before review.

---

## Cover note (paste into "comments to the editor")

> Dear editors,
>
> I am submitting **pyflam**, an open-source (AGPL-3.0), cross-platform Python library
> for wildland-fire-behavior modelling in support of fire management, planning, and
> suppression. It implements the Rothermel surface fire-spread model and the standard
> operational product set that desktop tools such as FlamMap established (surface rate
> of spread, fireline intensity, flame length, crown-fire potential, minimum-travel-time
> fire growth, and burn probability), and reads/writes the community `.lcp`, `.fms`,
> GeoTIFF and GeoJSON formats.
>
> Its distinctive contribution is to package this operational product set together with
> a weather-to-fire pipeline in one scriptable Python API, and to add — as transparent,
> selectable options next to the classical defaults — newer science that the closed
> desktop tools omit: weather-driven per-cell fuel-moisture conditioning, native
> terrain-wind solvers, an anisotropic-Eikonal (Finsler) spread solver alongside MTT, a
> literature-current crown-fire model, physics-based ember spotting, and coupled
> fire–atmosphere **pyroconvection and fire-danger assessment** driven by
> high-resolution regional forecasts (including the convection-permitting ICON-2I model
> over Italy).
>
> The deterministic surface core is cross-validated cell-by-cell against a real FlamMap
> run (rate of spread within ~3%, spread direction within ~1°, conditional fireline
> intensity within ~2%). The additional components are implemented and verified against
> their published equations and analytic benchmarks; broader validation against field
> observations and independent coupled models is described as ongoing future work. The
> software has a substantial automated test suite (~568 tests) run in CI on Python
> 3.11–3.13, with documentation, contribution and community guidelines, and an archived
> release with a DOI.
>
> I am the sole author and principal developer. The software has not been submitted
> elsewhere. I confirm the paper and metadata comply with JOSS requirements and that
> AI-assisted tooling was used in development, as disclosed in the repository.
>
> Suggested reviewer expertise: wildland-fire behavior modelling; scientific Python /
> geospatial software; coupled fire–atmosphere or numerical PDE methods. I have no
> conflicts to declare.
>
> Thank you for considering this submission.
>
> Cristiano Foderi (ORCID 0000-0002-5474-5433)

---

## Pre-submission checklist (mapped to JOSS's reviewer checklist)

### Software & repository
- [x] **Open-source, OSI-approved license** — AGPL-3.0-or-later (`LICENSE`, `NOTICE`).
- [x] **Public version-controlled repository** — `github.com/crifod/pyflam`.
- [x] **Substantial scholarly effort** — ~8,600 LOC across 25 modules; not a trivial
      utility. (JOSS looks for genuine research-software effort.)
- [x] **Author is a major contributor** — sole author/developer.
- [x] **Tagged release + archive with DOI** — `v0.1.0`, Zenodo `10.5281/zenodo.21129259`.
      ⚠️ Re-archive at acceptance to include the paper (see mechanics note).

### The paper (`paper/paper.md`) — current JOSS mandatory sections
- [x] **Summary** for a non-specialist audience.
- [x] **Statement of need** — who the audience is and what gap it fills.
- [x] **State of the field** — own section with the explicit build-vs-contribute
      justification, comparing to FlamMap, FARSITE, Cell2Fire, ELMFIRE, ForeFire,
      WindNinja, WRF-SFIRE.
- [x] **Software design** — trade-offs and architecture (SurfaceKernel, selectable
      back-ends, quasi-steady RANS coupling).
- [x] **Research impact statement** — near-term significance + realized use (daily
      ICON-2I operational mapping, the GUI, validated surface core).
- [x] **AI usage disclosure** — generative-AI assistance disclosed in the paper.
- [x] **Quality references** — `paper/paper.bib`, all entries verified (Crossref / FS
      Treesearch), full venue names, DOIs where they exist; includes related software.
- [x] **Length** — ~1,300 words (JOSS target 750–1750).
- [x] **Authors + affiliations + ORCID** present.
- [ ] **Final proofread** of the built PDF (Draft PDF workflow artifact).

### Documentation
- [x] **Statement of need** — also in `README.md`.
- [x] **Installation instructions** — `README.md` §Install (pip + extras).
- [x] **Example usage** — `README.md` §Quick start.
- [x] **Functionality documentation** — `README.md` + `docs/` + module docstrings.
- [x] **Automated tests** — 35 test files, ~568 tests, CI on 3 Python versions;
      how-to-run in `README.md` §Testing.
- [x] **Community guidelines** — `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, `SECURITY.md`.

### Integrity / disclosure
- [x] **AI-usage disclosure** — present in the report; note it in the cover letter.
- [x] **No competing prior submission**; no undisclosed conflicts.
- [ ] **Confirm claims are honestly scoped** — "verified vs validated" wording is
      deliberate; keep it (do not upgrade novel components to "validated").

### Optional strengthening before/ during review (not blockers)
- [ ] One worked, figure-bearing example (e.g. an ICON-2I pyroconvection/fire-danger map)
      in the docs, referenced from the paper.
- [ ] Consider whether to cite the ForeFire JOSS paper's model as the closest analogue.
- [ ] Add your institutional affiliation if it changes from "Independent researcher".

---

## Remaining bibliographic double-checks (author to confirm)
All five previously-uncertain refs were verified via Crossref / FS Treesearch and
corrected. No known outstanding errors. Re-confirm `castellnou2022` author list once
more against the published article if you want belt-and-suspenders certainty.
