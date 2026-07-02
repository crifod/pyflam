---
description: Cut a pyflam release (version bump + CHANGELOG + tag + GitHub release + Zenodo archive)
argument-hint: <version>  e.g. 0.2.0  (omit to be asked)
allowed-tools: Bash, Edit, Read
---

You are cutting a pyflam release. Target version: **$1** (if empty, ask the user for the
SemVer version first, and help them choose: PATCH = packaging/docs/bugfix only, MINOR =
new backwards-compatible features, MAJOR = breaking changes).

Follow this checklist exactly, pausing only where noted. Do not skip the CI gate.

## 0. Preflight
- Confirm on `main`, clean working tree, and up to date: `git checkout main && git pull`.
- Confirm the version argument is a release version (not `.devN`). Read the current
  version from `pyproject.toml`.
- Confirm `CHANGELOG.md` has content under `## [Unreleased]` (if it's empty, ask the
  user what changed — a release with no changelog entries is a smell).

## 1. Release branch + version bump
- `git checkout -b release-$1`
- Run `scripts/bump_version.sh $1` (edits pyproject.toml + both `__init__.py`, verifies).

## 2. CHANGELOG
- Promote `## [Unreleased]` → add a new dated section `## [$1] - <today>` beneath a fresh
  empty `## [Unreleased]`, moving the accumulated entries into `[$1]`.
- Update the reference links at the bottom:
  - `[Unreleased]: .../compare/v$1...HEAD`
  - add `[$1]: .../compare/v<prev>...v$1`
- Use today's date (get it with `date +%Y-%m-%d`).

## 3. Verify + PR
- Sanity: `PYTHONPATH=src python -c "import pyflam, pyflam_gui; print(pyflam.__version__)"`
  must print `$1`.
- Commit (`release: v$1`), push, open a PR titled `release: v$1`.
- **Gate:** watch CI with `gh pr checks <n> --watch` — all green before merging.
- Merge with `--merge --delete-branch`; `git checkout main && git pull`.

## 4. Tag + GitHub release
- Annotated tag on the merge commit: `git tag -a v$1 -m "pyflam v$1" HEAD && git push origin v$1`.
- `gh release create v$1 --title "pyflam v$1" --verify-tag --notes "<from the [$1] CHANGELOG section>"`.

## 5. Zenodo (automatic) + confirm
- The GitHub–Zenodo webhook is enabled, so the release auto-archives under concept DOI
  **10.5281/zenodo.21129259**. Tell the user to confirm the new version record at
  <https://zenodo.org/me/uploads> and that its metadata mapped from `.zenodo.json`.
- The concept DOI is unchanged, so `CITATION.cff` / README need no DOI edit. Do bump
  `version:` and `date-released:` in `CITATION.cff` to match `$1`.

## 6. Reopen the dev cycle
- New branch `post-release-dev-bump`; run `scripts/bump_version.sh <next>.dev0`
  (next minor, e.g. after 0.2.0 → 0.3.0.dev0; confirm the target with the user).
- Commit (`chore: reopen dev cycle at <next>.dev0`), PR, CI-green, merge.

## Done
Report: released version, tag URL, and the two follow-ups for the user (verify Zenodo
record; submit/refresh any JOSS archive DOI if a paper submission is in flight).
