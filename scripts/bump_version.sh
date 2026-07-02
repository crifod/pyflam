#!/usr/bin/env bash
# Bump pyflam's version string in the three places it lives, in lockstep:
#   pyproject.toml, src/pyflam/__init__.py, pyflam_gui/__init__.py
#
# Usage:   scripts/bump_version.sh <new-version>
# Example: scripts/bump_version.sh 0.2.0
#          scripts/bump_version.sh 0.3.0.dev0
#
# Does NOT touch git — it only edits files, so it is safe to run for both a
# release version and a post-release .devN bump. Verifies all three afterwards.
set -euo pipefail

NEW="${1:-}"
if [[ -z "$NEW" ]]; then
  echo "usage: $0 <new-version>   (e.g. 0.2.0 or 0.3.0.dev0)" >&2
  exit 2
fi

# PEP 440-ish sanity check (release or .devN / rcN / aN / bN).
if ! [[ "$NEW" =~ ^[0-9]+\.[0-9]+\.[0-9]+([.-]?(dev|a|b|rc)[0-9]+)?$ ]]; then
  echo "error: '$NEW' does not look like a PEP 440 version" >&2
  exit 2
fi

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

FILES=(pyproject.toml src/pyflam/__init__.py pyflam_gui/__init__.py)

# Current version (from pyproject) for the log line.
OLD="$(grep -E '^version = ' pyproject.toml | head -1 | sed -E 's/^version = "(.*)"/\1/')"

# pyproject uses:  version = "X"     ;  packages use:  __version__ = "X"
perl -0pi -e "s/^version = \"[^\"]+\"/version = \"$NEW\"/m" pyproject.toml
perl -0pi -e "s/^__version__ = \"[^\"]+\"/__version__ = \"$NEW\"/m" \
  src/pyflam/__init__.py pyflam_gui/__init__.py

# Verify every file now carries exactly the new version.
fail=0
for f in "${FILES[@]}"; do
  if ! grep -q "\"$NEW\"" "$f"; then
    echo "ERROR: $f was not updated to $NEW" >&2
    fail=1
  fi
done
[[ $fail -eq 0 ]] || exit 1

echo "bumped version: ${OLD:-?} -> $NEW"
printf '  %s\n' "${FILES[@]}"
