# SPDX-License-Identifier: AGPL-3.0-or-later
"""Paths and loaders for the vendored pyroconvection validation corpus.

Everything the §21 analysis chain needs lives in ``docs/pyroconv_validation/`` so the results
can be reproduced from a clone alone -- no scratchpad, no ``/tmp``, no re-download. Override the
root with ``PYFLAM_VALDATA`` to point at a different corpus.

Contents:

* ``portal_fires.json``    fire behaviour from the Wildfire Data Portal REST catalogue
* ``portal_all.json``      the same fires joined to the portal's fire-classification taxonomy
* ``fire_labels.json``     observed pyroconvection class per fire
* ``sonde_times.json``     sonde launch datetimes (Castellnou Ribau et al. 2025, Table S1)
* ``oos_fires.json``       the 8 labelled fires with no sounding
* ``oos_margin.json``      their computed margins (regenerate with scripts/oos_margin.py)
* ``era5_oos_index.json``  slug -> ERA5 file and hour for those 8
* ``sondes/``              27 ambient sonde profiles, one per labelled sonde-fire pair
* ``era5/``                ERA5 pressure-level profiles at sonde launch (era5s_) and at the
                           out-of-sample fires' peak-growth hour (era5o_)
* ``perims/``              perimeter KMZs, the source of every growth rate

See ``ATTRIBUTION.md`` in that directory for provenance and licensing of each source.
"""
from __future__ import annotations

import glob
import json
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.environ.get("PYFLAM_VALDATA", os.path.join(REPO, "docs", "pyroconv_validation"))

SONDES = os.path.join(DATA, "sondes")
ERA5 = os.path.join(DATA, "era5")
PERIMS = os.path.join(DATA, "perims")


def load(name: str):
    """Read one of the vendored JSON inputs by stem, e.g. ``load("fire_labels")``."""
    path = os.path.join(DATA, f"{name}.json")
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{path} missing. The validation corpus lives in docs/pyroconv_validation/; "
            f"set PYFLAM_VALDATA to override its location.")
    with open(path) as fh:
        return json.load(fh)


def portal_by_slug() -> dict:
    """Portal fire behaviour keyed by slug."""
    return {p["slug"]: p for p in load("portal_fires")}


def labels_by_slug() -> dict:
    """Observed pyroconvection class name keyed by slug."""
    return {f["slug"]: f["obs_class"] for f in load("fire_labels")}


def sonde_paths() -> list:
    """Every vendored ambient sonde CSV, sorted."""
    return sorted(glob.glob(os.path.join(SONDES, "*.csv")))


def era5_for(slug: str, when, prefix: str = "era5s") -> str | None:
    """ERA5 profile for ``slug`` nearest ``when``; None if the corpus has no match.

    Prefers the exact hour, then the closest hour on the same day -- the sonde-launch batch is
    keyed to the minute-rounded hour, which does not always match a caller's own rounding.
    """
    exact = os.path.join(ERA5, f"{prefix}_{slug}_{when:%Y%m%d}_{when.hour:02d}.nc")
    if os.path.exists(exact):
        return exact
    same_day = sorted(glob.glob(os.path.join(ERA5, f"era5*_{slug}_{when:%Y%m%d}_*.nc")))
    if not same_day:
        return None
    return min(same_day, key=lambda c: abs(int(c[-5:-3]) - when.hour))
