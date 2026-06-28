"""Session-state keys and typed accessors shared across the GUI pages.

``Home.py`` writes the user's choices (area of interest, landscape, atmosphere
source, output folder) into ``st.session_state``; every page reads them back
through the helpers here so the wiring is in one place and pages can show a clear
"set this up on the Home page first" message when a prerequisite is missing.
"""

from __future__ import annotations

from dataclasses import dataclass

import streamlit as st

# --- session_state keys -------------------------------------------------------
AOI = "pyflam_aoi"                 # AreaOfInterest | None
LANDSCAPE = "pyflam_landscape"     # pyflam.Landscape | None
LANDSCAPE_SRC = "pyflam_landscape_src"   # human-readable description
ATMO_CHOICE = "pyflam_atmo_choice"       # dict of atmosphere-picker selections
OUTPUT_BASE = "pyflam_output_base"       # str path


@dataclass
class AreaOfInterest:
    """A geographic bounding box selected on the map (WGS84 degrees)."""

    north: float
    west: float
    south: float
    east: float
    label: str = "custom"

    @property
    def bbox(self) -> tuple[float, float, float, float]:
        """``(north, west, south, east)`` -- the order ERA5 / fetch_era5 expects."""
        return (self.north, self.west, self.south, self.east)

    @property
    def lonlat_bbox(self) -> tuple[float, float, float, float]:
        """``(lon_min, lon_max, lat_min, lat_max)`` -- map/extent order."""
        return (self.west, self.east, self.south, self.north)

    @property
    def center(self) -> tuple[float, float]:
        """``(lat, lon)`` centre -- the point fed to single-column providers."""
        return (0.5 * (self.north + self.south), 0.5 * (self.west + self.east))


def get_aoi() -> AreaOfInterest | None:
    return st.session_state.get(AOI)


def set_aoi(aoi: AreaOfInterest | None) -> None:
    st.session_state[AOI] = aoi


def get_landscape():
    return st.session_state.get(LANDSCAPE)


def set_landscape(ls, src: str) -> None:
    st.session_state[LANDSCAPE] = ls
    st.session_state[LANDSCAPE_SRC] = src


def get_output_base() -> str:
    from . import REPO_ROOT
    import os
    return st.session_state.get(OUTPUT_BASE) or os.path.join(REPO_ROOT, "outputs")


def require(*, aoi: bool = False, landscape: bool = False) -> bool:
    """Guard a page: ``st.warning`` + ``False`` when a prerequisite is unmet."""
    ok = True
    if aoi and get_aoi() is None:
        st.warning("No area of interest yet -- pick one on the **Home** page first.")
        ok = False
    if landscape and get_landscape() is None:
        st.warning("No landscape loaded yet -- load one on the **Home** page first.")
        ok = False
    return ok
