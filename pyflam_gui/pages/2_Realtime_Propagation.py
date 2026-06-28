"""Page 2 — Real-time fire propagation (present-date GFS / ICON-2I).

Declare the run-duration scenario up front, place ignitions on the map, and spread
the fire through the evolving present-date weather with
:func:`pyflam.run_realtime` (``spatial=True``). Surfaces the whole energy-flux
dynamics over the declared window, the operative driving-force analysis, an
interactive perimeter map and a full output bundle. The shared body lives in
:mod:`pyflam_gui.core.propagation`.
"""

from __future__ import annotations

import pyflam_gui.core  # noqa: F401 -- path bootstrap for `pyflam`
import streamlit as st

from pyflam_gui.core import atmosphere as atmo
from pyflam_gui.core import propagation

st.set_page_config(page_title="Realtime Propagation", page_icon="🔥", layout="wide")
st.title("Real-time fire propagation")
st.info("Atmosphere is restricted to **present-date** sources (GFS / ICON-2I). "
        "For a past event use the **Event Reanalysis** page (ERA5).")

propagation.render_page("realtime", allowed_sources=atmo.PRESENT_SOURCES,
                        default_hour=10, default_total=240)
