"""Page 3 — Reanalysis of a past fire event (ERA5).

Re-runs the propagation engine against ERA5 reanalysis at the event's historical
date to assess observed fire behavior and operational choices. Same panel as the
real-time page (:mod:`pyflam_gui.core.propagation`); only the weather source (past
real weather) and the default start time differ.
"""

from __future__ import annotations

from datetime import date

import pyflam_gui.core  # noqa: F401 -- path bootstrap for `pyflam`
import streamlit as st

from pyflam_gui.core import atmosphere as atmo
from pyflam_gui.core import propagation

st.set_page_config(page_title="Event Reanalysis", page_icon="🕰️", layout="wide")
st.title("Reanalysis of a past event (ERA5)")
st.info("Atmosphere is restricted to **ERA5 reanalysis** (past real weather). "
        "Needs `cdsapi` + CDS credentials (`~/.cdsapirc`) for fetch, or a "
        "downloaded ERA5 file. ERA5 surface heat fluxes are converted to W/m² up.")

propagation.render_page("reanalysis", allowed_sources=atmo.PAST_SOURCES,
                        default_day=date(2021, 7, 25), default_hour=12,
                        default_total=360)
