"""pyflam_gui: a Streamlit GUI suite over the pyflam fire-behavior pipelines.

Four pages drive the existing pyflam pipelines from a map-driven UI:

1. Fire-weather / fire-danger / pyroconvection preview.
2. Real-time fire propagation (GFS / ICON-2I).
3. Reanalysis of a past event (ERA5).
4. Burn probability (drawn / random ignitions).

Launch with: ``streamlit run pyflam_gui/Home.py``.
"""

from __future__ import annotations

__version__ = "0.1.0"
