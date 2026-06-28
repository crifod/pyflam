"""Output-folder convention + writers shared by every page.

Each run gets a timestamped directory ``<base>/<kind>_<stamp>/`` into which the
page drops its rasters (GeoTIFF), vectors (GeoJSON), figures (PNG) and report
(PDF). ``zip_dir`` bundles it for a single Streamlit download button. The raster
and GeoJSON writers reuse pyflam's own I/O (``Landscape.to_geotiff``,
``OperativeReport.write_geojson``).
"""

from __future__ import annotations

import io
import os
import zipfile
from datetime import datetime


def run_stamp(when: datetime | None = None) -> str:
    """Filesystem-safe timestamp, e.g. ``20260628T143005``."""
    return (when or datetime.now()).strftime("%Y%m%dT%H%M%S")


def make_run_dir(base: str, kind: str, when: datetime | None = None) -> str:
    """Create and return ``<base>/<kind>_<stamp>/`` (parents created)."""
    path = os.path.join(base, f"{kind}_{run_stamp(when)}")
    os.makedirs(path, exist_ok=True)
    return path


def write_geotiff(ls, array, path: str, *, dtype: str = "float32", nodata=None):
    """Write one 2D array on the landscape grid to ``path`` (reuses pyflam I/O)."""
    import numpy as np
    kw = {"dtype": dtype}
    if nodata is not None:
        kw["nodata"] = nodata
    elif dtype.startswith("float"):
        kw["nodata"] = np.nan
    ls.to_geotiff(path, np.asarray(array), **kw)
    return path


def write_geojson(report, ls, path: str, **kwargs):
    """Write an :class:`OperativeReport` to GeoJSON (reuses pyflam I/O)."""
    report.write_geojson(ls, path, **kwargs)
    return path


def build_pdf(markdown_path: str, pdf_path: str, *, engine: str = "tectonic") -> str | None:
    """Render a Markdown report to PDF via pandoc; ``None`` (and a note) on failure.

    Mirrors the pattern in ``tests/pyroconv_daily.py``: pandoc + a self-contained
    LaTeX engine (tectonic / xelatex). Returns the PDF path, or ``None`` if pandoc
    or the engine is unavailable -- the PNG/GeoTIFF outputs are still written.
    """
    import subprocess
    try:
        subprocess.run(
            ["pandoc", markdown_path, "-o", pdf_path, f"--pdf-engine={engine}"],
            check=True, capture_output=True, timeout=300)
        return pdf_path
    except Exception:
        return None


def zip_dir(path: str) -> bytes:
    """Zip every file under ``path`` into an in-memory archive (for st.download_button)."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for root, _dirs, files in os.walk(path):
            for name in files:
                full = os.path.join(root, name)
                zf.write(full, os.path.relpath(full, path))
    buf.seek(0)
    return buf.getvalue()
