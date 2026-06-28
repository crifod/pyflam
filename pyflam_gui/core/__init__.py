"""Shared building blocks for the pyflam Streamlit GUI suite.

This package is the orchestration + presentation layer over the ``pyflam``
science library: an area-of-interest map, landscape/atmosphere loaders, an
output-folder writer and plotting helpers, all reused by the four pages in
``pyflam_gui/pages``.

Importing this package makes ``pyflam`` importable even when the project has not
been ``pip install``-ed, by adding the sibling ``src/`` directory to ``sys.path``
(the src-layout used by the repo). With an editable install this is a no-op.
"""

from __future__ import annotations

import os
import sys

# pyflam_gui/core/__init__.py -> repo root is two levels up from this file.
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_SRC = os.path.join(_REPO_ROOT, "src")
if os.path.isdir(_SRC) and _SRC not in sys.path:
    sys.path.insert(0, _SRC)

REPO_ROOT = _REPO_ROOT
