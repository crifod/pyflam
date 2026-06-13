"""Locate OpenFOAM and orchestrate the solver pipeline.

OpenFOAM is invoked through the ``openfoam`` wrapper (the native macOS app and
most installs expose it): ``openfoam -c 'cd <case> && <cmd>'`` runs a command
inside the OpenFOAM environment. Discovery order: an explicit path argument, the
``PYFLAM_OPENFOAM`` env var, then ``openfoam`` on ``PATH``. A clear
:class:`FileNotFoundError` is raised when none is found.
"""

from __future__ import annotations

import os
import shutil
import subprocess


def find_openfoam(explicit: str | None = None) -> str | None:
    """Locate the ``openfoam`` wrapper executable, or return ``None``."""
    for cand in (explicit, os.environ.get("PYFLAM_OPENFOAM"), "openfoam"):
        if not cand:
            continue
        path = cand if os.path.exists(cand) else shutil.which(cand)
        if path:
            return path
    return None


def require_openfoam(explicit: str | None = None) -> str:
    of = find_openfoam(explicit)
    if of is None:
        raise FileNotFoundError(
            "OpenFOAM not found. Install OpenFOAM (e.g. `brew install "
            "gerlero/openfoam/openfoam`) and ensure the `openfoam` wrapper is on "
            "PATH, or set PYFLAM_OPENFOAM."
        )
    return of


def foam_run(case_dir: str, command: str, *, openfoam: str,
             check: bool = True, capture: bool = False) -> subprocess.CompletedProcess:
    """Run one OpenFOAM command inside ``case_dir`` via the wrapper."""
    inner = f"cd {_q(case_dir)} && {command}"
    return subprocess.run(
        [openfoam, "-c", inner], check=check,
        capture_output=capture, text=True,
    )


def check_mesh(case_dir: str, *, openfoam: str | None = None) -> str:
    """Run ``checkMesh`` and return its output."""
    of = require_openfoam(openfoam)
    return foam_run(case_dir, "checkMesh -constant", openfoam=of,
                    capture=True).stdout


def run_solver(
    case_dir: str,
    *,
    solver: str,
    n_processors: int = 1,
    init_potential: bool = False,
    openfoam: str | None = None,
) -> None:
    """Run the (potentialFoam) -> (decompose) -> solve -> reconstruct pipeline.

    Note: ``renumberMesh`` is intentionally *not* run — the reader relies on the
    structured cell ordering (``cidx``) to map the solved field back to columns.
    """
    of = require_openfoam(openfoam)
    if init_potential:
        foam_run(case_dir, "potentialFoam -initialiseUBCs", openfoam=of, check=False)
    if n_processors > 1:
        foam_run(case_dir, "decomposePar -force", openfoam=of)
        foam_run(case_dir, f"mpirun -np {n_processors} {solver} -parallel",
                 openfoam=of)
        foam_run(case_dir, "reconstructPar -latestTime", openfoam=of)
    else:
        foam_run(case_dir, solver, openfoam=of)


def _q(path: str) -> str:
    return "'" + path.replace("'", "'\\''") + "'"
