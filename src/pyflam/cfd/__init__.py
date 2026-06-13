"""OpenFOAM-driven RANS terrain-wind solver for pyflam.

This subpackage builds, runs and reads an OpenFOAM atmospheric-boundary-layer
(ABL) CFD case to produce a momentum-conserving terrain wind field — the
high-fidelity counterpart to the fast diagnostic :mod:`pyflam.windsolver`.

It models, via ESI OpenFOAM's ``atmosphericModels`` suite:
    * momentum + k-epsilon turbulence (RANS),
    * Richards & Hoxey ABL inlet (``atmBoundaryLayerInlet*``) and rough-wall
      ground functions (``atmNutkWallFunction``),
    * optional non-neutral stability and diurnal slope flows via Boussinesq
      buoyancy (``buoyantBoussinesqSimpleFoam`` + ``atmBuoyancyTurbSource``) with
      a prescribed ground sensible-heat flux.

OpenFOAM is an external dependency, discovered at runtime through the ``openfoam``
wrapper (see :mod:`pyflam.cfd.run`); the pure-Python case/mesh/reader code is
usable and tested without it. The result is a :class:`pyflam.wind.WindField`,
identical in interface to the diagnostic solver's output.

Public API (added across implementation phases):
    solve_rans, wind_field_from_landscape
"""

from __future__ import annotations

import shutil
import tempfile

import numpy as np

from . import case, dictwriter, mesh, read, run
from .case import CaseConfig
from .mesh import build_terrain_mesh

__all__ = [
    "case", "dictwriter", "mesh", "read", "run",
    "CaseConfig", "solve_rans", "wind_field_from_landscape",
]


def solve_rans(
    elevation: np.ndarray,
    cellsize: float,
    *,
    speed: float,
    direction: float,
    z0=0.1,
    reference_height: float = 20.0,
    output_height: float = 6.1,
    buoyant: bool = False,
    surface_heat_flux: float = 0.0,
    reference_temperature: float = 300.0,
    nz: int = 20,
    domain_height: float | None = None,
    expansion_ratio: float = 1.2,
    iterations: int = 2000,
    n_processors: int = 1,
    init_potential: bool = False,
    west: float = 0.0,
    north: float = 0.0,
    crs: object | None = None,
    case_dir: str | None = None,
    keep_case: bool = False,
    openfoam: str | None = None,
):
    """Run an OpenFOAM RANS ABL simulation over a DEM and return a `WindField`.

    Builds a terrain-following mesh, writes the case, runs the solver
    (``simpleFoam`` neutral, or ``buoyantBoussinesqSimpleFoam`` when
    ``buoyant``), and reads the wind at ``output_height`` AGL. For ``buoyant``
    runs, ``surface_heat_flux`` (W/m^2) drives the buoyancy: a scalar gives a
    uniform flux (diurnal slope flows -- positive daytime/convective, negative
    nighttime/katabatic), or a **2D array on the DEM grid** gives a per-cell flux
    (a fire's convective heat -> a pyroconvective plume; see
    :mod:`pyflam.pyroconvection`).

    Requires OpenFOAM (discovered via the ``openfoam`` wrapper); raises
    ``FileNotFoundError`` if absent. The case is written to a temp dir (removed
    unless ``keep_case``) or to ``case_dir``.
    """
    elev = np.asarray(elevation, dtype=float)
    if domain_height is None:
        domain_height = max(3.0 * float(elev.max() - elev.min()), 500.0)

    m = build_terrain_mesh(elev, cellsize, nz=nz, domain_height=domain_height,
                           expansion_ratio=expansion_ratio)
    cfg = CaseConfig(
        speed=speed, direction=direction, z0=z0,
        reference_height=reference_height, buoyant=buoyant,
        surface_heat_flux=surface_heat_flux,
        reference_temperature=reference_temperature,
        iterations=iterations, n_processors=n_processors,
    )
    path = case_dir or tempfile.mkdtemp(prefix="pyflam_rans_")
    try:
        case.write_case(path, m, cfg)
        solver = "buoyantBoussinesqSimpleFoam" if buoyant else "simpleFoam"
        run.run_solver(path, solver=solver, n_processors=n_processors,
                       init_potential=init_potential, openfoam=openfoam)
        return read.read_wind_field(
            path, m, output_height=output_height, z0=z0,
            west=west, north=north, crs=crs)
    finally:
        if not keep_case and case_dir is None:
            shutil.rmtree(path, ignore_errors=True)


def wind_field_from_landscape(ls, *, speed: float, direction: float,
                              z0=None, per_cell_z0: bool = True, **kwargs):
    """RANS wind for a :class:`~pyflam.landscape.Landscape` (mirrors the
    diagnostic :func:`pyflam.windsolver.wind_field_from_landscape`).

    Pulls elevation and georeferencing from the landscape. If ``z0`` is omitted
    the roughness comes from the fuel grid: a **per-cell** field by default (the
    ground wall functions vary roughness cell by cell, the ABL inlet uses the
    median), or set ``per_cell_z0=False`` for a single median value. ``z0`` may
    also be passed explicitly as a scalar or a DEM-grid array.
    """
    if ls.elevation is None:
        raise ValueError("landscape has no elevation band for the RANS solver")
    if z0 is None:
        from ..windsolver import roughness_from_fuel
        rough = roughness_from_fuel(ls.fuel_model)
        z0 = np.asarray(rough, dtype=float) if per_cell_z0 \
            else float(np.median(rough))
    return solve_rans(
        np.asarray(ls.elevation, dtype=float), ls.cellsize_x,
        speed=speed, direction=direction, z0=z0,
        west=ls.west, north=ls.north, crs=ls.crs, **kwargs)
