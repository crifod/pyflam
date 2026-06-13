"""Generate an OpenFOAM ABL case (fields, constant, system) for a terrain mesh.

Two regimes, selected by ``CaseConfig.buoyant``:

* **Neutral** — ``simpleFoam``; fields ``U, p, k, epsilon, nut``. Richards & Hoxey
  ABL inlet (``atmBoundaryLayerInlet*``), rough-wall ground (``atmNutkWallFunction``).
* **Non-neutral / diurnal** — ``buoyantBoussinesqSimpleFoam``; adds ``T``,
  ``p_rgh``, ``alphat``, gravity and the ``atmBuoyancyTurbSource`` fvOption.
  Diurnal forcing is a prescribed ground sensible-heat flux ``q`` via
  ``atmTurbulentHeatFluxTemperature`` (q > 0 daytime/convective -> upslope;
  q < 0 nighttime/stable -> katabatic downslope).

Lateral patches are assigned inlet / outlet / slip per wind direction, so the
upwind face(s) get the ABL inlet and the rest let flow out or slip.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass

import numpy as np

from . import dictwriter
from .mesh import EAST, GROUND, LATERAL, NORTH, SOUTH, TOP, WEST, TerrainMesh

_ATM_LIB = '("libatmosphericModels.so")'


@dataclass
class CaseConfig:
    speed: float                       # input wind speed at reference_height (m/s)
    direction: float                   # direction wind blows FROM (deg from north)
    # surface roughness length (m): a scalar, or a 2D north-up array on the DEM
    # grid for per-cell roughness (e.g. from the fuel grid). The ABL inlet uses a
    # single representative value; the ground wall functions use the per-cell field.
    z0: "float | np.ndarray" = 0.1
    reference_height: float = 20.0     # Zref (m)
    buoyant: bool = False              # non-neutral (temperature + buoyancy)
    # ground sensible heat flux q (W/m^2): a scalar (uniform, e.g. diurnal) OR a
    # 2D north-up array on the DEM grid (a fire's per-cell convective flux ->
    # pyroconvective plume). See pyflam.pyroconvection.
    surface_heat_flux: "float | np.ndarray" = 0.0
    reference_temperature: float = 300.0   # TRef (K)
    iterations: int = 2000
    n_processors: int = 1
    nu: float = 1.5e-5                 # kinematic viscosity (m^2/s)
    beta: float = 3.0e-3              # thermal expansion (1/K)
    Pr: float = 0.9
    Prt: float = 0.74
    kappa: float = 0.4
    Cmu: float = 0.09


def flow_direction(direction_from: float) -> tuple[float, float, float]:
    """Unit wind velocity vector (east, north, up) for a from-direction."""
    phi = math.radians(direction_from)
    return (-math.sin(phi), -math.cos(phi), 0.0)


def lateral_roles(flow_dir, eps: float = 0.087) -> dict[str, str]:
    """Assign each lateral patch 'inlet' / 'outlet' / 'slip' from the flow dir.

    A patch whose outward normal opposes the flow is an inlet; aligned with it,
    an outlet; near-perpendicular (|component| < eps ~ 5 deg), a slip wall.
    """
    fx, fy, _ = flow_dir
    def role(comp):
        if comp > eps:
            return "inlet"
        if comp < -eps:
            return "outlet"
        return "slip"
    return {
        WEST: role(fx), EAST: role(-fx),
        SOUTH: role(fy), NORTH: role(-fy),
    }


# --- per-field boundary builders ---------------------------------------------

def representative_z0(cfg: CaseConfig) -> float:
    """Single roughness for the ABL inlet / initial profile (median if a field)."""
    z = cfg.z0
    return float(z) if np.isscalar(z) else float(np.median(np.asarray(z)))


def _z0_ground_entry(cfg: CaseConfig, mesh) -> str:
    """Ground-patch roughness: uniform scalar, or a per-face nonuniform list."""
    z = cfg.z0
    if np.isscalar(z):
        return f"uniform {float(z)}"
    vals = ground_face_values(np.asarray(z, dtype=float), mesh)
    body = "\n".join(repr(float(v)) for v in vals)
    return f"nonuniform List<scalar> \n{len(vals)}\n(\n{body}\n)"


def _abl_params(cfg: CaseConfig, flow_dir):
    return {
        "flowDir": f"({flow_dir[0]:.6f} {flow_dir[1]:.6f} 0)",
        "zDir": "(0 0 1)",
        "Uref": cfg.speed,
        "Zref": cfg.reference_height,
        "z0": f"uniform {representative_z0(cfg)}",
        "d": "uniform 0.0",
        "kappa": cfg.kappa,
        "Cmu": cfg.Cmu,
    }


def _u_field(cfg, mesh, flow_dir, roles, k_init):
    fx, fy, _ = flow_dir
    s = cfg.speed
    bf = {
        GROUND: {"type": "fixedValue", "value": "uniform (0 0 0)"},
        TOP: {"type": "slip"},
    }
    inlet = {"type": "atmBoundaryLayerInletVelocity", **_abl_params(cfg, flow_dir)}
    outlet = {"type": "inletOutlet", "inletValue": "uniform (0 0 0)",
              "value": "$internalField"}
    for p in LATERAL:
        bf[p] = dict(inlet) if roles[p] == "inlet" else (
            dict(outlet) if roles[p] == "outlet" else {"type": "slip"})
    return {
        "dimensions": "[0 1 -1 0 0 0 0]",
        "internalField": f"uniform ({fx * s:.4f} {fy * s:.4f} 0)",
        "boundaryField": bf,
    }


def _k_field(cfg, flow_dir, roles, k_init):
    bf = {
        GROUND: {"type": "kqRWallFunction", "value": f"uniform {k_init}"},
        TOP: {"type": "slip"},
    }
    inlet = {"type": "atmBoundaryLayerInletK", **_abl_params(cfg, flow_dir)}
    outlet = {"type": "inletOutlet", "inletValue": f"uniform {k_init}",
              "value": "$internalField"}
    for p in LATERAL:
        bf[p] = dict(inlet) if roles[p] == "inlet" else (
            dict(outlet) if roles[p] == "outlet" else {"type": "slip"})
    return {"dimensions": "[0 2 -2 0 0 0 0]",
            "internalField": f"uniform {k_init}", "boundaryField": bf}


def _epsilon_field(cfg, mesh, flow_dir, roles, eps_init):
    bf = {
        GROUND: {"type": "atmEpsilonWallFunction",
                 "z0": _z0_ground_entry(cfg, mesh),
                 "kappa": cfg.kappa, "Cmu": cfg.Cmu,
                 "value": f"uniform {eps_init}"},
        TOP: {"type": "slip"},
    }
    inlet = {"type": "atmBoundaryLayerInletEpsilon", **_abl_params(cfg, flow_dir)}
    outlet = {"type": "inletOutlet", "inletValue": f"uniform {eps_init}",
              "value": "$internalField"}
    for p in LATERAL:
        bf[p] = dict(inlet) if roles[p] == "inlet" else (
            dict(outlet) if roles[p] == "outlet" else {"type": "slip"})
    return {"dimensions": "[0 2 -3 0 0 0 0]",
            "internalField": f"uniform {eps_init}", "boundaryField": bf}


def _nut_field(cfg, mesh, roles):
    bf = {
        GROUND: {"type": "atmNutkWallFunction", "boundNut": False,
                 "z0": _z0_ground_entry(cfg, mesh), "value": "uniform 0"},
        TOP: {"type": "calculated", "value": "uniform 0"},
    }
    for p in LATERAL:
        bf[p] = ({"type": "slip"} if roles[p] == "slip"
                 else {"type": "calculated", "value": "uniform 0"})
    return {"dimensions": "[0 2 -1 0 0 0 0]", "internalField": "uniform 0",
            "boundaryField": bf}


def _p_field(roles, name="p"):
    bf = {GROUND: {"type": "zeroGradient"}, TOP: {"type": "zeroGradient"}}
    for p in LATERAL:
        if roles[p] == "outlet":
            bf[p] = {"type": "fixedValue", "value": "uniform 0"}
        elif roles[p] == "slip":
            bf[p] = {"type": "slip"}
        else:
            bf[p] = {"type": "zeroGradient"}
    return {"dimensions": "[0 2 -2 0 0 0 0]", "internalField": "uniform 0",
            "boundaryField": bf}


def _prgh_field(roles):
    bf = {GROUND: {"type": "fixedFluxPressure", "value": "uniform 0"},
          TOP: {"type": "fixedFluxPressure", "value": "uniform 0"}}
    for p in LATERAL:
        if roles[p] == "outlet":
            bf[p] = {"type": "fixedValue", "value": "uniform 0"}
        elif roles[p] == "slip":
            bf[p] = {"type": "slip"}
        else:
            bf[p] = {"type": "fixedFluxPressure", "value": "uniform 0"}
    return {"dimensions": "[0 2 -2 0 0 0 0]", "internalField": "uniform 0",
            "boundaryField": bf}


def ground_face_values(field: np.ndarray, mesh: TerrainMesh) -> np.ndarray:
    """Map a north-up DEM-grid field to the GROUND patch's face order.

    Ground faces are emitted with mesh column ``i`` outer, ``j`` inner (see
    :func:`pyflam.cfd.mesh.build_terrain_mesh`), and mesh ``(i, j)`` maps to DEM
    ``(row=ny-1-j, col=i)``. Returns one value per ground face.
    """
    arr = np.asarray(field, dtype=float)
    if arr.shape != (mesh.ny, mesh.nx):
        raise ValueError(
            f"heat-flux field {arr.shape} != DEM grid {(mesh.ny, mesh.nx)}")
    ii, jj = np.meshgrid(np.arange(mesh.nx), np.arange(mesh.ny), indexing="ij")
    ii, jj = ii.ravel(), jj.ravel()           # i outer, j inner (face order)
    return arr[mesh.ny - 1 - jj, ii]


def _q_entry(cfg, mesh) -> str:
    """The ``q`` boundary value: uniform scalar, or a per-face nonuniform list."""
    shf = cfg.surface_heat_flux
    if np.isscalar(shf):
        return f"uniform {float(shf)}"
    vals = ground_face_values(np.asarray(shf, dtype=float), mesh)
    body = "\n".join(repr(float(v)) for v in vals)
    return f"nonuniform List<scalar> \n{len(vals)}\n(\n{body}\n)"


def _t_field(cfg, mesh, roles):
    tref = cfg.reference_temperature
    bf = {
        GROUND: {"type": "atmTurbulentHeatFluxTemperature", "heatSource": "flux",
                 "alphaEff": "alphaEff", "Cp0": 1005.0,
                 "q": _q_entry(cfg, mesh),
                 "value": f"uniform {tref}"},
        TOP: {"type": "fixedValue", "value": f"uniform {tref}"},
    }
    for p in LATERAL:
        if roles[p] == "inlet":
            bf[p] = {"type": "fixedValue", "value": f"uniform {tref}"}
        elif roles[p] == "outlet":
            bf[p] = {"type": "inletOutlet", "inletValue": f"uniform {tref}",
                     "value": f"uniform {tref}"}
        else:
            bf[p] = {"type": "zeroGradient"}
    return {"dimensions": "[0 0 0 1 0 0 0]",
            "internalField": f"uniform {tref}", "boundaryField": bf}


def _alphat_field(cfg, mesh, roles):
    bf = {
        GROUND: {"type": "atmAlphatkWallFunction", "Cmu": cfg.Cmu,
                 "kappa": cfg.kappa, "Pr": cfg.Pr,
                 "z0": _z0_ground_entry(cfg, mesh),
                 "Prt": f"uniform {cfg.Prt}", "value": "uniform 0"},
        TOP: {"type": "calculated", "value": "uniform 0"},
    }
    for p in LATERAL:
        bf[p] = ({"type": "slip"} if roles[p] == "slip"
                 else {"type": "calculated", "value": "uniform 0"})
    return {"dimensions": "[0 2 -1 0 0 0 0]", "internalField": "uniform 0",
            "boundaryField": bf}


# --- top-level case writer ----------------------------------------------------

def write_case(case_dir: str, mesh: TerrainMesh, cfg: CaseConfig) -> None:
    """Write a complete OpenFOAM case (mesh + 0/ + constant/ + system/)."""
    for sub in ("0", "constant", "system"):
        os.makedirs(os.path.join(case_dir, sub), exist_ok=True)
    mesh.write(case_dir)

    flow_dir = flow_direction(cfg.direction)
    roles = lateral_roles(flow_dir)
    z0_ref = representative_z0(cfg)
    ustar = cfg.kappa * cfg.speed / math.log(
        (cfg.reference_height + z0_ref) / z0_ref)
    k_init = round(ustar ** 2 / math.sqrt(cfg.Cmu), 5)
    eps_init = round(ustar ** 3 / (cfg.kappa * (cfg.reference_height + z0_ref)), 6)

    def w(sub, cls, obj, entries):
        dictwriter.write(os.path.join(case_dir, sub, obj), cls, obj, entries,
                         location=sub)

    # 0/ fields
    w("0", "volVectorField", "U", _u_field(cfg, mesh, flow_dir, roles, k_init))
    w("0", "volScalarField", "k", _k_field(cfg, flow_dir, roles, k_init))
    w("0", "volScalarField", "epsilon",
      _epsilon_field(cfg, mesh, flow_dir, roles, eps_init))
    w("0", "volScalarField", "nut", _nut_field(cfg, mesh, roles))
    if cfg.buoyant:
        w("0", "volScalarField", "p_rgh", _prgh_field(roles))
        w("0", "volScalarField", "T", _t_field(cfg, mesh, roles))
        w("0", "volScalarField", "alphat", _alphat_field(cfg, mesh, roles))
    else:
        w("0", "volScalarField", "p", _p_field(roles))

    _write_constant(case_dir, cfg, w)
    _write_system(case_dir, cfg, w)


def _write_constant(case_dir, cfg: CaseConfig, w):
    w("constant", "dictionary", "turbulenceProperties", {
        "simulationType": "RAS",
        "RAS": {"RASModel": "kEpsilon", "turbulence": "on", "printCoeffs": "on"},
    })
    if cfg.buoyant:
        w("constant", "dictionary", "transportProperties", {
            "transportModel": "Newtonian",
            "nu": f"[0 2 -1 0 0 0 0] {cfg.nu}",
            "beta": f"[0 0 0 -1 0 0 0] {cfg.beta}",
            "TRef": f"[0 0 0 1 0 0 0] {cfg.reference_temperature}",
            "Pr": f"[0 0 0 0 0 0 0] {cfg.Pr}",
            "Prt": f"[0 0 0 0 0 0 0] {cfg.Prt}",
        })
        w("constant", "uniformDimensionedVectorField", "g", {
            "dimensions": "[0 1 -2 0 0 0 0]", "value": "(0 0 -9.81)"})
        w("constant", "dictionary", "fvOptions", {
            "atmAmbientTurbSource1": {
                "type": "atmAmbientTurbSource", "selectionMode": "all",
                "kAmb": 1.0e-04, "epsilonAmb": 7.208e-08},
            "atmBuoyancyTurbSource1": {
                "type": "atmBuoyancyTurbSource", "selectionMode": "all",
                "rho": "rhok", "Lmax": 41.575, "beta": cfg.beta},
        })
    else:
        w("constant", "dictionary", "transportProperties", {
            "transportModel": "Newtonian",
            "nu": f"[0 2 -1 0 0 0 0] {cfg.nu}"})


def _write_system(case_dir, cfg: CaseConfig, w):
    solver = "buoyantBoussinesqSimpleFoam" if cfg.buoyant else "simpleFoam"
    w("system", "dictionary", "controlDict", {
        "application": solver,
        "startFrom": "startTime", "startTime": 0,
        "stopAt": "endTime", "endTime": cfg.iterations,
        "deltaT": 1, "writeControl": "timeStep",
        "writeInterval": cfg.iterations, "purgeWrite": 0,
        "writeFormat": "ascii", "writePrecision": 7, "runTimeModifiable": "true",
        "libs": _ATM_LIB,
    })

    grad = "Gauss linear"
    div = {"default": "none",
           "div(phi,U)": "bounded Gauss linearUpwind grad(U)",
           "div(phi,k)": "bounded Gauss limitedLinear 1",
           "div(phi,epsilon)": "bounded Gauss limitedLinear 1",
           "div((nuEff*dev2(T(grad(U)))))": "Gauss linear",
           "div(phi,T)": "bounded Gauss limitedLinear 1",
           "div((nuEff*dev(2*symm(grad(U)))))": "Gauss linear"}
    w("system", "dictionary", "fvSchemes", {
        "ddtSchemes": {"default": "steadyState"},
        "gradSchemes": {"default": grad},
        "divSchemes": div,
        "laplacianSchemes": {"default": "Gauss linear corrected"},
        "interpolationSchemes": {"default": "linear"},
        "snGradSchemes": {"default": "corrected"},
    })

    fields = ["p", "U", "k", "epsilon"]
    if cfg.buoyant:
        fields = ["p_rgh", "U", "T", "k", "epsilon"]
    solvers = {}
    for f in fields:
        if f in ("p", "p_rgh"):
            solvers[f] = {"solver": "GAMG", "smoother": "GaussSeidel",
                          "tolerance": 1e-7, "relTol": 0.01}
        else:
            solvers[f] = {"solver": "smoothSolver", "smoother": "symGaussSeidel",
                          "tolerance": 1e-7, "relTol": 0.1}
    # Buoyant runs are far less stable over terrain -> stronger under-relaxation.
    if cfg.buoyant:
        p_field, p_relax = "p_rgh", 0.3
        relax_fields = {"U": 0.3, "k": 0.3, "epsilon": 0.3, "T": 0.5}
    else:
        p_field, p_relax = "p", 0.3
        relax_fields = {"U": 0.7, "k": 0.7, "epsilon": 0.7}
    w("system", "dictionary", "fvSolution", {
        "solvers": solvers,
        "SIMPLE": {"nNonOrthogonalCorrectors": 1, "pRefCell": 0, "pRefValue": 0,
                   "residualControl": {p_field: 1e-4, "U": 1e-4, "k": 1e-4,
                                       "epsilon": 1e-4}},
        "relaxationFactors": {
            "fields": {p_field: p_relax},
            "equations": relax_fields},
    })

    n = cfg.n_processors
    if n > 1:
        w("system", "dictionary", "decomposeParDict", {
            "numberOfSubdomains": n, "method": "scotch"})
