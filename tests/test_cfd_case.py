"""Tests for OpenFOAM ABL case generation (no OpenFOAM)."""

from __future__ import annotations

import os

import numpy as np
import pytest

from pyflam.cfd import case as C
from pyflam.cfd import mesh as M


def _mesh():
    return M.build_terrain_mesh(np.full((4, 4), 1000.0), 30.0, nz=4,
                                domain_height=300.0)


def test_flow_direction():
    assert C.flow_direction(270.0) == pytest.approx((1.0, 0.0, 0.0), abs=1e-9)  # from W
    assert C.flow_direction(0.0) == pytest.approx((0.0, -1.0, 0.0), abs=1e-9)   # from N


def test_lateral_roles_cardinal_and_diagonal():
    west_wind = C.lateral_roles(C.flow_direction(270.0))  # blows east
    assert west_wind["west"] == "inlet" and west_wind["east"] == "outlet"
    assert west_wind["south"] == "slip" and west_wind["north"] == "slip"

    sw_wind = C.lateral_roles(C.flow_direction(225.0))     # from SW, blows NE
    assert sw_wind["west"] == "inlet" and sw_wind["south"] == "inlet"
    assert sw_wind["east"] == "outlet" and sw_wind["north"] == "outlet"


def test_neutral_case_files(tmp_path):
    C.write_case(str(tmp_path), _mesh(),
                 C.CaseConfig(speed=6.0, direction=270.0, z0=0.1))
    assert set(os.listdir(tmp_path / "0")) == {"U", "k", "epsilon", "nut", "p"}
    u = (tmp_path / "0" / "U").read_text()
    assert "atmBoundaryLayerInletVelocity" in u
    assert "Uref" in u and "6.0" in u
    # West patch (upwind) is the inlet; east is the outlet.
    assert "inletOutlet" in u
    ctrl = (tmp_path / "system" / "controlDict").read_text()
    assert "application     simpleFoam;" in ctrl
    assert "libatmosphericModels" in ctrl


def test_buoyant_case_adds_temperature_and_buoyancy(tmp_path):
    C.write_case(str(tmp_path), _mesh(),
                 C.CaseConfig(speed=6.0, direction=270.0, buoyant=True,
                              surface_heat_flux=50.0))
    fields = set(os.listdir(tmp_path / "0"))
    assert {"T", "p_rgh", "alphat"} <= fields and "p" not in fields
    t = (tmp_path / "0" / "T").read_text()
    assert "atmTurbulentHeatFluxTemperature" in t
    assert "q               uniform 50.0" in t       # diurnal heat flux
    assert (tmp_path / "constant" / "g").exists()
    fv = (tmp_path / "constant" / "fvOptions").read_text()
    assert "atmBuoyancyTurbSource" in fv
    ctrl = (tmp_path / "system" / "controlDict").read_text()
    assert "application     buoyantBoussinesqSimpleFoam;" in ctrl


def test_per_cell_z0_ground_nonuniform_inlet_uniform(tmp_path):
    """A per-cell z0 field -> nonuniform ground wall functions, uniform inlet."""
    z0 = np.full((4, 4), 0.1)
    z0[1, 2] = 0.9
    cfg = C.CaseConfig(speed=6.0, direction=270.0, z0=z0)
    assert C.representative_z0(cfg) == pytest.approx(0.1)   # median for the inlet
    C.write_case(str(tmp_path), _mesh(), cfg)
    nut = (tmp_path / "0" / "nut").read_text()
    eps = (tmp_path / "0" / "epsilon").read_text()
    u = (tmp_path / "0" / "U").read_text()
    assert "nonuniform List<scalar>" in nut and "0.9" in nut
    assert "nonuniform List<scalar>" in eps
    # The ABL inlet (in U) takes the single representative z0, not a list.
    assert "uniform 0.1" in u and "nonuniform" not in u


def test_scalar_z0_stays_uniform(tmp_path):
    C.write_case(str(tmp_path), _mesh(),
                 C.CaseConfig(speed=6.0, direction=270.0, z0=0.25))
    nut = (tmp_path / "0" / "nut").read_text()
    assert "uniform 0.25" in nut and "nonuniform" not in nut


def test_parallel_writes_decompose_dict(tmp_path):
    C.write_case(str(tmp_path), _mesh(),
                 C.CaseConfig(speed=5.0, direction=180.0, n_processors=4))
    dd = (tmp_path / "system" / "decomposeParDict").read_text()
    assert "numberOfSubdomains 4;" in dd
