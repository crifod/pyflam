"""Tests for the fire->atmosphere (pyroconvection) coupling.

The OpenFOAM solve itself is exercised separately (needs the solver); these tests
cover the pure-Python pieces: the heat-flux physics, the DEM->ground-face mapping,
and that a spatially-varying flux is written into the case as a nonuniform list.
"""

from __future__ import annotations

import os
import tempfile

import numpy as np
import pytest

import pyflam
from pyflam import pyroconvection
from pyflam.cfd import case as cfd_case
from pyflam.cfd.mesh import build_terrain_mesh


# --- heat-flux physics --------------------------------------------------------

def test_fire_heat_flux_formula():
    # q = chi_c * I[W/m] / cellsize; 1 Btu/ft/s = 3464.14 W/m.
    q = pyroconvection.fire_heat_flux(100.0, 100.0, convective_fraction=0.7)
    expected = 0.7 * 100.0 * 3464.14 / 100.0
    assert float(q) == pytest.approx(expected)


def test_fire_heat_flux_scales_and_masks():
    intensity = np.array([[0.0, 50.0], [200.0, 0.0]])
    q = pyroconvection.fire_heat_flux(intensity, 30.0)
    assert q.shape == intensity.shape
    assert q[1, 0] > q[0, 1] > 0.0          # more intensity -> more flux
    assert q[0, 0] == 0.0                    # no intensity -> no flux
    mask = np.array([[True, True], [False, True]])
    qm = pyroconvection.fire_heat_flux(intensity, 30.0, active_mask=mask)
    assert qm[1, 0] == 0.0                    # masked out despite intensity


def test_smaller_cells_give_higher_flux():
    """Same per-length intensity over a smaller cell -> higher area flux."""
    coarse = pyroconvection.fire_heat_flux(500.0, 100.0)
    fine = pyroconvection.fire_heat_flux(500.0, 10.0)
    assert fine == pytest.approx(10.0 * coarse)


# --- DEM -> ground-face mapping -----------------------------------------------

def test_ground_face_values_mapping():
    elev = np.zeros((3, 4))                   # ny=3, nx=4
    mesh = build_terrain_mesh(elev, 30.0, nz=3, domain_height=300.0)
    field = np.arange(12, dtype=float).reshape(3, 4)   # field[row, col]
    faces = cfd_case.ground_face_values(field, mesh)
    # Ground faces: i outer (0..nx-1), j inner (0..ny-1); (i,j)->DEM(ny-1-j, i).
    assert faces.size == mesh.nx * mesh.ny
    assert faces[0] == field[mesh.ny - 1, 0]            # i=0, j=0
    assert faces[1] == field[mesh.ny - 2, 0]            # i=0, j=1
    assert faces[mesh.ny] == field[mesh.ny - 1, 1]      # i=1, j=0


def test_ground_face_values_shape_check():
    mesh = build_terrain_mesh(np.zeros((3, 4)), 30.0, nz=2, domain_height=200.0)
    with pytest.raises(ValueError):
        cfd_case.ground_face_values(np.zeros((4, 4)), mesh)


# --- case writing with a per-cell heat flux -----------------------------------

def test_buoyant_case_writes_nonuniform_q():
    elev = np.zeros((4, 5))
    mesh = build_terrain_mesh(elev, 30.0, nz=4, domain_height=400.0)
    q = np.zeros((4, 5))
    q[1, 2] = 5000.0                          # one hot cell
    cfg = cfd_case.CaseConfig(speed=5.0, direction=270.0, buoyant=True,
                              surface_heat_flux=q)
    with tempfile.TemporaryDirectory() as d:
        cfd_case.write_case(d, mesh, cfg)
        with open(os.path.join(d, "0", "T")) as fh:
            text = fh.read()
    assert "nonuniform List<scalar>" in text
    assert f"\n{mesh.nx * mesh.ny}\n(" in text   # one value per ground face
    assert "5000.0" in text                      # the hot cell's flux appears


def test_scalar_heat_flux_stays_uniform():
    mesh = build_terrain_mesh(np.zeros((3, 3)), 30.0, nz=3, domain_height=300.0)
    cfg = cfd_case.CaseConfig(speed=5.0, direction=270.0, buoyant=True,
                              surface_heat_flux=150.0)
    with tempfile.TemporaryDirectory() as d:
        cfd_case.write_case(d, mesh, cfg)
        with open(os.path.join(d, "0", "T")) as fh:
            text = fh.read()
    assert "uniform 150.0" in text
    assert "nonuniform" not in text


def test_couple_fire_wind_requires_elevation():
    ls = pyflam.Landscape(
        fuel_model=np.full((4, 4), 102, dtype=int), slope=np.zeros((4, 4)),
        cellsize_x=30.0, cellsize_y=30.0, west=0.0, north=120.0)
    with pytest.raises(ValueError):
        pyroconvection.couple_fire_wind(ls, np.ones((4, 4)), speed=5.0,
                                        direction=270.0)


# --- time-marched coupling (mock wind provider, no OpenFOAM) ------------------

def _flat_landscape(n=51, fuel=102):
    return pyflam.Landscape(
        fuel_model=np.full((n, n), fuel, dtype=int), slope=np.zeros((n, n)),
        aspect=np.zeros((n, n)), elevation=np.zeros((n, n)),
        cellsize_x=30.0, cellsize_y=30.0, west=0.0, north=n * 30.0,
        slope_units="degrees")


_SC = dict(m_1h=0.06, m_10h=0.07, m_100h=0.08, m_live_herb=0.6, m_live_woody=0.9)


def _const_wind(ls, intensity, active, spd, dirn):
    return pyroconvection._uniform_wind_field(ls, spd, dirn)


def _plume_wind(ls, intensity, active, spd, dirn):
    wf = pyroconvection._uniform_wind_field(ls, spd, dirn)
    wf.speed = wf.speed.copy()
    wf.speed[active] = spd * 2.5            # plume accelerates wind at the front
    return wf


def test_march_runs_and_grows_outward():
    ls = _flat_landscape()
    res = pyflam.fire_atmosphere_march(
        ls, [(25, 25)], total_time=30, dt=10, speed=4.0, direction=270.0,
        wind_provider=_const_wind, **_SC)
    at = res["arrival_time"]
    assert at[25, 25] == 0.0
    assert np.isfinite(at).sum() > 1
    assert at[25, 26] < at[25, 35]          # arrival grows with distance


def test_march_history_length_matches_increments():
    ls = _flat_landscape()
    res = pyflam.fire_atmosphere_march(
        ls, [(25, 25)], total_time=40, dt=10, speed=4.0, direction=270.0,
        wind_provider=_const_wind, return_history=True, **_SC)
    assert res["times"] == [10, 20, 30, 40]
    assert len(res["winds"]) == len(res["fields"]) == 4


def test_plume_feedback_grows_a_larger_fire():
    ls = _flat_landscape()
    base = pyflam.fire_atmosphere_march(
        ls, [(25, 25)], total_time=40, dt=10, speed=4.0, direction=270.0,
        wind_provider=_const_wind, **_SC)
    coupled = pyflam.fire_atmosphere_march(
        ls, [(25, 25)], total_time=40, dt=10, speed=4.0, direction=270.0,
        wind_provider=_plume_wind, **_SC)
    assert np.isfinite(coupled["arrival_time"]).sum() > \
        np.isfinite(base["arrival_time"]).sum()


def test_march_constant_wind_matches_plain_mtt():
    """With a constant wind the march must equal a single uncoupled MTT run."""
    ls = _flat_landscape(n=41)
    res = pyflam.fire_atmosphere_march(
        ls, [(20, 20)], total_time=30, dt=10, speed=4.0, direction=270.0,
        wind_provider=_const_wind, **_SC)
    wf = pyroconvection._uniform_wind_field(ls, 4.0, 270.0)
    field = pyroconvection._field_from_wind(
        ls, wf, wind_reduction_factor=0.4, moist=_SC, load_factor=1.0)
    plain = pyflam.minimum_travel_time(field, [(20, 20)], max_time=30)
    both = np.isfinite(res["arrival_time"]) & np.isfinite(plain)
    assert np.array_equal(np.isfinite(res["arrival_time"]), np.isfinite(plain))
    assert np.allclose(res["arrival_time"][both], plain[both], atol=1e-6)


# --- atmosphere-driven march --------------------------------------------------

def test_march_from_atmosphere_provider():
    from pyflam.atmosphere import AtmosphericState, ConstantAtmosphere
    ls = _flat_landscape()
    st = AtmosphericState(wind_speed=6.0, wind_direction=270.0, temperature=30.0,
                          relative_humidity=20.0, sensible_heat_flux=200.0)
    res = pyflam.fire_atmosphere_march(
        ls, [(25, 25)], total_time=30, dt=10, atmosphere=ConstantAtmosphere(st),
        location=(43.0, 11.0), m_live_herb=0.6, m_live_woody=0.9,
        wind_provider=_const_wind)
    assert res["arrival_time"][25, 25] == 0.0
    assert np.isfinite(res["arrival_time"]).sum() > 1


def test_march_requires_inputs_or_atmosphere():
    ls = _flat_landscape(n=11)
    with pytest.raises(ValueError):
        pyflam.fire_atmosphere_march(ls, [(5, 5)], total_time=20, dt=10,
                                     wind_provider=_const_wind)


def test_march_evolving_weather_changes_wind():
    from datetime import datetime
    from pyflam.atmosphere import AtmosphericState, AtmosphereProvider

    class Ramp(AtmosphereProvider):
        def state_at(self, lat, lon, time):
            h = time.hour + time.minute / 60.0
            return AtmosphericState(wind_speed=2.0 + (h - 12.0),
                                    wind_direction=270.0, temperature=25.0,
                                    relative_humidity=30.0)

    ls = _flat_landscape()
    res = pyflam.fire_atmosphere_march(
        ls, [(25, 25)], total_time=60, dt=15, atmosphere=Ramp(),
        location=(43.0, 11.0), start_time=datetime(2024, 8, 1, 12, 0),
        m_live_herb=0.6, m_live_woody=0.9, wind_provider=_const_wind,
        return_history=True)
    winds = [float(w.speed.mean()) for w in res["winds"]]
    assert winds[-1] > winds[0]          # wind ramps up over the run


def test_spatial_march_per_cell_weather():
    from pyflam.atmosphere import AtmosphericState, ConstantAtmosphere
    ls = _flat_landscape()
    st = AtmosphericState(wind_speed=6.0, wind_direction=270.0, temperature=30.0,
                          relative_humidity=25.0, sensible_heat_flux=150.0)
    res = pyflam.fire_atmosphere_march(
        ls, [(25, 25)], total_time=40, dt=20, atmosphere=ConstantAtmosphere(st),
        spatial=True, m_live_herb=0.6, m_live_woody=0.9, return_history=True)
    assert res["arrival_time"][25, 25] == 0.0
    assert res["winds"][0].shape == ls.shape       # per-cell wind field


# --- plume / ambient superposition (no OpenFOAM) ------------------------------

def test_superpose_plume_adds_perturbation():
    from pyflam.atmosphere import AtmosphericState, wind_field_from_state
    from pyflam.pyroconvection import superpose_plume
    ls = _flat_landscape(n=8)
    amb = wind_field_from_state(AtmosphericState(5, 270, 25, 30), ls)
    base = wind_field_from_state(AtmosphericState(3, 270, 25, 30), ls)
    fire = wind_field_from_state(AtmosphericState(6, 270, 25, 30), ls)
    merged = superpose_plume(amb, base, fire)
    # ambient 5 + (fire 6 - base 3) = 8 m/s, same (westerly) direction.
    assert merged.speed[0, 0] == pytest.approx(8.0, abs=1e-6)
    assert merged.direction[0, 0] == pytest.approx(270.0, abs=1e-6)


def test_uv_windfield_roundtrip():
    from pyflam.atmosphere import AtmosphericState, wind_field_from_state
    from pyflam.pyroconvection import _uv_to_windfield, _windfield_to_uv
    ls = _flat_landscape(n=5)
    wf = wind_field_from_state(AtmosphericState(7.0, 225.0, 25, 30), ls)
    u, v = _windfield_to_uv(wf)
    back = _uv_to_windfield(u, v, wf)
    assert np.allclose(back.speed, wf.speed, atol=1e-6)
    assert np.allclose(back.direction, wf.direction, atol=1e-6)


# --- end-to-end OpenFOAM solve (skipped when the solver is absent) ------------

@pytest.mark.skipif(pyflam.cfd.run.find_openfoam() is None,
                    reason="needs OpenFOAM")
def test_coupled_solve_changes_wind():
    """A hot fire patch must measurably change the near-surface wind vs no fire."""
    n, cs = 14, 30.0
    ls = pyflam.Landscape(
        fuel_model=np.full((n, n), 104, dtype=int), slope=np.zeros((n, n)),
        aspect=np.zeros((n, n)), elevation=np.zeros((n, n)),
        cellsize_x=cs, cellsize_y=cs, west=0.0, north=n * cs,
        slope_units="degrees")
    intensity = np.zeros((n, n))
    intensity[6:9, 6:9] = 400.0          # intense central fire
    base = pyflam.cfd.wind_field_from_landscape(
        ls, speed=3.0, direction=270.0, iterations=800, nz=14)
    fire = pyflam.couple_fire_wind(
        ls, intensity, speed=3.0, direction=270.0, iterations=800, nz=14)
    delta = np.abs(fire.speed - base.speed)
    assert np.nanmax(delta) > 0.2        # the plume perturbs the wind
    assert not np.allclose(fire.speed, base.speed)


@pytest.mark.skipif(pyflam.cfd.run.find_openfoam() is None,
                    reason="needs OpenFOAM")
def test_marched_coupling_runs_with_real_cfd():
    """The full march drives the real buoyant CFD solver each increment."""
    n, cs = 12, 30.0
    ls = pyflam.Landscape(
        fuel_model=np.full((n, n), 104, dtype=int), slope=np.zeros((n, n)),
        aspect=np.zeros((n, n)), elevation=np.zeros((n, n)),
        cellsize_x=cs, cellsize_y=cs, west=0.0, north=n * cs,
        slope_units="degrees")
    res = pyflam.fire_atmosphere_march(
        ls, [(6, 6)], total_time=20, dt=10, speed=3.0, direction=270.0,
        load_factor=1.3, iterations=500, nz=12,   # -> couple_fire_wind / OpenFOAM
        m_1h=0.06, m_10h=0.07, m_100h=0.08, m_live_herb=0.6, m_live_woody=0.9)
    at = res["arrival_time"]
    assert at[6, 6] == 0.0
    assert np.isfinite(at).sum() > 1


@pytest.mark.skipif(pyflam.cfd.run.find_openfoam() is None,
                    reason="needs OpenFOAM")
def test_merge_plume_wind_keeps_gradient_and_perturbs():
    """Merged wind keeps the spatial ambient gradient AND adds a plume perturbation."""
    from pyflam.pyroconvection import merge_plume_wind
    n, cs = 14, 30.0
    ls = pyflam.Landscape(
        fuel_model=np.full((n, n), 104, dtype=int), slope=np.zeros((n, n)),
        aspect=np.zeros((n, n)), elevation=np.zeros((n, n)),
        cellsize_x=cs, cellsize_y=cs, west=0.0, north=n * cs, slope_units="degrees")
    spd = np.tile(np.linspace(3.0, 6.0, n), (n, 1))         # ramps W->E
    amb = pyflam.wind.WindField(speed=spd, direction=np.full((n, n), 270.0),
                                cellsize=cs, west=0.0, north=n * cs, speed_units="m/s")
    intensity = np.zeros((n, n))
    intensity[6:9, 6:9] = 400.0
    merged = merge_plume_wind(ls, amb, intensity, active_mask=intensity > 0,
                              iterations=400, nz=12)
    assert merged.speed[7, 0] != merged.speed[7, -1]        # ambient gradient kept
    assert np.nanmax(np.abs(merged.speed - amb.speed)) > 0.1  # plume perturbs it


@pytest.mark.skipif(pyflam.cfd.run.find_openfoam() is None,
                    reason="needs OpenFOAM")
def test_march_spatial_plume_merged():
    from pyflam.atmosphere import AtmosphericState, ConstantAtmosphere
    n, cs = 12, 30.0
    ls = pyflam.Landscape(
        fuel_model=np.full((n, n), 104, dtype=int), slope=np.zeros((n, n)),
        aspect=np.zeros((n, n)), elevation=np.zeros((n, n)),
        cellsize_x=cs, cellsize_y=cs, west=0.0, north=n * cs, slope_units="degrees")
    st = AtmosphericState(wind_speed=4.0, wind_direction=270.0, temperature=30.0,
                          relative_humidity=20.0, sensible_heat_flux=150.0)
    res = pyflam.fire_atmosphere_march(
        ls, [(6, 6)], total_time=20, dt=10, atmosphere=ConstantAtmosphere(st),
        spatial=True, plume=True, m_live_herb=0.6, m_live_woody=0.9,
        iterations=400, nz=12)
    assert res["arrival_time"][6, 6] == 0.0
    assert np.isfinite(res["arrival_time"]).sum() > 1
