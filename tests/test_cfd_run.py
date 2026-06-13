"""Tests for OpenFOAM discovery / orchestration helpers (no real solve)."""

from __future__ import annotations

import pytest

from pyflam.cfd import run


def test_find_explicit_and_env(tmp_path, monkeypatch):
    fake = tmp_path / "openfoam"
    fake.write_text("#!/bin/sh\n")
    assert run.find_openfoam(str(fake)) == str(fake)

    monkeypatch.setenv("PYFLAM_OPENFOAM", str(fake))
    assert run.find_openfoam() == str(fake)


def test_find_missing(monkeypatch):
    monkeypatch.delenv("PYFLAM_OPENFOAM", raising=False)
    monkeypatch.setattr(run.shutil, "which", lambda *_: None)
    assert run.find_openfoam() is None


def test_require_raises_when_missing(monkeypatch):
    monkeypatch.delenv("PYFLAM_OPENFOAM", raising=False)
    monkeypatch.setattr(run.shutil, "which", lambda *_: None)
    with pytest.raises(FileNotFoundError, match="OpenFOAM not found"):
        run.require_openfoam()


def test_quote_helper():
    assert run._q("/tmp/a b") == "'/tmp/a b'"
    assert run._q("it's") == "'it'\\''s'"


# --- end-to-end solve with per-cell roughness (skipped without OpenFOAM) ------

@pytest.mark.skipif(run.find_openfoam() is None, reason="needs OpenFOAM")
def test_per_cell_z0_roughness_slows_wind():
    """A rougher half-domain must slow the near-surface wind more than a smooth one."""
    import numpy as np
    import pyflam

    n, cs = 16, 30.0
    fuel = np.full((n, n), 104, dtype=int)   # grass (smooth)
    fuel[:, n // 2:] = 165                    # timber-shrub (rough) east half
    ls = pyflam.Landscape(
        fuel_model=fuel, slope=np.zeros((n, n)), aspect=np.zeros((n, n)),
        elevation=np.zeros((n, n)), cellsize_x=cs, cellsize_y=cs,
        west=0.0, north=n * cs, slope_units="degrees")
    wf = pyflam.cfd.wind_field_from_landscape(
        ls, speed=5.0, direction=270.0, iterations=800, nz=14)
    west_grass = wf.speed[:, :n // 2].mean()
    east_timber = wf.speed[:, n // 2:].mean()
    assert east_timber < west_grass           # rougher surface -> slower wind
