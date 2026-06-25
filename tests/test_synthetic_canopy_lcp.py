"""Exercise the full crown pipeline on a synthetic canopy `.lcp` file.

Generates a synthetic landscape with the full canopy fuel stack, round-trips it
through the `.lcp` writer/reader, and runs crown classification, the crown-aware
spread field, and the plume-coupled crown march on the file-based landscape -- the
end-to-end path the bundled (canopy-band-less) Tuscany data can't drive.
"""

from __future__ import annotations

import sys

import numpy as np
import pytest

import pyflam
from pyflam import validate
from pyflam.units import mph_to_ft_per_min

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from make_synthetic_canopy_lcp import (  # noqa: E402
    RUN, make_synthetic_canopy_landscape,
)


def test_landscape_has_full_canopy_stack():
    ls = make_synthetic_canopy_landscape(80, seed=1)
    for band in ("canopy_cover", "canopy_base_height", "canopy_bulk_density",
                 "canopy_height"):
        a = getattr(ls, band)
        assert a is not None and np.any(np.asarray(a) > 0)


def test_lcp_round_trips_canopy_bands(tmp_path):
    ls = make_synthetic_canopy_landscape(80, seed=2)
    path = str(tmp_path / "canopy.lcp")
    ls.to_lcp(path)
    back = pyflam.Landscape.from_lcp(path)
    assert back.shape == ls.shape
    for band in ("canopy_base_height", "canopy_bulk_density", "canopy_height",
                 "canopy_cover", "fuel_model"):
        np.testing.assert_array_equal(
            np.asarray(getattr(back, band)), np.asarray(getattr(ls, band)))


def _crown(ls, model):
    return pyflam.crownfire.crown_fire_potential(
        ls, foliar_moisture=100.0, wind_20ft_ft_per_min=mph_to_ft_per_min(30),
        wind_midflame=pyflam.midflame_field(ls, mph_to_ft_per_min(30)),
        crown_spread=model, **RUN)


def test_synthetic_landscape_spans_surface_and_crown(tmp_path):
    ls = make_synthetic_canopy_landscape(100, seed=0)
    path = str(tmp_path / "c.lcp")
    ls.to_lcp(path)
    back = pyflam.Landscape.from_lcp(path)
    ft = _crown(back, "cruz2005")["fire_type"][back.fuel_model != 91]
    assert (ft == 0).any() and (ft == 2).any()        # both surface and active occur


def test_cruz_not_more_conservative_than_rothermel(tmp_path):
    """Cruz should classify at least as many active-crown cells as Rothermel
    (the operational stack under-predicts; Cruz 2005 corrects it)."""
    ls = make_synthetic_canopy_landscape(100, seed=0)
    burnable = ls.fuel_model != 91
    cruz = _crown(ls, "cruz2005")["fire_type"]
    roth = _crown(ls, "rothermel1991")["fire_type"]
    assert int((cruz[burnable] == 2).sum()) >= int((roth[burnable] == 2).sum())
    cmp = validate.compare_categories(cruz, roth, labels=[0, 1, 2], mask=burnable)
    assert 0.0 <= cmp.overall_agreement <= 1.0        # the harness runs on the file


def test_crown_march_runs_on_the_lcp(tmp_path):
    ls = make_synthetic_canopy_landscape(80, seed=3)
    path = str(tmp_path / "m.lcp")
    ls.to_lcp(path)
    back = pyflam.Landscape.from_lcp(path)
    caf = pyflam.crown_spread_field(
        back, wind_midflame=pyflam.midflame_field(back, mph_to_ft_per_min(30)),
        wind_direction=45.0, wind_20ft_ft_per_min=mph_to_ft_per_min(30),
        foliar_moisture=100.0, crown_spread="cruz2005", **RUN)
    r, c = np.unravel_index(int(np.argmax(caf.field.ros_max)), caf.field.shape)
    res = pyflam.fire_atmosphere_march(
        back, [(int(r), int(c))], total_time=30, dt=10, speed=8.0, direction=45.0,
        wind_provider=lambda l, i, a, s, d:
            pyflam.pyroconvection._uniform_wind_field(l, s, d),
        crown=True, foliar_moisture=100.0, max_wind_factor=4.0, **RUN)
    assert np.isfinite(res["arrival_time"]).sum() > 1
    assert "fire_type" in res
