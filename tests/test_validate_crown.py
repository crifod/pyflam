"""Tests for the crown-fire validation harness: the categorical (confusion-matrix)
comparison and its integration with crown_fire_potential."""

from __future__ import annotations

import numpy as np
import pytest

import pyflam
from pyflam import validate


# --- compare_categories -------------------------------------------------------

def test_perfect_agreement_is_diagonal():
    a = np.array([[0, 1, 2], [2, 1, 0]])
    cmp = validate.compare_categories(a, a, labels=[0, 1, 2])
    assert cmp.overall_agreement == 1.0
    assert cmp.n == 6
    assert np.array_equal(np.diag(cmp.matrix), [2, 2, 2])
    assert (cmp.matrix.sum() == np.trace(cmp.matrix))      # all on diagonal


def test_confusion_counts_and_recall():
    py = np.array([0, 0, 1, 1, 2, 2])
    ref = np.array([0, 1, 1, 2, 2, 2])     # 2 surface->{0,1}, passive->{1,2}, active->{2,2}
    cmp = validate.compare_categories(py, ref, labels=[0, 1, 2])
    # matrix[i, j] = #(pyflam==i & ref==j)
    assert cmp.matrix[0, 0] == 1 and cmp.matrix[0, 1] == 1
    assert cmp.matrix[1, 1] == 1 and cmp.matrix[1, 2] == 1
    assert cmp.matrix[2, 2] == 2
    assert cmp.overall_agreement == pytest.approx(4 / 6)
    assert cmp.recall(2) == pytest.approx(2 / 3)           # 2 of 3 ref-active matched


def test_mask_and_off_label_cells_ignored():
    py = np.array([[0, 1], [2, 9]])        # 9 is off-label
    ref = np.array([[0, 1], [2, 2]])
    mask = np.array([[True, True], [False, True]])         # exclude (1,0)
    cmp = validate.compare_categories(py, ref, labels=[0, 1, 2], mask=mask)
    assert cmp.n == 2                       # (0,0) and (1,1); (1,0) masked, (1,1) off-label
    assert cmp.overall_agreement == 1.0


def test_shape_mismatch_raises():
    with pytest.raises(ValueError):
        validate.compare_categories(np.zeros((2, 2)), np.zeros((2, 3)), labels=[0])


def test_summary_is_readable():
    cmp = validate.compare_categories(
        np.array([0, 1, 2]), np.array([0, 1, 1]), labels=[0, 1, 2])
    s = cmp.summary(names=["surface", "passive", "active"])
    assert "categorical agreement" in s and "confusion" in s


# --- integration with crown_fire_potential ------------------------------------

def _canopy_landscape(n=24):
    """A timber landscape with graded CBH/CBD spanning surface/passive/active."""
    cbh = np.tile(np.linspace(2.0, 80.0, n), (n, 1))         # CBH*10 -> 0.2..8 m
    cbd = np.tile(np.linspace(2.0, 40.0, n), (n, 1))         # CBD*100 -> 0.02..0.40
    return pyflam.Landscape(
        fuel_model=np.full((n, n), 10, dtype=int),
        slope=np.full((n, n), 25.0, dtype=float),
        aspect=np.full((n, n), 180.0, dtype=float),
        canopy_base_height=cbh, canopy_bulk_density=cbd,
        canopy_height=np.full((n, n), 150.0, dtype=float),
        cellsize_x=30.0, cellsize_y=30.0, west=0.0, north=n * 30.0,
        slope_units="degrees")


def test_crown_potential_round_trips_through_comparison():
    ls = _canopy_landscape()
    out = pyflam.crownfire.crown_fire_potential(
        ls, foliar_moisture=90.0, wind_20ft_ft_per_min=pyflam.units.mph_to_ft_per_min(20),
        wind_midflame=pyflam.midflame_field(ls, pyflam.units.mph_to_ft_per_min(20)),
        m_1h=0.04, m_10h=0.05, m_100h=0.06, m_live_herb=0.5, m_live_woody=0.7)
    ftype = np.asarray(out["fire_type"], dtype=np.int64)
    assert set(np.unique(ftype)).issubset({0, 1, 2})
    assert ftype.max() >= 1                              # some crowning is produced
    cmp = validate.compare_categories(ftype, ftype, labels=[0, 1, 2])
    assert cmp.overall_agreement == 1.0                 # self-diff is perfect
