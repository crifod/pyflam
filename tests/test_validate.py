"""Tests for the FlamMap-comparison machinery (synthetic data, no big rasters)."""

from __future__ import annotations

import numpy as np
import pytest

from pyflam import validate
from pyflam.units import (
    chains_per_hour_to_ft_per_min,
    ft_per_min_to_chains_per_hour,
)


def test_chains_per_hour_roundtrip():
    assert chains_per_hour_to_ft_per_min(60.0) == pytest.approx(66.0)  # 1 chain/min
    assert ft_per_min_to_chains_per_hour(66.0) == pytest.approx(60.0)
    for r in (0.0, 3.3, 100.0, 1713.8):
        assert ft_per_min_to_chains_per_hour(
            chains_per_hour_to_ft_per_min(r)) == pytest.approx(r)


def test_identical_fields_score_perfectly():
    rng = np.random.default_rng(0)
    a = rng.uniform(1, 100, size=(20, 20))
    c = validate.compare_fields(a, a.copy())
    assert c.n == 400
    assert c.bias == pytest.approx(0.0, abs=1e-12)
    assert c.rmse == pytest.approx(0.0, abs=1e-12)
    assert c.median_ratio == pytest.approx(1.0)
    assert c.pearson_r == pytest.approx(1.0)
    assert c.log_rmse == pytest.approx(0.0, abs=1e-12)
    assert c.within_10pct == 1.0
    assert c.fit_slope == pytest.approx(1.0)


def test_scaled_field_detected():
    a = np.full((10, 10), 10.0)
    b = np.full((10, 10), 5.0)          # pyflam is 2x the reference
    c = validate.compare_fields(a, b)
    assert c.median_ratio == pytest.approx(2.0)
    assert c.bias == pytest.approx(5.0)
    assert c.log_rmse == pytest.approx(np.log10(2.0))
    assert c.within_25pct == 0.0


def test_classification_counts():
    # 2x2: agree-burn, agree-zero, pyflam-only, reference-only
    a = np.array([[10.0, 0.0], [5.0, 0.0]])
    b = np.array([[8.0, 0.0], [0.0, 4.0]])
    c = validate.compare_fields(a, b)
    assert c.n_mask == 4
    assert c.both_burn == 1
    assert c.both_zero == 1
    assert c.pyflam_only == 1
    assert c.reference_only == 1
    assert c.n == 1                      # only the agree-burn cell is in numerics


def test_mask_excludes_nodata():
    a = np.array([[10.0, np.nan], [20.0, 30.0]])
    b = np.array([[10.0, 5.0], [np.nan, 30.0]])
    mask = np.ones((2, 2), dtype=bool)
    c = validate.compare_fields(a, b, mask=mask)
    # The two NaN cells drop out; only (0,0) and (1,1) remain, both exact.
    assert c.n == 2
    assert c.rmse == pytest.approx(0.0)


def test_scan_parameter_finds_best():
    ref = np.full((8, 8), 50.0)
    # model(v) = v * 10; the reference equals model(5).
    best, best_cmp, table = validate.scan_parameter(
        lambda v: np.full((8, 8), v * 10.0), [3, 4, 5, 6, 7], ref,
    )
    assert best == 5
    assert best_cmp.log_rmse == pytest.approx(0.0, abs=1e-12)
    assert len(table) == 5


def test_summary_is_string():
    a = np.full((4, 4), 12.0)
    b = np.full((4, 4), 10.0)
    s = validate.compare_fields(a, b).summary()
    assert "bias" in s and "classification" in s
    assert isinstance(s, str)


# --- Direction comparison (circular) ------------------------------------------

def test_identical_directions():
    a = np.array([[0.0, 90.0], [180.0, 270.0]])
    c = validate.compare_directions(a, a.copy())
    assert c.mean_abs_deg == pytest.approx(0.0)
    assert c.mean_cos == pytest.approx(1.0)
    assert c.within_5deg == 1.0


def test_direction_wraps_around_north():
    # 359 vs 1 deg is a 2-deg error, not 358.
    a = np.array([[359.0]])
    b = np.array([[1.0]])
    c = validate.compare_directions(a, b)
    assert c.mean_abs_deg == pytest.approx(2.0)


def test_direction_opposite_is_180():
    a = np.array([[0.0, 90.0]])
    b = np.array([[180.0, 270.0]])
    c = validate.compare_directions(a, b)
    assert c.mean_abs_deg == pytest.approx(180.0)
    assert c.mean_cos == pytest.approx(-1.0)
    assert c.within_30deg == 0.0


# --- perimeter / arrival-time comparison --------------------------------------

def _disk(n, cr, cc, radius):
    rr, cc_ = np.mgrid[0:n, 0:n]
    return (rr - cr) ** 2 + (cc_ - cc) ** 2 <= radius ** 2


def test_identical_perimeters_score_perfectly():
    a = _disk(40, 20, 20, 8)
    c = validate.compare_perimeters(a, a.copy(), cellsize_x=30.0)
    assert c.jaccard == pytest.approx(1.0)
    assert c.dice == pytest.approx(1.0)
    assert c.area_ratio == pytest.approx(1.0)
    assert c.mean_perimeter_distance == pytest.approx(0.0)
    assert c.hausdorff == pytest.approx(0.0)


def test_disjoint_perimeters_zero_overlap():
    a = _disk(60, 15, 15, 6)
    b = _disk(60, 45, 45, 6)
    c = validate.compare_perimeters(a, b, cellsize_x=30.0)
    assert c.jaccard == 0.0
    assert c.mean_perimeter_distance > 0.0


def test_overlapping_perimeters_partial():
    a = _disk(60, 30, 28, 10)
    b = _disk(60, 30, 32, 10)
    c = validate.compare_perimeters(a, b, cellsize_x=30.0)
    assert 0.0 < c.jaccard < 1.0
    assert c.dice > c.jaccard            # Dice >= Jaccard always
    assert c.hausdorff >= c.mean_perimeter_distance > 0.0


def test_compare_arrival_times_overlap_only():
    a = np.full((10, 10), np.inf)
    b = np.full((10, 10), np.inf)
    a[:5, :5] = 10.0
    b[:5, :5] = 10.0                      # same timing where both burn
    b[5:, 5:] = 99.0                      # only in B, beyond max_time anyway
    c = validate.compare_arrival_times(a, b, max_time=60.0)
    assert c.n == 25
    assert c.rmse == pytest.approx(0.0)


def test_direction_within_tolerances():
    a = np.array([10.0, 10.0, 10.0, 10.0])
    b = np.array([12.0, 18.0, 35.0, 200.0])   # errors 2, 8, 25, 170 deg
    c = validate.compare_directions(a, b)
    assert c.within_5deg == pytest.approx(0.25)
    assert c.within_10deg == pytest.approx(0.50)
    assert c.within_30deg == pytest.approx(0.75)
