"""Tests for the operative sector / driving-force analysis."""

from __future__ import annotations

import numpy as np
import pytest

import pyflam
from pyflam import operative as op
from pyflam.units import mph_to_ft_per_min

SC = dict(m_1h=0.06, m_10h=0.07, m_100h=0.08, m_live_herb=0.6, m_live_woody=0.9)


def _landscape(n=81, fuel=104, slope=None, aspect=180.0):
    return pyflam.Landscape(
        fuel_model=(fuel if isinstance(fuel, np.ndarray)
                    else np.full((n, n), fuel, dtype=int)),
        slope=(slope if slope is not None else np.zeros((n, n))),
        aspect=np.full((n, n), aspect), cellsize_x=30.0, cellsize_y=30.0,
        west=0.0, north=n * 30.0, slope_units="degrees")


# --- driver fields ------------------------------------------------------------

def test_driver_fields_keys_and_signs():
    ls = _landscape(n=10)
    d = op.driver_fields(ls, wind_midflame=mph_to_ft_per_min(8),
                         wind_direction=270.0, **SC)
    assert set(d) >= {"phi_w", "phi_s", "r0", "i0", "ros_max", "heading"}
    assert (d["phi_w"] > 0).all()                # wind present
    assert (d["phi_s"] == 0).all()               # flat
    assert (d["r0"] > 0).all() and (d["i0"] > 0).all()


# --- perimeter extraction -----------------------------------------------------

def test_perimeter_is_the_edge_ring():
    n = 21
    arrival = np.full((n, n), np.inf)
    rr, cc = np.mgrid[0:n, 0:n]
    arrival[(rr - 10) ** 2 + (cc - 10) ** 2 <= 25] = 1.0   # a burned disk
    edge = op.perimeter_cells(arrival, 10.0)
    assert edge[10, 10] == False                 # interior not on the edge
    assert edge.sum() > 0
    # every edge cell is burned and touches an unburned neighbour
    burned = np.isfinite(arrival)
    er, ec = np.where(edge)
    assert burned[er, ec].all()


# --- sector analysis ----------------------------------------------------------

def _run(ls, wind_mph, wind_dir, max_time=80):
    field = pyflam.spread_field(ls, wind_midflame=mph_to_ft_per_min(wind_mph),
                                wind_direction=wind_dir, **SC)
    arr = pyflam.minimum_travel_time(field, [(ls.shape[0] // 2, ls.shape[1] // 2)],
                                     max_time=max_time)
    rep = pyflam.analyze_perimeter(
        ls, arr, max_time, wind_midflame=mph_to_ft_per_min(wind_mph),
        wind_direction=wind_dir, **SC)
    return rep


def test_four_sectors_and_heading():
    rep = _run(_landscape(), wind_mph=8, wind_dir=270.0)   # wind from west
    assert {s.name for s in rep.sectors} == set(op.SECTORS)
    assert rep.heading == pytest.approx(90.0, abs=8)        # heads east


def test_flat_uniform_is_wind_driven_everywhere():
    rep = _run(_landscape(), wind_mph=8, wind_dir=270.0)
    for s in rep.sectors:
        assert s.dominant == "wind"
        assert s.slope.magnitude == pytest.approx(0.0, abs=1e-6)


def test_slope_drives_the_head_upslope():
    n = 81
    slope = np.zeros((n, n))
    for r in range(n):
        slope[r, :] = max(0.0, 40 - r)            # steeper to the north
    ls = _landscape(n=n, slope=slope, aspect=180.0)   # faces S -> upslope = N
    rep = _run(ls, wind_mph=3, wind_dir=180.0, max_time=90)  # light S wind
    head = next(s for s in rep.sectors if s.name == "head")
    assert head.slope.magnitude > 0.0
    assert head.slope.bearing == pytest.approx(0.0, abs=1.0)   # pushes north
    assert head.dominant == "slope"


def test_fuel_vector_points_to_more_dangerous_fuel():
    n = 81
    fuel = np.full((n, n), 104, dtype=int)        # grass west
    fuel[:, 40:] = 147                            # high-load shrub SH7 east
    ls = _landscape(n=n, fuel=fuel)
    rep = _run(ls, wind_mph=4, wind_dir=225.0, max_time=90)
    # somewhere on the perimeter the fuel push has an eastward (toward shrub) sign
    east_comp = [op.ForceVector("fuel", s.fuel.magnitude, s.fuel.bearing).components[0]
                 for s in rep.sectors if s.fuel.magnitude > 0]
    assert any(e > 0 for e in east_comp)


def test_resultant_is_vector_sum_of_components():
    rep = _run(_landscape(), wind_mph=8, wind_dir=270.0)
    for s in rep.sectors:
        we, wn = s.wind.components
        se, sn = s.slope.components
        fe, fn = s.fuel.components
        re, rn = s.resultant.components
        assert re == pytest.approx(we + se + fe, abs=1e-6)
        assert rn == pytest.approx(wn + sn + fn, abs=1e-6)


def test_subsectors_split_each_sector():
    rep = pyflam.analyze_perimeter(
        _landscape(), *_arr_args(), subsectors=2,
        wind_midflame=mph_to_ft_per_min(8), wind_direction=270.0, **SC)
    names = [s.name for s in rep.sectors]
    assert all("/2" in n for n in names)
    assert len(rep.sectors) <= 8


def _arr_args():
    ls = _landscape()
    field = pyflam.spread_field(ls, wind_midflame=mph_to_ft_per_min(8),
                                wind_direction=270.0, **SC)
    arr = pyflam.minimum_travel_time(field, [(40, 40)], max_time=80)
    return arr, 80


# --- GeoJSON export -----------------------------------------------------------

def _geo_landscape_run():
    n = 61
    fuel = np.full((n, n), 104, dtype=int)
    fuel[:, 40:] = 147
    ls = pyflam.Landscape(
        fuel_model=fuel, slope=np.zeros((n, n)), aspect=np.full((n, n), 180.0),
        cellsize_x=100.0, cellsize_y=100.0, west=4_300_000.0, north=2_370_000.0,
        slope_units="degrees", crs="EPSG:3035")
    field = pyflam.spread_field(ls, wind_midflame=mph_to_ft_per_min(6),
                                wind_direction=225.0, **SC)
    arr = pyflam.minimum_travel_time(field, [(30, 30)], max_time=120)
    rep = pyflam.analyze_perimeter(ls, arr, 120, wind_midflame=mph_to_ft_per_min(6),
                                   wind_direction=225.0, **SC)
    return ls, rep


def test_geojson_structure():
    ls, rep = _geo_landscape_run()
    gj = rep.to_geojson(ls, arrow_length=2000.0)
    assert gj["type"] == "FeatureCollection"
    kinds = [f["properties"]["kind"] for f in gj["features"]]
    assert kinds.count("sector") == 4
    assert kinds.count("arrow") == 16          # 4 sectors x 4 forces
    assert kinds.count("perimeter") == 1
    assert gj["crs"]["properties"]["name"] == "EPSG:3035"


def test_geojson_sector_and_arrow_properties():
    ls, rep = _geo_landscape_run()
    gj = rep.to_geojson(ls)
    sector = next(f for f in gj["features"]
                  if f["properties"]["kind"] == "sector")
    p = sector["properties"]
    assert {"wind_mag", "slope_mag", "fuel_mag", "resultant_mag",
            "dominant", "mean_ros_ft_min"} <= set(p)
    arrow = next(f for f in gj["features"]
                 if f["properties"].get("force") == "resultant")
    assert arrow["geometry"]["type"] == "LineString"
    assert len(arrow["geometry"]["coordinates"]) == 2
    # coordinates fall within the landscape's world extent
    (x0, y0), (x1, y1) = arrow["geometry"]["coordinates"]
    assert 4_300_000 - 1e4 < x0 < 4_300_000 + 61 * 100 + 1e4


def test_geojson_reproject_to_wgs84():
    pytest.importorskip("pyproj")
    ls, rep = _geo_landscape_run()
    gj = rep.to_geojson(ls, to_wgs84=True)
    assert gj["crs"]["properties"]["name"] == "EPSG:4326"
    pt = next(f for f in gj["features"]
              if f["properties"]["kind"] == "sector")
    lon, lat = pt["geometry"]["coordinates"]
    assert -180 <= lon <= 180 and 30 < lat < 60      # plausible Europe lat/lon


def test_perimeter_polygon_ring_is_ordered_and_closed():
    pytest.importorskip("skimage")
    ls, rep = _geo_landscape_run()
    gj = rep.to_geojson(ls, perimeter_geometry="polygon")
    perim = next(f for f in gj["features"]
                 if f["properties"]["kind"] == "perimeter")
    assert perim["geometry"]["type"] in ("Polygon", "MultiPolygon")
    ring = (perim["geometry"]["coordinates"][0] if perim["geometry"]["type"]
            == "Polygon" else perim["geometry"]["coordinates"][0][0])
    assert len(ring) >= 4
    assert ring[0] == ring[-1]            # closed ring


def test_perimeter_points_fallback():
    ls, rep = _geo_landscape_run()
    gj = rep.to_geojson(ls, perimeter_geometry="points")
    perim = next(f for f in gj["features"]
                 if f["properties"]["kind"] == "perimeter")
    assert perim["geometry"]["type"] == "MultiPoint"


def test_perimeter_rings_helper():
    pytest.importorskip("skimage")
    ls, rep = _geo_landscape_run()
    rings = op.perimeter_rings(rep.burned, ls)
    assert rings and all(r[0] == r[-1] for r in rings)   # each closed


def test_write_geojson_file(tmp_path):
    import json
    ls, rep = _geo_landscape_run()
    path = tmp_path / "operative.geojson"
    rep.write_geojson(ls, str(path))
    gj = json.loads(path.read_text())
    assert gj["type"] == "FeatureCollection" and gj["features"]


def test_no_perimeter_raises():
    ls = _landscape(n=11)
    arrival = np.full((11, 11), np.inf)           # nothing burned
    with pytest.raises(ValueError):
        pyflam.analyze_perimeter(ls, arrival, 30.0,
                                 wind_midflame=mph_to_ft_per_min(5),
                                 wind_direction=270.0, **SC)
