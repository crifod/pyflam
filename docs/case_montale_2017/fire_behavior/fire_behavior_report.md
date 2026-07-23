# Montale 2017 — Module 2: Surface fire behavior (Rothermel) + arrival

Ignition cell (300, 300) (43.96095 N, 11.04026 E). Anisotropic-Eikonal front,
alpha_samples=16, per-cell midflame wind (WAF from canopy).

| metric | median | 90th pct | max |
|---|---|---|---|
| Rate of spread (m/min) | 17.8 | 52.3 | 121.9 |
| Fireline intensity (kW/m) | 3483 | 14117 | 103924 |
| Flame length (m) | 3.3 | 6.3 | 15.7 |

- Burnable area reached in the sim window: **3484 ha** (vs ~300 ha reported).
- Max arrival time in domain: inf min.

Rasters: ros_m_per_min, fireline_intensity_kW_m, flame_length_m, arrival_time_min (GeoTIFF, EPSG:3035).
