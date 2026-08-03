# Montale 2017 — Module 5: Pyroconvection (real ERA5 sounding)

Profile: **ERA5 reanalysis** at Montale (43.96 N, 11.04 E), 16 Jul 2017 **15:00 UTC
(17:00 local, peak burn)**, 21 pressure levels + single-level surface.

## The environment
| | value |
|---|---|
| Surface T / dewpoint | 24.1 °C / 7.9 °C |
| Surface CAPE | **0 J/kg** |
| ERA5 boundary-layer height | 1943 m |
| Bulk-Richardson ABL (pyflam) | 1943 m |
| LCL (Bolton) | 2011 m |
| LCL / ABL ratio | 1.03 |
| ML dθ/dz | 6.63e-04 K/m |
| Cap γ-θ (ABL+200..1200) | 6.49e-03 K/m |
| RH at ABL top | 60 % |
| Mid-trop RH (min 700–650 hPa) | ~10–14 % (very dry) |
| Continuous Haines | 8.8 |
| Inverted-V (dry unstable) | False (sfc dewpoint depr 18.5 °C, mid RH 17 %) |

## Pyroconvection assessment
- **Classification:** `overshooting_pyrocu`
- Peak fireline intensity 45922 kW/m (99th pct) → mixed-layer flux 11814 W/m²
- **Fire-induced boundary layer (fireABL): 16391 m**; decoupling ratio 8.4×

This is the classic **dry-pyroconvection / plume-dominated** environment the literature
identifies: a **deep dry boundary layer** with **near-zero CAPE** and a **very dry
mid-troposphere** — pyroconvection here is a *vertical* (LCL/ABL/dryness-aloft)
problem, not a surface-CAPE one. The strong fire grows a fireABL well above the
ambient ABL (decoupling 8.4×). fireABL magnitude is a diagnostic (read the
ratio); C-Haines 8.8 independently flags the dry/unstable atmosphere.
