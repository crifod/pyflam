# Tuscany — Meteo · Pyroconvection · Fire-danger — GFS 2026-06-26 00Z run

*3-hourly (fxx 0–21, valid 26 Jun 00Z–21Z UTC; CEST = UTC+2). Live GFS 0.25° forcing through pyflam 0.1.0: surface state + GFS 850/700/600/500 hPa soundings → meteo_report, pyroconvection_potential, Continuous-Haines, PFT; Rothermel surface + Scott&Reinhardt/Cruz crown fire.*

## Regional pyroCb watch (peak heating ≈14:00 CEST / 12Z)

| Province | C-Haines | Inverted-V | LCL m | Favourable |
|----------|----------|-----------|-------|-----------|
| Grosseto | 6.4 | no | 1323 | no |
| Siena | 6.2 | no | 3243 | no |
| Arezzo | 5.7 | no | 2967 | no |
| Firenze | 4.8 | no | 3222 | no |
| Pistoia | 4.6 | no | 2541 | no |
| Prato | 4.6 | no | 2541 | no |
| Lucca | 4.5 | no | 2717 | no |
| Pisa | 4.5 | no | 2717 | no |
| Livorno | 4.5 | no | 1216 | no |
| Massa-Carrara | 3.9 | no | 1622 | no |

## Regional fire-danger (peak heating, crown model = Cruz 2005)

| Province | Fuel | 1h FM % | ROS m/min | Flame m | Intensity kW/m | Crown type | Torch mph | Crown mph |
|----------|------|---------|-----------|---------|----------------|-----------|-----------|-----------|
| Arezzo | SH5 | 5 | 42.9 | 4.1 | 11932 | active | 0 | 26 |
| Firenze | SH5 | 4 | 31.0 | 3.3 | 8138 | active | 0 | 26 |
| Grosseto | SH5 | 9 | 21.7 | 4.0 | 7725 | passive | 0 | 30 |
| Siena | SH5 | 4 | 16.1 | 2.9 | 3824 | passive | 0 | 26 |
| Prato | SH5 | 6 | 7.2 | 2.2 | 1799 | passive | 0 | 27 |
| Massa-Carrara | TU5 | 8 | 1.6 | 1.6 | 743 | surface | 13 | 21 |
| Pistoia | TU5 | 6 | 1.2 | 1.5 | 610 | surface | 11 | 20 |
| Lucca | GS2 | 5 | 2.6 | 0.9 | 224 | surface | 16 | 24 |
| Pisa | GS2 | 5 | 2.5 | 0.9 | 216 | surface | 16 | 24 |
| Livorno | GS2 | 10 | 2.0 | 0.5 | 63 | surface | 65 | 27 |

## Per-province detail

### Massa-Carrara (MS) — mountain · fuel TU5 · 44.037N 10.140E

| Valid UTC | T °C | RH % | Wind m/s (dir) | FM 1/10/100h % | CAPE | BLH m | Stab |
|-----------|------|------|----------------|----------------|------|-------|------|
| 00Z | 23 | 61 | 2.5 NE | 11/11/11 | 357 | - | neut |
| 03Z | 23 | 59 | 2.2 NE | 10/11/11 | 26 | - | neut |
| 06Z | 25 | 65 | 1.4 NE | 11/11/11 | 603 | - | unst |
| 09Z | 29 | 48 | 1.8 W | 9/10/11 | 852 | - | unst |
| 12Z | 30 | 46 | 2.3 W | 8/10/11 | 1174 | - | unst |
| 15Z | 30 | 48 | 2.5 NW | 9/9/11 | 1452 | - | unst |
| 18Z | 27 | 67 | 2.6 N | 12/10/11 | 2376 | - | unst |
| 21Z | 24 | 63 | 2.2 NE | 11/10/11 | 575 | - | unst |

**Vertical / pyroCb @12Z** — LCL **1622 m**, mixed-layer depth -, C-Haines **3.9**, inverted-V **no** (mid-600 hPa RH 41%), cap N² 1.5e-04 s⁻², plume ×2.14. PFT ≈ **16.4 GW**; a 2 GW fireline lifts ≈ 954 m. plume-dominated favourable: **no**.

**Fire-danger @12Z** — 1-h FM 8%, midflame wind 2.1 mph; surface ROS 1.6 m/min, flame 1.6 m, surface intensity 743 kW/m. Crown: **surface**, total ROS 1.6 m/min, intensity 743 kW/m. Torching index 13 mph / crowning index 21 mph.

_pyroCb evolution_: LOW — capped/moist column; pyroCb not expected.  _Drivers 00→21Z_: RH +2% (rising), T +1°C, wind -0.3 m/s, 1-h FM +0.3%.

### Lucca (LU) — coast · fuel GS2 · 43.843N 10.508E

| Valid UTC | T °C | RH % | Wind m/s (dir) | FM 1/10/100h % | CAPE | BLH m | Stab |
|-----------|------|------|----------------|----------------|------|-------|------|
| 00Z | 24 | 64 | 1.3 E | 11/11/11 | 163 | - | neut |
| 03Z | 22 | 64 | 1.4 E | 11/11/11 | 0 | - | neut |
| 06Z | 26 | 54 | 1.7 E | 10/11/11 | 0 | - | neut |
| 09Z | 33 | 36 | 1.1 N | 7/10/11 | 144 | - | neut |
| 12Z | 36 | 27 | 2.5 NW | 5/9/11 | 362 | - | neut |
| 15Z | 34 | 36 | 4.4 NW | 6/8/11 | 992 | - | unst |
| 18Z | 28 | 61 | 2.6 W | 10/9/11 | 1461 | - | unst |
| 21Z | 25 | 74 | 0.7 NE | 13/10/11 | 1032 | - | unst |

**Vertical / pyroCb @12Z** — LCL **2717 m**, mixed-layer depth -, C-Haines **4.5**, inverted-V **no** (mid-600 hPa RH 48%), cap N² 1.5e-04 s⁻², plume ×1.24. PFT ≈ **77.1 GW**; a 2 GW fireline lifts ≈ 954 m. plume-dominated favourable: **no**.

**Fire-danger @12Z** — 1-h FM 5%, midflame wind 2.2 mph; surface ROS 2.6 m/min, flame 0.9 m, surface intensity 224 kW/m. Crown: **surface**, total ROS 2.6 m/min, intensity 224 kW/m. Torching index 16 mph / crowning index 24 mph.

_pyroCb evolution_: MODERATE — buoyant plume, lofting/spotting, deep pyroCb not favoured.  _Drivers 00→21Z_: RH +9% (rising), T +1°C, wind -0.6 m/s, 1-h FM +2.1%.

### Pistoia (PT) — mountain · fuel TU5 · 43.933N 10.918E

| Valid UTC | T °C | RH % | Wind m/s (dir) | FM 1/10/100h % | CAPE | BLH m | Stab |
|-----------|------|------|----------------|----------------|------|-------|------|
| 00Z | 22 | 51 | 2.7 N | 10/10/10 | 0 | - | neut |
| 03Z | 22 | 53 | 2.9 N | 10/10/10 | 0 | - | neut |
| 06Z | 26 | 48 | 1.5 N | 9/9/10 | 4 | - | neut |
| 09Z | 31 | 32 | 1.4 E | 6/9/9 | 57 | - | neut |
| 12Z | 33 | 29 | 1.5 E | 6/8/9 | 472 | - | neut |
| 15Z | 33 | 32 | 2.3 N | 6/7/9 | 748 | - | unst |
| 18Z | 27 | 48 | 2.7 NE | 9/8/9 | 369 | - | neut |
| 21Z | 24 | 46 | 3.1 N | 8/8/9 | 0 | - | neut |

**Vertical / pyroCb @12Z** — LCL **2541 m**, mixed-layer depth -, C-Haines **4.6**, inverted-V **no** (mid-600 hPa RH 22%), cap N² 1.6e-04 s⁻², plume ×1.31. PFT ≈ **42.3 GW**; a 2 GW fireline lifts ≈ 1089 m. plume-dominated favourable: **no**.

**Fire-danger @12Z** — 1-h FM 6%, midflame wind 1.3 mph; surface ROS 1.2 m/min, flame 1.5 m, surface intensity 610 kW/m. Crown: **surface**, total ROS 1.2 m/min, intensity 610 kW/m. Torching index 11 mph / crowning index 20 mph.

_pyroCb evolution_: MODERATE — buoyant plume, lofting/spotting, deep pyroCb not favoured.  _Drivers 00→21Z_: RH -6% (falling), T +1°C, wind +0.4 m/s, 1-h FM -1.2%.

### Prato (PO) — basin · fuel SH5 · 43.881N 11.097E

| Valid UTC | T °C | RH % | Wind m/s (dir) | FM 1/10/100h % | CAPE | BLH m | Stab |
|-----------|------|------|----------------|----------------|------|-------|------|
| 00Z | 22 | 51 | 2.7 N | 10/10/10 | 0 | - | neut |
| 03Z | 22 | 53 | 2.9 N | 10/10/10 | 0 | - | neut |
| 06Z | 26 | 48 | 1.5 N | 9/9/10 | 4 | - | neut |
| 09Z | 31 | 32 | 1.4 E | 6/9/9 | 57 | - | neut |
| 12Z | 33 | 29 | 1.5 E | 6/8/9 | 472 | - | neut |
| 15Z | 33 | 32 | 2.3 N | 6/7/9 | 748 | - | unst |
| 18Z | 27 | 48 | 2.7 NE | 9/8/9 | 369 | - | neut |
| 21Z | 24 | 46 | 3.1 N | 8/8/9 | 0 | - | neut |

**Vertical / pyroCb @12Z** — LCL **2541 m**, mixed-layer depth -, C-Haines **4.6**, inverted-V **no** (mid-600 hPa RH 22%), cap N² 1.6e-04 s⁻², plume ×1.31. PFT ≈ **42.3 GW**; a 2 GW fireline lifts ≈ 1089 m. plume-dominated favourable: **no**.

**Fire-danger @12Z** — 1-h FM 6%, midflame wind 1.3 mph; surface ROS 4.9 m/min, flame 2.2 m, surface intensity 1471 kW/m. Crown: **passive**, total ROS 7.2 m/min, intensity 1799 kW/m. Torching index 0 mph / crowning index 27 mph.

_pyroCb evolution_: MODERATE — buoyant plume, lofting/spotting, deep pyroCb not favoured.  _Drivers 00→21Z_: RH -6% (falling), T +1°C, wind +0.4 m/s, 1-h FM -1.2%.

### Firenze (FI) — basin · fuel SH5 · 43.769N 11.256E

| Valid UTC | T °C | RH % | Wind m/s (dir) | FM 1/10/100h % | CAPE | BLH m | Stab |
|-----------|------|------|----------------|----------------|------|-------|------|
| 00Z | 25 | 48 | 1.3 E | 9/9/9 | 0 | - | neut |
| 03Z | 23 | 51 | 1.6 E | 10/9/9 | 0 | - | neut |
| 06Z | 27 | 36 | 0.8 NE | 7/8/9 | 0 | - | neut |
| 09Z | 35 | 24 | 1.1 NE | 5/7/9 | 0 | - | neut |
| 12Z | 37 | 21 | 3.1 NE | 4/7/8 | 257 | - | neut |
| 15Z | 38 | 21 | 4.7 NE | 4/6/8 | 484 | - | neut |
| 18Z | 34 | 27 | 3.0 E | 5/6/8 | 63 | - | neut |
| 21Z | 28 | 40 | 1.1 NE | 7/6/8 | 0 | - | neut |

**Vertical / pyroCb @12Z** — LCL **3222 m**, mixed-layer depth -, C-Haines **4.8**, inverted-V **no** (mid-600 hPa RH 18%), cap N² 1.6e-04 s⁻², plume ×1.17. PFT ≈ **184.3 GW**; a 2 GW fireline lifts ≈ 846 m. plume-dominated favourable: **no**.

**Fire-danger @12Z** — 1-h FM 4%, midflame wind 2.8 mph; surface ROS 11.1 m/min, flame 3.3 m, surface intensity 3488 kW/m. Crown: **active**, total ROS 31.0 m/min, intensity 8138 kW/m. Torching index 0 mph / crowning index 26 mph.

_pyroCb evolution_: MODERATE — buoyant plume, lofting/spotting, deep pyroCb not favoured.  _Drivers 00→21Z_: RH -8% (falling), T +3°C, wind -0.3 m/s, 1-h FM -1.5%.

### Pisa (PI) — coast · fuel GS2 · 43.716N 10.401E

| Valid UTC | T °C | RH % | Wind m/s (dir) | FM 1/10/100h % | CAPE | BLH m | Stab |
|-----------|------|------|----------------|----------------|------|-------|------|
| 00Z | 24 | 64 | 1.3 E | 11/11/11 | 163 | - | neut |
| 03Z | 22 | 64 | 1.4 E | 11/11/11 | 0 | - | neut |
| 06Z | 26 | 54 | 1.7 E | 10/11/11 | 0 | - | neut |
| 09Z | 33 | 36 | 1.1 N | 7/10/11 | 144 | - | neut |
| 12Z | 36 | 27 | 2.5 NW | 5/9/11 | 362 | - | neut |
| 15Z | 34 | 36 | 4.4 NW | 6/8/11 | 992 | - | unst |
| 18Z | 28 | 61 | 2.6 W | 10/9/11 | 1461 | - | unst |
| 21Z | 25 | 74 | 0.7 NE | 13/10/11 | 1032 | - | unst |

**Vertical / pyroCb @12Z** — LCL **2717 m**, mixed-layer depth -, C-Haines **4.5**, inverted-V **no** (mid-600 hPa RH 48%), cap N² 1.5e-04 s⁻², plume ×1.24. PFT ≈ **77.1 GW**; a 2 GW fireline lifts ≈ 954 m. plume-dominated favourable: **no**.

**Fire-danger @12Z** — 1-h FM 5%, midflame wind 2.2 mph; surface ROS 2.5 m/min, flame 0.9 m, surface intensity 216 kW/m. Crown: **surface**, total ROS 2.5 m/min, intensity 216 kW/m. Torching index 16 mph / crowning index 24 mph.

_pyroCb evolution_: MODERATE — buoyant plume, lofting/spotting, deep pyroCb not favoured.  _Drivers 00→21Z_: RH +9% (rising), T +1°C, wind -0.6 m/s, 1-h FM +2.1%.

### Livorno (LI) — coast · fuel GS2 · 43.548N 10.311E

| Valid UTC | T °C | RH % | Wind m/s (dir) | FM 1/10/100h % | CAPE | BLH m | Stab |
|-----------|------|------|----------------|----------------|------|-------|------|
| 00Z | 26 | 78 | 1.4 N | 15/15/15 | 1965 | - | unst |
| 03Z | 26 | 69 | 1.4 NE | 12/14/15 | 428 | - | neut |
| 06Z | 26 | 66 | 1.9 SE | 12/14/15 | 180 | - | neut |
| 09Z | 27 | 60 | 1.8 W | 10/13/15 | 243 | - | neut |
| 12Z | 28 | 56 | 4.6 NW | 10/12/15 | 176 | - | neut |
| 15Z | 28 | 69 | 6.4 NW | 12/12/14 | 1907 | - | unst |
| 18Z | 27 | 74 | 5.1 NW | 14/13/14 | 2151 | - | unst |
| 21Z | 26 | 81 | 2.7 NW | 16/13/15 | 2049 | - | unst |

**Vertical / pyroCb @12Z** — LCL **1216 m**, mixed-layer depth -, C-Haines **4.5**, inverted-V **no** (mid-600 hPa RH 53%), cap N² 1.4e-04 s⁻², plume ×1.12. PFT ≈ **12.4 GW**; a 2 GW fireline lifts ≈ 784 m. plume-dominated favourable: **no**.

**Fire-danger @12Z** — 1-h FM 10%, midflame wind 4.1 mph; surface ROS 2.0 m/min, flame 0.5 m, surface intensity 63 kW/m. Crown: **surface**, total ROS 2.0 m/min, intensity 63 kW/m. Torching index 65 mph / crowning index 27 mph.

_pyroCb evolution_: LOW — capped/moist column; pyroCb not expected.  _Drivers 00→21Z_: RH +3% (rising), T +0°C, wind +1.2 m/s, 1-h FM +1.0%.

### Arezzo (AR) — basin · fuel SH5 · 43.463N 11.879E

| Valid UTC | T °C | RH % | Wind m/s (dir) | FM 1/10/100h % | CAPE | BLH m | Stab |
|-----------|------|------|----------------|----------------|------|-------|------|
| 00Z | 19 | 55 | 2.1 NE | 10/10/10 | 0 | - | neut |
| 03Z | 18 | 55 | 1.8 NE | 10/10/10 | 0 | - | neut |
| 06Z | 25 | 45 | 1.4 NE | 8/10/10 | 0 | - | neut |
| 09Z | 32 | 25 | 4.2 NE | 5/8/10 | 0 | - | neut |
| 12Z | 34 | 24 | 4.9 NE | 5/7/10 | 39 | - | neut |
| 15Z | 33 | 26 | 5.5 NE | 5/7/10 | 136 | - | neut |
| 18Z | 28 | 38 | 3.4 NE | 7/7/9 | 0 | - | neut |
| 21Z | 21 | 51 | 2.5 NE | 10/8/9 | 0 | - | neut |

**Vertical / pyroCb @12Z** — LCL **2967 m**, mixed-layer depth -, C-Haines **5.7**, inverted-V **no** (mid-600 hPa RH 16%), cap N² 1.7e-04 s⁻², plume ×1.09. PFT ≈ **226.8 GW**; a 2 GW fireline lifts ≈ 727 m. plume-dominated favourable: **no**.

**Fire-danger @12Z** — 1-h FM 5%, midflame wind 4.4 mph; surface ROS 17.8 m/min, flame 4.1 m, surface intensity 5500 kW/m. Crown: **active**, total ROS 42.9 m/min, intensity 11932 kW/m. Torching index 0 mph / crowning index 26 mph.

_pyroCb evolution_: MODERATE — buoyant plume, lofting/spotting, deep pyroCb not favoured.  _Drivers 00→21Z_: RH -4% (falling), T +1°C, wind +0.4 m/s, 1-h FM -0.5%.

### Siena (SI) — basin · fuel SH5 · 43.318N 11.331E

| Valid UTC | T °C | RH % | Wind m/s (dir) | FM 1/10/100h % | CAPE | BLH m | Stab |
|-----------|------|------|----------------|----------------|------|-------|------|
| 00Z | 22 | 62 | 1.0 NE | 11/11/11 | 0 | - | neut |
| 03Z | 20 | 67 | 1.0 E | 12/11/11 | 0 | - | neut |
| 06Z | 24 | 56 | 0.8 E | 10/11/11 | 0 | - | neut |
| 09Z | 31 | 34 | 1.9 NE | 7/10/11 | 0 | - | neut |
| 12Z | 36 | 21 | 2.4 NE | 4/8/11 | 5 | - | neut |
| 15Z | 36 | 21 | 3.5 NE | 4/7/10 | 138 | - | neut |
| 18Z | 28 | 50 | 1.7 NW | 9/8/10 | 558 | - | unst |
| 21Z | 24 | 52 | 0.4 NE | 10/8/10 | 0 | - | neut |

**Vertical / pyroCb @12Z** — LCL **3243 m**, mixed-layer depth -, C-Haines **6.2**, inverted-V **no** (mid-600 hPa RH 14%), cap N² 1.6e-04 s⁻², plume ×1.12. PFT ≈ **141.4 GW**; a 2 GW fireline lifts ≈ 930 m. plume-dominated favourable: **no**.

**Fire-danger @12Z** — 1-h FM 4%, midflame wind 2.2 mph; surface ROS 8.6 m/min, flame 2.9 m, surface intensity 2673 kW/m. Crown: **passive**, total ROS 16.1 m/min, intensity 3824 kW/m. Torching index 0 mph / crowning index 26 mph.

_pyroCb evolution_: MODERATE — buoyant plume, lofting/spotting, deep pyroCb not favoured.  _Drivers 00→21Z_: RH -10% (falling), T +2°C, wind -0.6 m/s, 1-h FM -1.4%.

### Grosseto (GR) — basin · fuel SH5 · 42.760N 11.113E

| Valid UTC | T °C | RH % | Wind m/s (dir) | FM 1/10/100h % | CAPE | BLH m | Stab |
|-----------|------|------|----------------|----------------|------|-------|------|
| 00Z | 24 | 73 | 1.7 E | 14/14/14 | 401 | - | neut |
| 03Z | 23 | 75 | 2.0 E | 14/14/14 | 357 | - | neut |
| 06Z | 25 | 63 | 2.3 NE | 11/13/14 | 2 | - | neut |
| 09Z | 29 | 47 | 0.8 S | 9/12/13 | 0 | - | neut |
| 12Z | 30 | 53 | 6.0 W | 9/11/13 | 505 | - | unst |
| 15Z | 29 | 58 | 5.6 W | 10/11/13 | 1027 | - | unst |
| 18Z | 28 | 64 | 3.0 W | 11/11/13 | 1037 | - | unst |
| 21Z | 26 | 66 | 2.0 N | 12/11/13 | 364 | - | neut |

**Vertical / pyroCb @12Z** — LCL **1323 m**, mixed-layer depth -, C-Haines **6.4**, inverted-V **no** (mid-600 hPa RH 14%), cap N² 1.5e-04 s⁻², plume ×1.82. PFT ≈ **22.0 GW**; a 2 GW fireline lifts ≈ 706 m. plume-dominated favourable: **no**.

**Fire-danger @12Z** — 1-h FM 9%, midflame wind 5.3 mph; surface ROS 19.0 m/min, flame 4.0 m, surface intensity 5226 kW/m. Crown: **passive**, total ROS 21.7 m/min, intensity 7725 kW/m. Torching index 0 mph / crowning index 30 mph.

_pyroCb evolution_: MODERATE — buoyant plume, lofting/spotting, deep pyroCb not favoured.  _Drivers 00→21Z_: RH -7% (falling), T +2°C, wind +0.3 m/s, 1-h FM -2.0%.

---
**Source:** GFS 0.25° operational forecast, 2026-06-26 00Z run (AWS Open-Data via Herbie), fxx 0–21. Surface: 10 m wind, 2 m T/RH, surface CAPE, PBL height. Sounding: GFS T/RH at 850/700/600/500 hPa + the GFS surface point. Fuel models: Scott&Burgan SB40 (SH5 interior shrub/maquis, GS2 coastal grass-shrub, TU5 Apennine timber-shrub); canopy structure and slope are representative per regime (swap in an .lcp for cell-resolved structure). Live-fire moisture 70%/90%.
