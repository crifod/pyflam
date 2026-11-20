# GRAF/WUR 2021 reference wildfire campaign — complete corpus

The calibration data behind the pyflam dry-pyrocloud fireABL work comes from the
Bombers de la Generalitat de Catalunya (GRAF) / Wageningen University & Research
2021 field campaign (Castellnou et al. 2022, *JGR: Atmospheres*; datasets
`zenodo_6433389`, `zenodo_10823337`, `zenodo_12946249`). It is **7 fires,
144 fire-hours**, with three tiers of data completeness.

## Data-availability matrix

| # | Fire | Region | Dates 2021 | Behaviour (144 h) | Spread + fireABL | ERA5 sounding | Fire sondes | LES |
|---|---|---|---|---:|---:|:--:|:--:|:--:|
| 1 | **Santa Coloma de Queralt** (SCQ) | Conca de Barberà, Tarragona (Cat) | 24–25 Jul | 44 h | 19 h | ✅ | 3 | ✅ (Eghdami 2023) |
| 2 | **Martorell** | Baix Llobregat, Barcelona (Cat) | 11–13 Jul | 18 h | 7 h | ✅ | 3 | — |
| 3 | **Torroella de Montgrí** | Baix Empordà, Girona (Cat) | 20–22 Jul | 18 h | 5 h | ✅ | 1 | — |
| 4 | **Pobla de Massaluca** | Terra Alta, Tarragona (Cat) | 12 Aug | 8 h | 5 h | ✅ | 1 | — |
| 5 | **Sierra Bermeja** (≡ Jubrique) | Málaga (**Andalusia**) | 8–13 Sep | 23 h | — | ✅ (as Jubrique) | 3 (SB7/SB8) | — |
| 6 | **Alfarràs** | Segrià, Lleida (Cat) | 27–29 Jun | 11 h | — | — | — | — |
| 7 | **Llançà** | Alt Empordà, Girona (Cat) | 17–18 Jul | 22 h | — | — | — | — |

- **Tier 1 — full (SCQ):** behaviour + 733-level profile + observed fireABL + 3 sondes + coupled WRF-Fire LES. The anchor case.
- **Tier 2 — profile-usable (Martorell, Torroella, Pobla, Sierra Bermeja/Jubrique):** ERA5 sounding + sondes; the first four also have observed fireABL rows, forming the **36-hour multi-fire calibration set** (Sierra Bermeja has ERA5+sondes but no fireABL rows).
- **Tier 3 — classification-only (Alfarràs, Llançà):** Ctype/ABL/LCL/ROS in the 144-hour table; no sounding or fireABL.

## Per-fire detail

### 1. Santa Coloma de Queralt (SCQ) — the anchor
41.521 °N, 1.374 °E · elev 535–853 m · 448 ha · 24 Jul 14:00 Z → 25 Jul 19:00 Z.
Prototypes (144 h): 16 surface, 15 convective, 8 pyroCu, 4 oPyroCu, **1 pyroCb**.
Peak intensity **32 859 kW/m** (24 Jul 18:00 Z, deepest hour; fireABL_sonde 3346 m).
Evening sonde (25 Jul 20:19 Z, 7120 m) shows **RH 90 % at top** — the moist-aloft,
pyroCb-favourable profile. Only fire with a coupled-LES cross-check (cloud top ~8 km,
pyroCu onset ~18:28 Z). See `dry_pyroconvection_fireabl.md`.

### 2. Martorell — moderate, one pyroCu burst
41.465 °N, 1.939 °E · elev 166–272 m · 21 ha · 11–13 Jul. Ctype: 13 surface, 3
convective, 2 PyroCu. ABL 270–3850 m. Model fireABL tracks the sonde tightly
(e.g. 3066 vs 3062 m). Three sondes (13 Jul 18:05/18:18/19:02, tops 3720/3460/2280 m).

### 3. Torroella de Montgrí — the decisive decoupling case
42.061 °N, 3.147 °E (coastal) · elev 57–115 m · 28 ha · 22 Jul. Ctype: all
surface/convective, yet high intensity (**7459–8534 kW/m**) over a **very shallow
ABL (110–315 m)** → strong decoupling (fireABL_sonde 1133 m ≈ 3.6× ABL). This is the
fire where the fire-front flux under-predicts (780 m) and the **mixed-layer flux
nails it (+27 m)** — the single most important evidence for the flux reformulation.

### 4. Pobla de Massaluca — deep, hot, convective
41.224 °N, 0.347 °E · elev 109–166 m · 22 ha · 12 Aug. Ctype: 4 convective,
4 surface. Very hot (**avg θ 42 °C**), **deep ABL ~4.85 km**; both flux forms agree
here (model fireABL 3548 vs sonde 3474 m). Anchors the deep-ABL end of the range.

### 5. Sierra Bermeja (≡ Jubrique) — the out-of-region case
~36.5 °N, −5.1 °E, Málaga, Andalusia · 8–13 Sep (multi-week fire). The behaviour
table (23 h, 9–10 Sep) names it **Sierra Bermeja** with sondes **SB7/SB8**; the ERA5
(`Jubrique_20210905.nc`) and 3 sonde files (13 Sep, launched 36.45–36.53 °N /
−5.19 to −5.27 °E) use the launch-village name **Jubrique** — **verified same fire**.
Ctype: 14 convective, 9 surface; ABL 730–1950 m over the 9–10 Sep window (a clean
daytime CBL-growth cycle). The only non-Catalan, non-summer-peak reference.

### 6. Alfarràs — classification-only
Segrià, Lleida · 27–29 Jun · 11 h. Ctype: 8 surface, 2 convective, **1 PyroCu**
(27 Jun 19:00 Z, a ROS 7.1 km/h burst; ABL 1540 m, LCL 2235 m → overshooting). No
sounding/fireABL; profile pipeline cannot be run.

### 7. Llançà — classification-only
Alt Empordà, Girona (coast) · 17–18 Jul · 22 h. Ctype: 12 surface, 9 convective,
**1 oPyroCu** (17 Jul 14:00 Z; ABL 1920 m, LCL 2276 m). No sounding/fireABL.

## Role in the pyflam calibration
- **Port fidelity & LES referee:** SCQ (Tier 1) — machine-precision port check; peak
  fireABL ≈ 8 km LES cloud top; pyroCu onset timing.
- **Multi-fire magnitude confirmation (36 h):** SCQ + Martorell + Torroella + Pobla —
  mixed-layer flux bias −15 m, MAE 400 m, r 0.91; the decisive shallow-ABL case is
  Torroella.
- **Population-scale label separation (Stage A):** the full 90-fire / 795-hour
  `Fire_spread_events` set (a superset that includes these campaign fires) — decoupling
  ratio adds AUC 0.556 → 0.725 over LCL/ABL.

Sources: Castellnou, M. et al. (2022) *JGR: Atmospheres* 127(22), e2022JD036920;
Eghdami, M. et al. (2023) *JAMES* 15(5), e2022MS003288. Datasets CC BY 4.0.
