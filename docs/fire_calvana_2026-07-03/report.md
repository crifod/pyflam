# Calvana wildfire — pyflam propagation run

**Ignition** 43.93561, 11.09635 · **start** 2026-07-03 14:00 CEST (2026-07-03T12:00:00+00:00) · **duration** 120 min · dt 30 min
**Physics** ICON-2I + Cruz2005 crown + OpenFOAM RANS plume + pyroconvection + spotting
**AOI** 180×180 @100 m · fuel model at ignition SB40 #165

**Total burned area:** 1211.8 ha · **crown-involved:** 639.9 ha

## Time-of-arrival / behaviour every 30 min

| t (min) | local | area ha | ROS head p95 (m/min) | FLI max (kW/m) | flame (m) | wind m/s | RH % | 1-h % | plume× | spots |
|--:|:--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 30 | 14:30 | 136.0 | 47.2 | 27006 | 8.5 | 4.8 | 36 | 6.3 | 1.05 | 0 |
| 60 | 15:00 | 491.9 | 38.7 | 22943 | 7.9 | 4.8 | 36 | 6.5 | 1.05 | 0 |
| 90 | 15:30 | 815.8 | 33.6 | 17165 | 6.9 | 4.6 | 33 | 6.6 | 1.05 | 0 |
| 120 | 16:00 | 1211.8 | 30.8 | 18181 | 7.1 | 4.6 | 33 | 6.4 | 1.05 | 2 |

## Convective / pyro behaviour

- Continuous-Haines **5.5**, inverted-V False, surface dewpoint depression 18.9 °C, mid-level RH 40%, LCL 2357 m
- Plume loft factor (max) **1.05**
- Fire-type cell counts: {'0': 132, '1': 71, '2': 100}

## Fire Weather Index (Canadian FWI System)

Multi-week spin-up (SIR Toscana rain + ERA5) → fire-day gridded FWI. **Danger class: very high.**

| | FFMC | DMC | DC | ISI | BUI | FWI |
|:--|--:|--:|--:|--:|--:|--:|
| spin-up codes | 90.1 | 101.3 | 254.6 | – | – | – |
| fire-day AOI mean | 91.6 | 105.9 | 263.6 | 12.8 | 106.0 | 37.8 |

- Codes are fuel-layer moisture (FFMC fast litter → DMC duff → DC deep drought); indexes ISI/BUI/FWI are behaviour.
- FFMC fine fuel moisture 9.2% seeds the dead-fuel model; BUI 106.0 cures the live herb.
- Rasters + 6-panel map: `fwi/  (ffmc/dmc/dc/isi/bui/fwi/fine_fuel_moisture .tif + fwi_map.png)`

## Operative sector analysis (final)

```
Operative analysis at t=120 min (fire heading 235 deg):
  head             n=  28 ROS   85.3 ft/min | resultant 11.51@239 | driver: wind
      wind 12.11@235  slope  0.73@ 43  fuel  0.80@319
  right flank      n=  30 ROS   58.0 ft/min | resultant 12.65@248 | driver: wind
      wind 12.90@235  slope  0.81@ 45  fuel  2.77@320
  tail             n=  47 ROS   31.6 ft/min | resultant  8.74@234 | driver: wind
      wind  9.64@235  slope  0.86@ 55  fuel  0.12@130
  left flank       n=  45 ROS   36.0 ft/min | resultant 10.20@239 | driver: wind
      wind 10.37@235  slope  0.49@ 49  fuel  0.80@304
```