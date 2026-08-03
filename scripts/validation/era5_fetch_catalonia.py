"""Fetch ERA5 pressure-level + surface fields over Catalonia (for era5_regime_relationship.py).

ERA5 is the rich-level reference used to referee the ICON-EU pyroconvection product against
the Catalan fire-service (Bombers) product: 22 pressure levels (200-1000 hPa), far more than
ICON-2I's 5 or a coarse degrade, and an independent model (ECMWF reanalysis, not ICON).

Needs a Copernicus CDS account: `pip install cdsapi` and a ~/.cdsapirc with your key
(https://cds.climate.copernicus.eu/how-to-api). ERA5(T) lags ~5 days, so pick past dates --
the script errors with the latest available date if you ask for something too recent.

Usage:  python era5_fetch_catalonia.py [YYYY-MM] [day1,day2,...]
        (defaults: 2026-07, days 06..10 -- the most recent available near the Bombers period)
Writes era5_pl_cat.nc and era5_sfc_cat.nc into the current directory.
"""
import sys

import cdsapi

MONTH = sys.argv[1] if len(sys.argv) > 1 else "2026-07"
YEAR, MON = MONTH.split("-")
DAYS = sys.argv[2].split(",") if len(sys.argv) > 2 else ["06", "07", "08", "09", "10"]
AREA = [43.0, 0.0, 40.4, 3.4]                 # N, W, S, E (Catalonia)
TIMES = [f"{h:02d}:00" for h in (6, 9, 12, 15, 18)]
PLEV = ["1000", "975", "950", "925", "900", "875", "850", "825", "800", "775", "750",
        "700", "650", "600", "550", "500", "450", "400", "350", "300", "250", "200"]

c = cdsapi.Client()
print("pressure levels...", flush=True)
c.retrieve("reanalysis-era5-pressure-levels", {
    "product_type": "reanalysis",
    "variable": ["temperature", "relative_humidity", "u_component_of_wind",
                 "v_component_of_wind", "geopotential"],
    "pressure_level": PLEV, "year": YEAR, "month": MON, "day": DAYS,
    "time": TIMES, "area": AREA, "format": "netcdf",
}, "era5_pl_cat.nc")
print("single levels...", flush=True)
c.retrieve("reanalysis-era5-single-levels", {
    "product_type": "reanalysis",
    "variable": ["2m_temperature", "2m_dewpoint_temperature", "surface_pressure",
                 "10m_u_component_of_wind", "10m_v_component_of_wind",
                 "geopotential", "land_sea_mask"],
    "year": YEAR, "month": MON, "day": DAYS, "time": TIMES, "area": AREA, "format": "netcdf",
}, "era5_sfc_cat.nc")
print("DONE", flush=True)
