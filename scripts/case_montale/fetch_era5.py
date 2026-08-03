"""Fetch ERA5 for the Montale fire (16 Jul 2017) — surface + vertical profile.

Local time is CEST (UTC+2); the propagation afternoon 13:00-18:00 local = 11:00-16:00
UTC. Area box around Montale (43.96 N, 11.04 E). Two datasets:
  * single-levels  (via pyflam.fetch_era5): 2 m T/Td, 10 m wind, CAPE, BLH, heat fluxes
  * pressure-levels (via cdsapi): T, RH, geopotential, u, v, q -> the sounding for
    the pyroconvection fireABL / decoupling / classification.
"""
import os, sys, warnings
warnings.simplefilter("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from pyflam.atmosphere import fetch_era5

OUT = os.path.join(os.path.dirname(__file__), "..", "..", "docs", "case_montale_2017", "atmosphere")
os.makedirs(OUT, exist_ok=True)
AREA = (44.2, 10.8, 43.7, 11.3)          # N, W, S, E around Montale
TIMES = ["11:00", "13:00", "15:00"]       # UTC (= 13:00/15:00/17:00 local)


def main():
    # --- single levels (pyflam) ---
    sl = os.path.join(OUT, "era5_montale_20170716_single.nc")
    print("fetching ERA5 single-levels ...")
    try:
        fetch_era5(sl, date="2017-07-16", time=TIMES, area=AREA)
        print(f"  OK -> {sl}")
    except Exception as e:
        print(f"  single-level fetch failed: {str(e)[:150]}")

    # --- pressure levels (cdsapi direct) ---
    pl = os.path.join(OUT, "era5_montale_20170716_pressure.nc")
    print("fetching ERA5 pressure-levels ...")
    try:
        import cdsapi
        c = cdsapi.Client()
        c.retrieve("reanalysis-era5-pressure-levels", {
            "product_type": "reanalysis", "format": "netcdf",
            "variable": ["temperature", "relative_humidity", "geopotential",
                         "u_component_of_wind", "v_component_of_wind", "specific_humidity"],
            "pressure_level": ["1000", "975", "950", "925", "900", "850", "800", "750",
                               "700", "650", "600", "550", "500", "450", "400", "350",
                               "300", "250", "200", "150", "100"],
            "year": "2017", "month": "07", "day": "16", "time": TIMES,
            "area": [AREA[0], AREA[1], AREA[2], AREA[3]],
        }, pl)
        print(f"  OK -> {pl}")
    except Exception as e:
        print(f"  pressure-level fetch failed: {str(e)[:150]}")


if __name__ == "__main__":
    main()
