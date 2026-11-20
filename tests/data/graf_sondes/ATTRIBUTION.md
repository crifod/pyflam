# GRAF in-plume radiosonde environment sondes — test fixtures

Vendored ambient (environment) sonde profiles used as real-data regression fixtures
for the pyflam atmosphere pipeline (`tests/test_atmosphere_sondes.py`).

**Source:** Castellnou Ribau, M., et al. (2025). *In-Plume Radiosonde Fireline
Observations.* Zenodo. https://doi.org/10.5281/zenodo.15264835
**License:** CC BY 4.0. Redistributed here with attribution; only the ambient
"Environment" flight histories are included (a small subset for testing).

Companion to: Castellnou Ribau et al. (2025), *Atmospheric Measurement Techniques*
18(24), 7805–7831, doi:10.5194/amt-18-7805-2025.

Each CSV is a raw balloon flight history: UTC time, altitude (m MSL / m AGL),
pressure (Pa), GPS speed/heading (= wind), temperature (°C), relative humidity (%).
