# Security Policy

## Supported Versions

pyflam is pre-1.0 and under active development; security fixes are applied to the
latest `main`. Pin a commit if you need stability.

| Version | Supported          |
| ------- | ------------------ |
| `main`  | :white_check_mark: |
| < 0.1   | :x:                |

## Reporting a Vulnerability

Please **do not open a public issue** for security-sensitive reports.

- Preferred: open a private advisory via GitHub
  ([Security → Report a vulnerability](https://github.com/crifod/pyflam/security/advisories/new)).
- Alternatively, email **cristiano.foderi@gmail.com** with details and, if
  possible, a minimal reproduction.

We aim to acknowledge reports within a few days and will keep you informed of
progress toward a fix. Once resolved, we are happy to credit you in the release
notes unless you prefer to remain anonymous.

## Scope notes

pyflam can download forecast/reanalysis data (GFS via Herbie, ICON-2I from the
MISTRAL open archive, ERA5 from the Copernicus CDS) and can read landscape files
(`.lcp`, GeoTIFF). Treat downloaded data and third-party landscape files as
untrusted input. Credentials (e.g. the CDS token in `~/.cdsapirc`) are read from
your environment and never transmitted by pyflam to anywhere other than the
respective data provider.
