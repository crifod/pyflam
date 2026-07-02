---
title: 'pyflam: open, multiplatform wildfire-behavior modelling with weather-driven and coupled fire–atmosphere capabilities'
tags:
  - Python
  - wildfire
  - fire behavior
  - fire spread
  - Rothermel
  - crown fire
  - fire–atmosphere coupling
  - pyroconvection
authors:
  - name: Cristiano Foderi
    orcid: 0000-0002-5474-5433
    affiliation: 1
affiliations:
  - name: Independent researcher
    index: 1
date: 2 July 2026
bibliography: paper.bib
---

# Summary

`pyflam` is an open-source, cross-platform Python library for wildland-fire-behavior
modelling in support of fire management, planning, and suppression. It implements the
Rothermel [@rothermel1972] surface fire-spread model with Albini's [@albini1976]
refinements and the two standard operational fuel-model sets [@anderson1982;
@scottburgan2005], turns per-cell fire behavior into a directional spread field, and
grows fires with the minimum-travel-time (MTT) method [@finney2002]. From a landscape
and a set of ignitions it produces the per-cell product set that desktop systems such
as FlamMap [@finney2006] established as operationally useful — surface rate of spread,
fireline intensity and flame length; crown-fire potential; fire arrival time; and
random-ignition burn probability — and reads and writes the community `.lcp`
landscape, `.fms` moisture, GeoTIFF and GeoJSON formats so it fits existing workflows.

`pyflam` delivers this as scriptable, automatable, dependency-light software that runs
on any operating system (a pure-Python core on NumPy and SciPy, with optional
geospatial, atmospheric and JIT-acceleration extras). On top of the established core it
adds capabilities as *selectable* components implemented directly from the peer-reviewed
literature: weather-driven, per-cell dead-fuel-moisture conditioning from live forecast
or reanalysis data — global GFS, the convection-permitting regional forecast ICON-2I
(2.2 km, Italy), and ERA5 reanalysis; native terrain-wind solvers following the
WindNinja science [@forthofer2014]; an anisotropic-Eikonal (Randers–Finsler) spread
solver offered alongside MTT [@sethian2003; @mirebeau2014]; a literature-current
crown-fire model [@cruz2005; @cruz2004] provided as an alternative to the classic
operational stack; and physics-based ember spotting. Where canopy fuel inputs are
unavailable, `pyflam` can assemble a crown-fire-ready landscape from remotely sensed,
GEDI-calibrated canopy height (the Meta/WRI High-Resolution Canopy Height product;
@tolan2024, open data), deriving the canopy base height and bulk density the crown model
requires as clearly-flagged estimates.

A central emphasis of `pyflam` is **pyroconvection and fire-danger assessment**. The
coupling is built on a buoyant plume that modifies the wind driving the fire, closing a
crown → plume → wind → crown feedback. On top of it, `pyflam` implements vertical-profile
pyroconvection diagnostics keyed on the atmospheric geometry that the pyrocumulonimbus
(pyroCb) literature identifies as decisive — a deep, dry, well-mixed boundary layer (high
lifting condensation level) capped by moisture aloft, rather than high surface CAPE
[@peterson2017; @castellnou2022]. These include the lifting condensation level, the
Continuous Haines index [@mills2010], an "inverted-V" sounding detector, a
pyroconvection-type classifier implementing the @castellnou2022 procedure, and a pyroCb
firepower threshold [@tory2021]. Driven per cell by GFS or the high-resolution ICON-2I fields,
they produce spatial pyroconvection and fire-danger maps and feed a profile-aware plume
factor back into the coupled spread. This mapping pipeline is packaged to run unattended
on a daily schedule over Italy, driven by the operational ICON-2I forecast, and is also
exposed through an interactive graphical interface. Established models are retained as
trusted defaults; the additions are exposed as options so users can compare them.

# Statement of need

Landscape fire-behavior modelling is central to fuel-treatment planning, risk
assessment, and incident decision support, but the most widely used desktop tools
(FlamMap, FARSITE [@finney1998]) are closed-source, Windows-only binaries that are
difficult to script, automate, embed in larger systems, or run in the cloud or on
non-Windows platforms. This limits reproducible research and operational automation.

Open-source alternatives exist and target parts of this space — for example Cell2Fire
[@pais2021] and ELMFIRE [@lautenberger2013] for landscape spread, ForeFire
[@filippi2018] for spread and coupling, WindNinja [@forthofer2014] for diagnostic
wind, and WRF-SFIRE [@mandel2011] for research-grade coupled fire–atmosphere
simulation. `pyflam`'s contribution is to package the *operational FlamMap-style product
set* and a *weather-to-fire pipeline* together in one permissively installable Python
API, and to make the newer science available as transparent, selectable options next to
the classical defaults. Its most distinctive capability is **pyroconvection and
fire-danger assessment**: the operational tools listed above do not diagnose the
plume-driven, potentially pyroCb-forming behavior that governs the most dangerous fires,
and the research-grade coupled models that do (e.g. WRF-SFIRE) are heavy to run and not
oriented toward routine danger mapping. `pyflam` fills that gap with vertical-profile
pyroconvection diagnostics grounded in the current literature [@peterson2017;
@castellnou2022; @tory2021] and driven per cell by high-resolution regional forecasts —
notably the convection-permitting ICON-2I model over Italy — so an analyst can produce
spatial pyroconvection-potential and fire-danger maps and couple them into fire spread
within one reproducible pipeline. Together with the Finsler-Eikonal spread solver,
per-cell weather-driven moisture, and a literature-current crown-fire model, this lowers
the barrier to reproducible fire-behavior analysis and provides a research platform for
the fire–atmosphere coupling that operational tools approximate or omit.

The deterministic surface core is cross-validated cell-by-cell against a FlamMap run on
a real 1.6-million-cell landscape (surface rate of spread within ~3%, maximum-spread
direction within ~1°, conditional fireline intensity within ~2%). The added components
are implemented and verified against their published equations and, where applicable,
analytic benchmarks; broader validation of these components against field observations
and independent coupled models is ongoing and is the subject of future work. `pyflam`
is tested with a large automated suite run in continuous integration across Python
3.11–3.13.

# Acknowledgements

We acknowledge the USDA Forest Service / Missoula Fire Sciences Laboratory, whose
FlamMap system established the operational paradigm that inspired this work; `pyflam` is
an independent implementation built from published, peer-reviewed models and contains no
FlamMap code. This work received no specific funding.

# References
