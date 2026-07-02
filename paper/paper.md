---
title: 'pyflam: open, multiplatform wildfire-behavior modelling with coupled fire–atmosphere pyroconvection and fire-danger assessment'
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
assessment, and incident decision support. In practice this work is dominated by
closed-source, Windows-only desktop tools (FlamMap, FARSITE [@finney1998]) that are
difficult to script, automate, embed in larger systems, or run in the cloud or on
non-Windows platforms, which limits reproducible research and operational automation.
There is a need for an open, cross-platform, scriptable implementation of the established
operational product set that also connects fire behavior to live weather and to the
fire–atmosphere coupling the desktop tools omit — in particular the pyroconvection and
fire-danger assessment that governs the most dangerous fires. `pyflam` addresses this
need: it reproduces the FlamMap-style per-cell products through a Python API, drives them
from forecast and reanalysis data per cell, and adds coupled-atmosphere pyroconvection
diagnostics as selectable, literature-grounded components. Its intended users are
fire-behavior researchers and operational analysts in planning, incident management, and
suppression support.

# State of the field

Several open-source tools address parts of this space. For landscape spread, Cell2Fire
[@pais2021] and ELMFIRE [@lautenberger2013] provide fast fire-growth models; ForeFire
[@filippi2018; @filippi2025] offers front-tracking spread and, coupled to Meso-NH,
research-grade fire–atmosphere simulation; WindNinja [@forthofer2014] provides diagnostic
terrain winds; and WRF-SFIRE [@mandel2011] embeds a spread model in a full mesoscale
atmospheric model. Each is strong in its niche, but none packages the *operational
FlamMap product set* (surface behavior, crown-fire potential, minimum-travel-time growth,
and random-ignition burn probability, on community `.lcp`/`.fms` inputs) together with a
*weather-to-fire pipeline* and *routine pyroconvection/fire-danger mapping* in one
permissively licensed, scriptable Python library.

pyflam's build-versus-contribute rationale is therefore twofold. Where a trusted
operational standard exists (Rothermel surface spread, the classical crown-fire stack,
minimum-travel-time growth), pyflam re-implements it faithfully and cross-validates it
against the reference tool, rather than asking users to change paradigms. Where the
operational tools stop — weather-driven per-cell moisture, native terrain winds, a
Finsler-Eikonal spread alternative, a literature-current crown-fire model, physics-based
spotting, and vertically resolved pyroconvection — it contributes new, selectable
components built directly from the peer-reviewed literature. The distinctive contribution
is pyroconvection and fire-danger assessment: the operational desktop tools do not
diagnose plume-driven, potentially pyroCb-forming behavior, and the research-grade coupled
models that can (e.g. WRF-SFIRE) are heavy to run and not oriented toward routine danger
mapping. pyflam fills that gap with vertical-profile diagnostics [@peterson2017;
@castellnou2022; @tory2021] driven per cell by high-resolution regional forecasts such as
the convection-permitting ICON-2I model over Italy.

# Software design

`pyflam` is a pure-Python core on NumPy and SciPy, with heavier capabilities gated behind
optional extras (geospatial I/O, atmospheric forcing, and Numba JIT) and external engines
(OpenFOAM, Herbie) discovered at runtime and self-skipped when absent, so the core
installs and runs anywhere. A central design choice is a vectorized `SurfaceKernel`: the
wind- and slope-independent Rothermel terms are computed once per fuel and moisture state
and then applied to per-cell array inputs, which makes whole-landscape moisture
conditioning, weather forcing, and crown classification tractable at multi-million-cell
scale. Interchangeable back-ends are exposed as options rather than replacements — the
fire-growth engine (lattice minimum-travel-time Dijkstra, or a semi-Lagrangian
anisotropic-Eikonal solver), the wind model (mass-consistent diagnostic vs. buoyant RANS),
and the crown-spread model are all selectable. Two trade-offs are worth noting. First, the
fire–atmosphere coupling uses a *quasi-steady* buoyant-RANS plume re-solved on a march
interval rather than a transient large-eddy simulation; this trades some fidelity for
tractability and keeps the plume feedback usable operationally (and testable without
OpenFOAM through an injectable solver). Second, established models remain trusted defaults
while the novel methods are opt-in, so users can reproduce the reference tools exactly and
adopt the newer science by choice.

# Research impact statement

`pyflam` is a young project; its near-term significance rests on three points. First, it
lowers the barrier to reproducible fire-behavior analysis by making the operational
product set scriptable and cross-platform, with a validated surface core: cell-by-cell
agreement with a real FlamMap run reaches ~3% on surface rate of spread, ~1° on
maximum-spread direction, and ~2% on conditional fireline intensity over a 1.6-million-cell
landscape. Second, it provides a research platform for the fire–atmosphere coupling and
pyroconvection that the operational desktop tools omit, with each component traceable to
the peer-reviewed literature. Third, it is already deployed as an operational product: an
unattended daily pipeline produces pyroconvection and fire-danger maps over Italy from the
convection-permitting ICON-2I forecast, and an interactive graphical interface exposes the
same pipelines to non-programmer analysts. Broader validation of the novel components
against field observations and independent coupled models is ongoing and is the subject of
planned work. `pyflam` is developed openly with a substantial automated test suite run in
continuous integration across Python 3.11–3.13.

# Acknowledgements

We acknowledge the USDA Forest Service / Missoula Fire Sciences Laboratory, whose FlamMap
system established the operational paradigm that inspired this work; `pyflam` is an
independent implementation built from published, peer-reviewed models and contains no
FlamMap code. This work received no specific funding.

# AI usage disclosure

Parts of the `pyflam` implementation, and this paper, were produced with the assistance of
a generative-AI coding and research assistant (Claude, Anthropic). The fire-science models
are independent implementations of the cited peer-reviewed publications; the empirical
coefficients were checked against the original papers before implementation, and every
cited reference was verified against its authoritative record.

# References
