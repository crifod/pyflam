"""pyflam: a Python reimplementation of FlamMap's fire-behavior science.

Covers the Rothermel (1972) surface fire spread model and the standard fuel
models (steps 1/1b), landscape I/O and whole-landscape behavior (step 2), the
mass-consistent and RANS terrain-wind solvers (steps 2b/2c), crown fire (step 4)
and the directional Minimum Travel Time spread/perimeter engine (steps 3/5). See
README.md for the roadmap.
"""

from __future__ import annotations

from . import (
    atmosphere, cfd, crownfire, fuel_conditioning, fuel_models, io_lcp,
    landscape, meteo_report, mtt, nrt, operative, pyroconvection, rothermel,
    spotting, units, validate, wind, wind_reduction, windsolver,
)
from .nrt import RunProduct, run_realtime
from .meteo_report import MeteoReport, meteo_report as build_meteo_report
from .operative import (
    OperativeReport, analyze_perimeter, driver_fields, to_geojson as
    operative_geojson,
)
from .atmosphere import (
    AtmosphericState, AtmosphericProfile, AtmosphereProvider, ConstantAtmosphere,
    DeadFuelMoistureModel, GriddedAtmosphere, PyroconvectionPotential,
    continuous_haines, inverted_v, lcl_height_m, open_atmosphere,
    pyroconvection_potential, spread_inputs_from_state,
)
from .fuel_conditioning import (
    canopy_transmission, condition_dead_fuel_moisture, condition_from_weather,
    dead_fuel_moisture_vpd, equation_of_time, solar_position, sun_exposure,
    terrain_insolation_factor, vapour_pressure_deficit,
)
from .pyroconvection import (
    couple_fire_wind, fire_atmosphere_march, fire_heat_flux, merge_plume_wind,
    superpose_plume,
)
from .spotting import FirebrandPhysics, SpottingModel
from .wind_reduction import (
    midflame_field,
    waf_field,
    wind_adjustment_factor,
)
from .fuel_models import (
    ALL as ALL_FUEL_MODELS,
    FuelModel,
    STANDARD_13,
    STANDARD_40,
    get as get_fuel_model,
)
from .landscape import Landscape, basic_fire_behavior
from .rothermel import FireBehavior, SurfaceKernel, spread, surface_kernel
from .wind import WindField, read_esri_ascii
from .windsolver import solve_mass_consistent, wind_field_from_landscape
from .cfd import solve_rans
from .crownfire import (
    CrownAwareField,
    CrownFireBehavior,
    active_crown_ros_cruz,
    crown_fire_behavior,
    crown_fire_potential,
    crown_fire_probability,
    crown_spread_field,
    crowning_index,
    torching_index,
    wind_20ft_to_u10_kmh,
)
from .mtt import (
    BurnProbabilityResult,
    SpreadField,
    anisotropic_eikonal,
    burn_probability,
    ignition_from_xy,
    minimum_travel_time,
    perimeter_mask,
    spread_field,
    spread_perimeter,
    spread_with_spotting,
)

__all__ = [
    "atmosphere",
    "meteo_report",
    "MeteoReport",
    "build_meteo_report",
    "operative",
    "OperativeReport",
    "analyze_perimeter",
    "driver_fields",
    "operative_geojson",
    "nrt",
    "run_realtime",
    "RunProduct",
    "AtmosphericState",
    "AtmosphereProvider",
    "ConstantAtmosphere",
    "GriddedAtmosphere",
    "DeadFuelMoistureModel",
    "AtmosphericProfile",
    "PyroconvectionPotential",
    "pyroconvection_potential",
    "continuous_haines",
    "inverted_v",
    "lcl_height_m",
    "open_atmosphere",
    "spread_inputs_from_state",
    "fuel_conditioning",
    "condition_dead_fuel_moisture",
    "condition_from_weather",
    "sun_exposure",
    "solar_position",
    "equation_of_time",
    "terrain_insolation_factor",
    "canopy_transmission",
    "vapour_pressure_deficit",
    "dead_fuel_moisture_vpd",
    "fuel_models",
    "io_lcp",
    "landscape",
    "rothermel",
    "units",
    "wind",
    "windsolver",
    "cfd",
    "crownfire",
    "mtt",
    "pyroconvection",
    "couple_fire_wind",
    "fire_heat_flux",
    "fire_atmosphere_march",
    "merge_plume_wind",
    "superpose_plume",
    "spotting",
    "validate",
    "wind_reduction",
    "wind_adjustment_factor",
    "waf_field",
    "midflame_field",
    "SpottingModel",
    "FirebrandPhysics",
    "spread_with_spotting",
    "FuelModel",
    "STANDARD_13",
    "STANDARD_40",
    "ALL_FUEL_MODELS",
    "get_fuel_model",
    "FireBehavior",
    "SurfaceKernel",
    "spread",
    "surface_kernel",
    "Landscape",
    "basic_fire_behavior",
    "WindField",
    "read_esri_ascii",
    "solve_mass_consistent",
    "wind_field_from_landscape",
    "solve_rans",
    "CrownFireBehavior",
    "CrownAwareField",
    "crown_fire_behavior",
    "crown_fire_potential",
    "crown_fire_probability",
    "crown_spread_field",
    "active_crown_ros_cruz",
    "wind_20ft_to_u10_kmh",
    "crowning_index",
    "torching_index",
    "SpreadField",
    "spread_field",
    "minimum_travel_time",
    "anisotropic_eikonal",
    "spread_perimeter",
    "burn_probability",
    "BurnProbabilityResult",
    "perimeter_mask",
    "ignition_from_xy",
]

__version__ = "0.1.0"
