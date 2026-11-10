"""
ClimaOcean runs the Ocean Model Intercomparison Project (OMIP) simulations of NumericalEarth's ocean--sea ice model
at one, one quarter and one twelfth of a degree, and evaluates them against observations.
"""
module ClimaOcean

export OneDegree, QuarterDegree, TwelfthDegree, OMIPConfiguration,
       omip_simulation, add_omip_diagnostics!,
       TimeStepSchedule, TotalWaterNormalization, henyey_diffusivity,
       OMIPOutput, omip_figures,
       BoyerMontegutMonthly, HadISSTMonthly, NCEPMonthly

# Figures, implemented in the Makie extension
export sea_surface_temperature_bias, sea_surface_salinity_bias, sea_surface_height, mixed_layer_depth,
       sea_ice_concentration, sea_ice_concentration_bias, surface_fluxes, wind_stress,
       sea_surface_height_variability, sea_ice_extent, sea_ice_area, arctic_sea_ice_volume,
       sea_ice_area_time_series, arctic_sea_ice_volume_time_series, kinetic_energy, global_drift,
       mean_profiles, zonal_mean_sections, zonal_mean_drift, zonal_mean_mixed_layer_depth, drift_profiles,
       strait_transports, atlantic_overturning, near_surface_currents, equatorial_undercurrent,
       atlantic_overturning_at_rapid, sea_level_drift, content_conservation, barotropic_streamfunction,
       arctic_freshwater_transports, zonal_rms_bias, rms_bias_profiles

using Oceananigans
using Oceananigans.Units
using Oceananigans: fields
using Oceananigans.Architectures: architecture, cpu_architecture, on_architecture
using Oceananigans.Fields: ConstantField
using Oceananigans.Grids: λnode, φnode, znode, static_column_depthᶜᶜᵃ
using Oceananigans.ImmersedBoundaries: inactive_node
using Oceananigans.Models: buoyancy_operation
using Oceananigans.Models.HydrostaticFreeSurfaceModels: update_zstar_scaling!
using Oceananigans.Operators: Az, Axᶠᶜᶜ, Ayᶜᶠᶜ, Δzᶜᶜᶜ, Δz⁻¹ᶜᶠᶜ, δzᵃᵃᶜ, ℑxᶜᵃᵃ, ℑxᶠᵃᵃ, ℑyᵃᶜᵃ, ℑyᵃᶠᵃ
using Oceananigans.TurbulenceClosures: IsopycnalSkewSymmetricDiffusivity, κ_ϵSyᶜᶠᶠ
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEMixingLength
using Oceananigans.Utils: KernelParameters, launch!, prettysummary

using NumericalEarth
using NumericalEarth.DataWrangling: SurfaceFluxRestoring, ConservativeSurfaceFluxRestoring,
                                    ConservativeSurfaceFluxRestoringCallback
using NumericalEarth.DataWrangling.SeaWiFS: SeaWiFSMonthly
using NumericalEarth.Diagnostics: MixedLayerDepthField
using NumericalEarth.EarthSystemModels: AtmosphericThermodynamics
using NumericalEarth.EarthSystemModels.InterfaceComputations: LargeYeagerTransferCoefficients,
                                                              COARELogarithmicSimilarityProfile,
                                                              ConvectiveGustiness,
                                                              ImpureSaturationSpecificHumidity,
                                                              GillSaturationEnhancement,
                                                              ThreeEquationHeatFlux,
                                                              MomentumBasedFrictionVelocity,
                                                              atmosphere_sea_ice_stability_functions
using NumericalEarth.SeaIces: sea_ice_dynamics, default_snow_thermodynamics

using ClimaSeaIce: ConductiveFlux, IncrementalRemapping, UniformThicknessDistribution
using ConservativeRegridding: ConservativeRegridding
using SeawaterPolynomials.TEOS10: Θ_from_T
using WorldOceanAtlasTools: WorldOceanAtlasTools

using Dates: DateTime
using JLD2: ZstdFilter
using KernelAbstractions: @index, @kernel
using LazyArtifacts: @artifact_str
using Printf: @sprintf

include("Observations/Observations.jl")
using .Observations: BoyerMontegutMonthly, HadISSTMonthly, NCEPMonthly

include("absolute_salinity.jl")
include("initial_conditions.jl")
include("time_step_schedule.jl")
include("total_water.jl")
include("configurations.jl")
include("omip_simulation.jl")
include("diagnostics.jl")
include("omip_output.jl")
include("figures.jl")

end # module
