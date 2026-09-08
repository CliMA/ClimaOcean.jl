module OMIPConfigurations

using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: znode, Face, static_column_depthᶜᶜᵃ
using Oceananigans.Operators: Δxᶜᶜᶜ, Δyᶜᶜᶜ
using Dates
using NCDatasets
using CUDA

using NumericalEarth
using NumericalEarth.Oceans: ocean_simulation, default_ocean_closure, TwoColorRadiation
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities:
    CATKEVerticalDiffusivity, CATKEMixingLength, CATKEEquation
using NumericalEarth.SeaIces: sea_ice_simulation, LandfastBasalStress, ThicknessDependentConductivity
using NumericalEarth.EarthSystemModels: OceanSeaIceModel,
    SimilarityTheoryFluxes,
    LinearStableStabilityFunction,
    MomentumBasedFrictionVelocity,
    ThreeEquationHeatFlux

using NumericalEarth.EarthSystemModels.InterfaceComputations:
    ComponentInterfaces,
    CoefficientBasedFluxes,
    COARELogarithmicSimilarityProfile,
    LargeYeagerTransferCoefficients,
    MomentumRoughnessLength,
    ScalarRoughnessLength,
    WindDependentWaveFormulation,
    TemperatureDependentAirViscosity,
    SimilarityScales,
    FixedIterations,
    large_yeager_stability_functions,
    atmosphere_sea_ice_stability_functions

using NumericalEarth.Radiations: SeaIceAlbedo, SurfaceRadiationProperties

using NumericalEarth.Bathymetry: regrid_bathymetry, ORCAGrid
using NumericalEarth.DataWrangling.SeaWiFS: SeaWiFSMonthly
using NumericalEarth.DataWrangling: Metadatum, Metadata, DatasetRestoring,
                                    SurfaceFluxRestoring,
                                    EN4Monthly, ECCO4Monthly
using NumericalEarth.DataWrangling.WOA: WOAMonthly
using NumericalEarth.DataWrangling.ORCA: ORCADataset, ORCAOne, ORCAQuarter, ORCATwelfth
using NumericalEarth.DataWrangling.JRA55: MultiYearJRA55, RepeatYearJRA55,
                                          JRA55PrescribedAtmosphere,
                                          JRA55PrescribedRadiation,
                                          JRA55PrescribedLand
using NumericalEarth.Diagnostics: MixedLayerDepthField

export omip_simulation,
       ThicknessDependentConductivity,
       add_omip_diagnostics!,
       add_ke_spectrum_diagnostic!,
       compute_report_fields,
       compute_woa_bias,
       strait_transports,
       strait_freshwater_transports,
       strait_overflow_transports,
       strait_sections,
       StraitSection,
       woa_to_teos10!,
       woa_salinity_fts_to_teos10!,
       KPPVerticalDiffusivity, KPPParameters,
       NEMOTKEVerticalDiffusivity, NEMOTKEParameters,
       NORiBaseVerticalDiffusivity,
       BoundaryValueTransport

# Patches to Oceananigans
include("oceananigans_patches.jl")

include("KPP/KPP.jl")

using .KPP: KPPVerticalDiffusivity, KPPParameters

include("NEMOTKE/NEMOTKE.jl")

using .NEMOTKE: NEMOTKEVerticalDiffusivity, NEMOTKEParameters

include("atmosphere.jl")
include("jra55_data_staging.jl")
include("omip_diagnostics.jl")
include("ke_spectrum_diagnostic.jl")
include("nemo_eddy_coefficients.jl")
include("cesm_eddy_coefficients.jl")
include("hybrid_eddy_coefficients.jl")
include("boundary_value_transport.jl")
include("mixed_layer_tapering.jl")
include("triad_slope_tapering.jl")
include("bottom_boundary_layer.jl")
include("advective_bottom_boundary_layer.jl")
include("overflow_restoring.jl")
include("labrador_restoring.jl")
include("omip_simulation.jl")
include("report_fields.jl")
include("strait_transports.jl")

end # module
