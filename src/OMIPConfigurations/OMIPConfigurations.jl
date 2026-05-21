module OMIPConfigurations

using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: znode, Face
using Dates
using NCDatasets
using CUDA

using NumericalEarth
using NumericalEarth.Oceans: ocean_simulation, default_ocean_closure
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities:
    CATKEVerticalDiffusivity, CATKEMixingLength, CATKEEquation
using NumericalEarth.SeaIces: sea_ice_simulation
using NumericalEarth.EarthSystemModels: OceanSeaIceModel,
    SimilarityTheoryFluxes,
    MomentumBasedFrictionVelocity,
    ThreeEquationHeatFlux

using NumericalEarth.EarthSystemModels.InterfaceComputations:
    ComponentInterfaces,
    CoefficientBasedFluxes,
    COARELogarithmicSimilarityProfile,
    LargeYeagerTransferCoefficients,
    LinearStableStabilityFunction,
    MomentumRoughnessLength,
    ScalarRoughnessLength,
    WindDependentWaveFormulation,
    TemperatureDependentAirViscosity,
    SimilarityScales,
    FixedIterations,
    large_yeager_stability_functions,
    atmosphere_sea_ice_stability_functions

using NumericalEarth.Radiations: SeaIceAlbedo, SurfaceRadiationProperties

using NumericalEarth.DataWrangling: Metadatum, Metadata, DatasetRestoring,
                                    SurfaceFluxRestoring,
                                    EN4Monthly, ECCO4Monthly, WOAAnnual
using NumericalEarth.DataWrangling.WOA: WOAMonthly
using NumericalEarth.DataWrangling.JRA55: MultiYearJRA55, RepeatYearJRA55,
                                          JRA55PrescribedAtmosphere,
                                          JRA55PrescribedRadiation,
                                          JRA55PrescribedLand
using NumericalEarth.Diagnostics: MixedLayerDepthField

using ..OceanConfigurations: half_degree_tripolar_ocean,
                             tenth_degree_tripolar_ocean,
                             orca_ocean

export omip_simulation,
       add_omip_diagnostics!,
       add_ke_spectrum_diagnostic!,
       compute_report_fields,
       compute_woa_bias,
       strait_transports,
       strait_sections,
       StraitSection,
       woa_to_teos10!,
       woa_salinity_fts_to_teos10!,
       KPPVerticalDiffusivity, KPPParameters,
       NEMOTKEVerticalDiffusivity, NEMOTKEParameters,
       NORiBaseVerticalDiffusivity

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
include("omip_simulation.jl")
include("report_fields.jl")
include("strait_transports.jl")

end # module
