module OMIPConfigurations

using Oceananigans
using Oceananigans.Units
using Dates
using JLD2
using Statistics: mean

using NumericalEarth.Oceans: ocean_simulation
using NumericalEarth.SeaIces: sea_ice_simulation
using NumericalEarth.Atmospheres: JRA55PrescribedAtmosphere
using NumericalEarth.EarthSystemModels: OceanSeaIceModel, Radiation
using NumericalEarth.DataWrangling: Metadatum, Metadata, DatasetRestoring,
                                    EN4Monthly, ECCO4Monthly, WOAMonthly,
                                    MultiYearJRA55, JRA55NetCDFBackend

using ..OceanConfigurations: half_degree_tripolar_ocean, orca_ocean
using ..SeaIceConfigurations: half_degree_tripolar_sea_ice, orca_sea_ice
using ..Diagnostics: MixedLayerDepthField, compute_zonal_averages

export omip_simulation, add_omip_diagnostics!

include("atmosphere.jl")
include("diagnostics.jl")
include("transport_diagnostics.jl")
include("omip_simulation.jl")

end # module
