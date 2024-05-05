module ClimaOcean

export
    OceanSeaIceModel,
    MinimumTemperatureSeaIce,
    Radiation,
    JRA55_prescribed_atmosphere,
    NetCDFBackend,
    ecco4_field,
    regrid_bathymetry,
    stretched_vertical_faces,
    PowerLawStretching, LinearStretching,
    jra55_field_time_series,
    ecco4_field, ECCOMetadata,
    initialize!

using Oceananigans
using Oceananigans.Operators: ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ
using DataDeps

include("OceanSeaIceModels/OceanSeaIceModels.jl")
include("VerticalGrids.jl")
include("InitialConditions/InitialConditions.jl")
include("DataWrangling/DataWrangling.jl")
include("Bathymetry.jl")
include("Diagnostics.jl")
include("OceanSimulations/OceanSimulations.jl")

using .VerticalGrids
using .Bathymetry
using .InitialConditions
using .OceanSeaIceModels: OceanSeaIceModel
using .OceanSimulations
using .DataWrangling
using ClimaOcean.DataWrangling: NetCDFBackend
using ClimaOcean.DataWrangling: JRA55_prescribed_atmosphere
using ClimaOcean.DataWrangling: ecco4_field

using .OceanSeaIceModels: OceanSeaIceModel, Radiation

end # module

