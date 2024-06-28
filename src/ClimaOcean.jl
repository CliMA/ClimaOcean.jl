module ClimaOcean

export
    OceanSeaIceModel,
    MinimumTemperatureSeaIce,
    Radiation,
    SimilarityTheoryTurbulentFluxes,
    JRA55_prescribed_atmosphere,
    JRA55NetCDFBackend,
    ecco2_field,
    regrid_bathymetry,
    retrieve_bathymetry,
    stretched_vertical_faces,
    exponential_z_faces,
    PowerLawStretching, LinearStretching,
    jra55_field_time_series,
    ocean_simulation,
    ecco2_field, ECCO2Metadata,
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
using .DataWrangling: JRA55
using .DataWrangling: ECCO2
using .InitialConditions
using .OceanSeaIceModels: OceanSeaIceModel
using .OceanSimulations
using .DataWrangling: JRA55, ECCO2
using ClimaOcean.DataWrangling.JRA55: JRA55_prescribed_atmosphere, JRA55NetCDFBackend
using ClimaOcean.DataWrangling.ECCO2: ecco2_field

using .OceanSeaIceModels: OceanSeaIceModel, Radiation

end # module

