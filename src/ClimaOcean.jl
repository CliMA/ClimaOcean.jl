module ClimaOcean

# Use the README as the module docs
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end ClimaOcean

export
    OceanSeaIceModel,
    FreezingLimitedOceanTemperature,
    Radiation,
    LatitudeDependentAlbedo,
    SimilarityTheoryFluxes,
    CoefficientBasedFluxes,
    MomentumRoughnessLength,
    ScalarRoughnessLength,
    ComponentInterfaces,
    SkinTemperature,
    BulkTemperature,
    PrescribedAtmosphere,
    JRA55PrescribedAtmosphere,
    JRA55NetCDFBackend,
    regrid_bathymetry,
    retrieve_bathymetry,
    Metadata,
    Metadatum,
    ECCOMetadatum,
    EN4Metadatum,
    ETOPO2022,
    ECCO2Daily, ECCO2Monthly, ECCO4Monthly,
    ECCO2DarwinMonthly, ECCO4DarwinMonthly,
    EN4Monthly,
    GLORYSDaily, GLORYSMonthly, GLORYSStatic,
    RepeatYearJRA55, MultiYearJRA55,
    first_date,
    last_date,
    all_dates,
    JRA55FieldTimeSeries,
    LinearlyTaperedPolarMask,
    DatasetRestoring,
    ocean_simulation,
    sea_ice_simulation,
    atmosphere_simulation,
    sea_ice_dynamics,
    initialize!

using Oceananigans
using Oceananigans.Operators: ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ
using DataDeps

@reexport using NumericalEarth
@reexport using NumericalEarth.DataWrangling
@reexport using NumericalEarth.EarthSystemModels
@reexport using NumericalEarth.EarthSystemModels.InterfaceComputations

#####
##### Source code
#####

include("Oceans/Oceans.jl")
include("SeaIces/SeaIces.jl")
include("InitialConditions/InitialConditions.jl")
include("Bathymetry/Bathymetry.jl")
include("Diagnostics/Diagnostics.jl")

using NumericalEarth.DataWrangling: ETOPO, ECCO, GLORYS, EN4, JRA55
using NumericalEarth.Bathymetry
using .InitialConditions
using NumericalEarth.EarthSystemModels
using NumericalEarth.Atmospheres
using .Oceans
using .SeaIces

end # module
