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
    SkinTemperature,
    BulkTemperature,
    PrescribedAtmosphere,
    JRA55PrescribedAtmosphere,
    JRA55NetCDFBackend,
    regrid_bathymetry,
    retrieve_bathymetry,
    stretched_vertical_faces,
    exponential_z_faces,
    PowerLawStretching, LinearStretching,
    exponential_z_faces,
    Metadata,
    all_dates,
    JRA55FieldTimeSeries,
    ECCO_field, 
    ECCORestoring,
    LinearlyTaperedPolarMask,
    ocean_simulation,
    sea_ice_simulation,
    initialize!

using Oceananigans
using Oceananigans.Operators: ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ
using DataDeps

using Oceananigans.OutputReaders: GPUAdaptedFieldTimeSeries, FieldTimeSeries
using Oceananigans.Grids: node

# Does not really matter if there is no model
reference_density(::Nothing) = 0
heat_capacity(::Nothing) = 0

reference_density(unsupported) =
    throw(ArgumentError("Cannot extract reference density from $(typeof(unsupported))"))

heat_capacity(unsupported) =
    throw(ArgumentError("Cannot deduce the heat capacity from $(typeof(unsupported))"))

const SomeKindOfFieldTimeSeries = Union{FieldTimeSeries,
                                        GPUAdaptedFieldTimeSeries}

const SKOFTS = SomeKindOfFieldTimeSeries

@inline stateindex(a::Number, i, j, k, args...) = a
@inline stateindex(a::AbstractArray, i, j, k, args...) = @inbounds a[i, j, k]
@inline stateindex(a::SKOFTS, i, j, k, grid, time, args...) = @inbounds a[i, j, k, time]

@inline function stateindex(a::Function, i, j, k, grid, time, loc)
    LX, LY, LZ = loc 
    λ, φ, z = node(i, j, k, grid, LX(), LY(), LZ())
    return a(λ, φ, z, time)
end

@inline function stateindex(a::Tuple, i, j, k, grid, time)
    N = length(a)
    ntuple(Val(N)) do n
        stateindex(a[n], i, j, k, grid, time)
    end
end

@inline function stateindex(a::NamedTuple, i, j, k, grid, time)
    vals = stateindex(values(a), i, j, k, grid, time)
    names = keys(a)
    return NamedTuple{names}(vals)
end

include("OceanSimulations/OceanSimulations.jl")
include("SeaIceSimulations.jl")
include("OceanSeaIceModels/OceanSeaIceModels.jl")
include("VerticalGrids.jl")
include("InitialConditions/InitialConditions.jl")
include("DataWrangling/DataWrangling.jl")
include("Bathymetry.jl")
include("Diagnostics/Diagnostics.jl")

using .VerticalGrids
using .Bathymetry
using .DataWrangling
using .InitialConditions
using .OceanSeaIceModels
using .OceanSimulations
using .SeaIceSimulations
using .DataWrangling: JRA55, ECCO

using ClimaOcean.OceanSeaIceModels: PrescribedAtmosphere
using ClimaOcean.DataWrangling.JRA55: JRA55PrescribedAtmosphere, JRA55NetCDFBackend
using ClimaOcean.DataWrangling.ECCO

end # module

