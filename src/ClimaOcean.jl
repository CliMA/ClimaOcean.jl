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
    JRA55_field_time_series,
    ECCO_field,
    ECCOMetadata,
    ECCORestoring,
    LinearlyTaperedPolarMask,
    Checkpointer,
    ocean_simulation,
    sea_ice_simulation,
    initialize!,
    @root,
    @onrank,
    @distribute,
    @handshake

using Oceananigans
using Oceananigans.Operators: ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ
using DataDeps

using Oceananigans.OutputReaders: GPUAdaptedFieldTimeSeries, FieldTimeSeries
using Oceananigans.Grids: node

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

include("DistributedUtils.jl")
include("OceanSimulations/OceanSimulations.jl")
include("SeaIceSimulations.jl")
include("OceanSeaIceModels/OceanSeaIceModels.jl")
include("VerticalGrids.jl")
include("InitialConditions/InitialConditions.jl")
include("DataWrangling/DataWrangling.jl")
include("Bathymetry.jl")
include("Diagnostics/Diagnostics.jl")
include("OutputWriters.jl")

using .DistributedUtils
using .VerticalGrids
using .Bathymetry
using .DataWrangling
using .InitialConditions
using .OceanSeaIceModels
using .OceanSimulations
using .SeaIceSimulations
using .DataWrangling: JRA55, ECCO
using .OutputWriters

using ClimaOcean.OceanSeaIceModels: PrescribedAtmosphere
using ClimaOcean.DataWrangling.JRA55: JRA55PrescribedAtmosphere, JRA55NetCDFBackend
using ClimaOcean.DataWrangling.ECCO

end # module

