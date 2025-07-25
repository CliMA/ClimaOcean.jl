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
    initialize!

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

@inline function stateindex(a::Function, i, j, k, grid, time, (LX, LY, LZ), args...)
    λ, φ, z = node(i, j, k, grid, LX(), LY(), LZ())
    return a(λ, φ, z, time)
end

@inline function stateindex(a::Tuple, i, j, k, grid, time, args...)
    N = length(a)
    ntuple(Val(N)) do n
        stateindex(a[n], i, j, k, grid, time, args...)
    end
end

@inline function stateindex(a::NamedTuple, i, j, k, grid, time, args...)
    vals = stateindex(values(a), i, j, k, grid, time, args...)
    names = keys(a)
    return NamedTuple{names}(vals)
end

include("OceanSimulations/OceanSimulations.jl")
include("SeaIceSimulations.jl")
include("OceanSeaIceModels/OceanSeaIceModels.jl")
include("InitialConditions/InitialConditions.jl")
include("DataWrangling/DataWrangling.jl")
include("Bathymetry.jl")
include("Diagnostics/Diagnostics.jl")

using .DataWrangling
using .DataWrangling: ETOPO, ECCO, Copernicus, EN4, JRA55
using .Bathymetry
using .InitialConditions
using .OceanSeaIceModels
using .OceanSimulations
using .SeaIceSimulations

using ClimaOcean.OceanSeaIceModels: PrescribedAtmosphere, ComponentInterfaces, MomentumRoughnessLength, ScalarRoughnessLength
using ClimaOcean.DataWrangling.ETOPO
using ClimaOcean.DataWrangling.ECCO
using ClimaOcean.DataWrangling.Copernicus
using ClimaOcean.DataWrangling.EN4
using ClimaOcean.DataWrangling.JRA55
using ClimaOcean.DataWrangling.JRA55: JRA55NetCDFBackend

using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    Nx, Ny, Nz = 32, 32, 10
    @compile_workload begin
        depth = 6000
        z = Oceananigans.Grids.ExponentialCoordinate(Nz, -depth, 0)
        grid = Oceananigans.OrthogonalSphericalShellGrids.TripolarGrid(CPU(); size=(Nx, Ny, Nz), halo=(7, 7, 7), z)
        grid = ImmersedBoundaryGrid(grid, GridFittedBottom((x, y) -> -5000))
        # ocean = ocean_simulation(grid)
        # model = OceanSeaIceModel(ocean)
    end
end

end # module
