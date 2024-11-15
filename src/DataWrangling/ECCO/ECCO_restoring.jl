using Oceananigans.Grids: node, on_architecture
using Oceananigans.Fields: interpolate!, interpolate, location, instantiated_location
using Oceananigans.OutputReaders: Cyclical, TotallyInMemory, AbstractInMemoryBackend, FlavorOfFTS, time_indices
using Oceananigans.Utils: Time

using Base

using NCDatasets
using JLD2 

using Dates: Second
using ClimaOcean: stateindex
using ClimaOcean.DataWrangling: NearestNeighborInpainting

import Oceananigans.Fields: set!
import Oceananigans.Forcings: regularize_forcing
import Oceananigans.OutputReaders: new_backend, update_field_time_series!

@inline instantiate(T::DataType) = T()
@inline instantiate(T) = T

struct ECCONetCDFBackend{N, C, I, M} <: AbstractInMemoryBackend{Int}
    start :: Int
    length :: Int
    inpainting :: I
    metadata :: M

    function ECCONetCDFBackend{N, C}(start::Int, length::Int, inpainting, metadata) where {N, C}
        M = typeof(metadata)
        I = typeof(inpainting)
        return new{N, C, I, M}(start, length, inpainting, metadata)
    end
end

Adapt.adapt_structure(to, b::ECCONetCDFBackend{N, C}) where {N, C} = ECCONetCDFBackend{N, C}(b.start, b.length, nothing, nothing)

"""
    ECCONetCDFBackend(length; on_native_grid = false, inpainting = NearestNeighborInpainting(Inf))

Represent an ECCO FieldTimeSeries backed by ECCO native netCDF files.
Each time instance is stored in an individual file.
"""
function ECCONetCDFBackend(length, metadata;
                           on_native_grid = false, 
                           cache_inpainted_data = false,
                           inpainting = NearestNeighborInpainting(Inf))

    return ECCONetCDFBackend{on_native_grid, cache_inpainted_data}(1, length, inpainting, metadata)
end

Base.length(backend::ECCONetCDFBackend)  = backend.length
Base.summary(backend::ECCONetCDFBackend) = string("ECCONetCDFBackend(", backend.start, ", ", backend.length, ")")

const ECCONetCDFFTS{N} = FlavorOfFTS{<:Any, <:Any, <:Any, <:Any, <:ECCONetCDFBackend{N}} where N

new_backend(b::ECCONetCDFBackend{native, cache_data}, start, length) where {native, cache_data} =
    ECCONetCDFBackend{native}(start, length, b.inpainting, b.metadata)

on_native_grid(::ECCONetCDFBackend{native}) where native = native
cache_inpainted_data(::ECCONetCDFBackend{native, cache_data}) where {native, cache_data} = cache_data

function set!(fts::ECCONetCDFFTS) 
    backend = fts.backend
    start   = backend.start
    inpainting = backend.inpainting
    cache_data = cache_inpainted_data(backend)
    len = backend.length

    for t in start:start+len-1
        # Set each element of the time-series to the associated file
        metadatum = @inbounds backend.metadata[t] 
        set!(fts[t], metadatum; inpainting, cache_inpainted_data = cache_data)
    end

    fill_halo_regions!(fts)

    return nothing
end

"""
    ECCO_times(metadata; start_time = metadata.dates[1])

Extract the time values from the given metadata and calculates the time difference
from the start time.

Arguments
=========
- `metadata`: The metadata containing the date information.
- `start_time`: The start time for calculating the time difference. Defaults to the first date in the metadata.

Returns
=======
An array of time differences in seconds.
"""
function ECCO_times(metadata; start_time = first(metadata).dates)
    times = zeros(length(metadata))
    for (t, data) in enumerate(metadata)
        date = data.dates
        time = date - start_time
        time = Second(time).value
        times[t] = time
    end

    return times
end

"""
    ECCO_field_time_series(metadata::ECCOMetadata;
                           grid = nothing,
                           architecture = isnothing(grid) ? CPU() : architecture(grid),
                           time_indices_in_memory = 2,
                           time_indexing = Cyclical(),
                           inpainting_iterations = prod(size(metadata)),

Create a field time series object for ECCO data.

Arguments
=========

- `metadata`: `ECCOMetadata` containing information about the ECCO dataset.

Keyword Arguments
=================

- `grid`: where ECCO data is interpolated. If `nothing`, the native `ECCO` grid is used.

- `architecture`: where data is stored. Should only be set if `isnothing(grid)`.

- `time_indices_in_memory`: The number of time indices to keep in memory. Default: 2.

- `time_indexing`: The time indexing scheme to use. Default: `Cyclical()`.

- `inpainting`: The inpainting algorithm to use for ECCO interpolation.
                The only option is `NearestNeighborInpainting(maxiter)`, 
                where an average of the valid surrounding values is used `maxiter` times.

- `cache_inpainted_data`: If `true`, the data is cached to disk after inpainting for later retrieving. 
                          Default: `false`.

"""
function ECCO_field_time_series(metadata::ECCOMetadata;	
                                architecture = CPU(),	
                                time_indices_in_memory = 2,	
                                time_indexing = Cyclical(),
                                inpainting = NearestNeighborInpainting(prod(size(metadata))),
                                cache_inpainted_data = false,
                                grid = nothing)	

    # Make sure all the required individual files are downloaded
    download_dataset(metadata)

    inpainting isa Int && (inpainting = NearestNeighborInpainting(inpainting))

    ftmp = empty_ECCO_field(first(metadata); architecture)
    on_native_grid = isnothing(grid)
    on_native_grid && (grid = ftmp.grid)

    times = ECCO_times(metadata)
    loc = LX, LY, LZ = location(metadata)
    boundary_conditions = FieldBoundaryConditions(grid, loc)

    backend = ECCONetCDFBackend(time_indices_in_memory, metadata; on_native_grid, inpainting, cache_inpainted_data)

    fts = FieldTimeSeries{LX, LY, LZ}(grid, times; backend, time_indexing, boundary_conditions)
    set!(fts)	

    return fts	
end

ECCO_field_time_series(variable_name::Symbol, version=ECCO4Monthly(); kw...) = 
    ECCO_field_time_series(ECCOMetadata(variable_name, all_ECCO_dates(version), version); kw...)

# Variable names for restoreable data
struct Temperature end
struct Salinity end
struct UVelocity end
struct VVelocity end

oceananigans_fieldname = Dict(
    :temperature => Temperature(), 
    :salinity    => Salinity(), 
    :u_velocity  => UVelocity(), 
    :v_velocity  => VVelocity())

@inline Base.getindex(fields, i, j, k, ::Temperature) = @inbounds fields.T[i, j, k]
@inline Base.getindex(fields, i, j, k, ::Salinity)    = @inbounds fields.S[i, j, k]
@inline Base.getindex(fields, i, j, k, ::UVelocity)   = @inbounds fields.u[i, j, k]
@inline Base.getindex(fields, i, j, k, ::VVelocity)   = @inbounds fields.v[i, j, k]

Base.summary(::Temperature) = "temperature"
Base.summary(::Salinity)    = "salinity"
Base.summary(::UVelocity)   = "u_velocity"
Base.summary(::VVelocity)   = "v_velocity"

struct ECCORestoring{FTS, G, M, V, N}
    field_time_series :: FTS
    grid :: G
    mask :: M
    variable_name :: V
    rate :: N
end

Adapt.adapt_structure(to, p::ECCORestoring) = ECCORestoring(Adapt.adapt(to, p.field_time_series),
                                                            Adapt.adapt(to, p.grid),
                                                            Adapt.adapt(to, p.mask),
                                                            Adapt.adapt(to, p.variable_name),
                                                            Adapt.adapt(to, p.rate))

@inline function (p::ECCORestoring)(i, j, k, grid, clock, fields)

    # Figure out all the inputs: time, location, and node
    time = Time(clock.time)
    loc = location(p.field_time_series)

    # Possibly interpolate ECCO data from the ECCO grid to simulation grid.
    # Otherwise, simply extract the pre-interpolated data from p.field_time_series.
    backend = p.field_time_series.backend
    interpolating = on_native_grid(backend) 
    ψ_ecco = maybe_interpolate(Val(interpolating), p.field_time_series, i, j, k, p.grid, grid, time)

    ψ = @inbounds fields[i, j, k, p.variable_name]
    μ = stateindex(p.mask, i, j, k, grid, clock.time, loc)
    ω = p.rate

    return ω * μ * (ψ_ecco - ψ)
end

@inline maybe_interpolate(::Val{false}, fts, i, j, k, native_grid, grid, time) = @inbounds fts[i, j, k, time]

@inline function maybe_interpolate(::Val{true}, fts, i, j, k, native_grid, grid, time)
    times = fts.times
    data = fts.data
    time_indexing = fts.time_indexing
    backend = fts.backend
    loc = instantiated_location(fts)
    X = node(i, j, k, grid, loc...)

    # Interpolate field time series data onto the current node and time
    return interpolate(X, time, data, loc, native_grid, times, backend, time_indexing)
end    

"""
    ECCORestoring([arch=CPU(),]
                  variable_name::Symbol;
                  version=ECCO4Monthly(),
                  dates=all_ECCO_dates(version),
                  dates = all_ECCO_dates(version), 
                  time_indices_in_memory = 2, 
                  time_indexing = Cyclical(),
                  mask = 1,
                  rate = 1,
                  grid = nothing,
                  inpainting = NearestNeighborInpainting(prod(size(metadata))),
                  cache_inpainted_data = false)

Create a forcing term that restores to values stored in an ECCO field time series.
The restoring is applied as a forcing on the right hand side of the evolution equations calculated as

```julia
F = mask ⋅ rate ⋅ (ECCO_variable - simulation_variable[i, j, k])
```
where `ECCO_variable` is linearly interpolated in space and time from the ECCO dataset of choice to the 
simulation grid and time.

Arguments
=========

- `arch`: The architecture. Typically `CPU()` or `GPU()`. Default: `CPU()`.

- `variable_name`: The name of the variable to restore. Choices include:
  * `:temperature`,
  * `:salinity`,
  * `:u_velocity`,
  * `:v_velocity`,
  * `:sea_ice_thickness`,
  * `:sea_ice_area_fraction`.

Keyword Arguments
=================

- `version`: The version of the ECCO dataset. Default: `ECCO4Monthly()`.

- `dates`: The dates to use for the ECCO dataset. Default: `all_ECCO_dates(version)`.

- `time_indices_in_memory`: The number of time indices to keep in memory; trade-off between performance
                            and memory footprint.    

- `time_indexing`: The time indexing scheme for the field time series≥

- `mask`: The mask value. Can be a function of `(x, y, z, time)`, an array, or a number.

- `rate`: The restoring rate, i.e., the inverse of the restoring timescale (in s⁻¹).

- `time_indices_in_memory:` how many time instances are loaded in memory; the remaining are loaded lazily.

- `grid`: if not a `nothing`, the ECCO data is directly interpolated on the `grid`.

- `inpainting`: inpainting algorithm, see [`inpaint_mask!`](@ref). Default: `NearestNeighborInpainting(Inf)`.

- `grid`: If `isnothing(grid)`, ECCO data is interpolated on-the-fly to the simulation grid.
          If `!isnothing(grid)`, ECCO data is pre-interpolated to `grid`.
          Default: nothing.

- `cache_inpainted_data`: If `true`, the data is cached to disk after inpainting for later retrieving. 
                          Default: `false`.
                          
It is possible to also pass an `ECCOMetadata` type as the first argument without the need for the 
`variable_name` argument and the `version` and `dates` keyword arguments.
"""
function ECCORestoring(arch::AbstractArchitecture,
                       variable_name::Symbol;
                       version=ECCO4Monthly(),
                       dates=all_ECCO_dates(version),
                       kw...)

    metadata = ECCOMetadata(variable_name, dates, version)
    return ECCORestoring(arch, metadata; kw...)
end

function ECCORestoring(architecture::AbstractArchitecture,
                       metadata::ECCOMetadata;
                       rate,
                       mask = 1,
                       grid = nothing,
                       time_indices_in_memory = 2, # Not more than this if we want to use GPU!
                       time_indexing = Cyclical(),
                       inpainting = NearestNeighborInpainting(Inf),
                       cache_inpainted_data = false)

    fts = ECCO_field_time_series(metadata; 
                                 grid, 
                                 architecture, 
                                 time_indices_in_memory, 
                                 time_indexing, 
                                 inpainting,
                                 cache_inpainted_data)

    # Grab the correct Oceananigans field to restore
    variable_name = metadata.name
    field_name = oceananigans_fieldname[variable_name]

    return ECCORestoring(fts, fts.grid, mask, field_name, rate)
end

# Make sure we can call ECCORestoring with architecture as the first positional argument
ECCORestoring(variable_name::Symbol; kw...) = ECCORestoring(CPU(), variable_name; kw...)
ECCORestoring(metadata::ECCOMetadata; kw...) = ECCORestoring(CPU(), metadata; kw...)

Base.show(io::IO, p::ECCORestoring) = 
    print(io, "Three-dimensional restoring to ECCO data:", '\n',
              "├── restored variable: ", summary(p.variable_name), '\n',
              "├── restoring dataset: ", summary(p.field_time_series.backend.metadata), '\n',
              "├── restoring rate: ", p.rate, '\n',
              "├── mask: ", summary(p.mask), '\n',
              "└── grid: ", summary(p.grid))

regularize_forcing(forcing::ECCORestoring, field, field_name, model_field_names) = forcing
