using Oceananigans: location
using Oceananigans.Grids: AbstractGrid, node, on_architecture
using Oceananigans.Fields: interpolate!, interpolate, location, instantiated_location
using Oceananigans.OutputReaders: Cyclical, TotallyInMemory, AbstractInMemoryBackend, FlavorOfFTS, time_indices
using Oceananigans.Utils: Time

using Base
using NCDatasets
using JLD2

using Dates: Second
using ClimaOcean: stateindex
using ClimaOcean.DataWrangling: NearestNeighborInpainting, native_times

import Oceananigans.Fields: set!
import Oceananigans.Forcings: regularize_forcing
import Oceananigans.OutputReaders: new_backend, update_field_time_series!

@inline instantiate(T::DataType) = T()
@inline instantiate(T) = T

struct EN4NetCDFBackend{N, C, I, M} <: AbstractInMemoryBackend{Int}
    start :: Int
    length :: Int
    inpainting :: I
    metadata :: M

    function EN4NetCDFBackend{N, C}(start::Int, length::Int, inpainting, metadata) where {N, C}
        M = typeof(metadata)
        I = typeof(inpainting)
        return new{N, C, I, M}(start, length, inpainting, metadata)
    end
end

Adapt.adapt_structure(to, b::EN4NetCDFBackend{N, C}) where {N, C} = EN4NetCDFBackend{N, C}(b.start, b.length, nothing, nothing)

"""
    EN4NetCDFBackend(length, metadata;
                      on_native_grid = false, 
                      cache_inpainted_data = false,
                      inpainting = NearestNeighborInpainting(Inf))

Represent an EN4 FieldTimeSeries backed by EN4 native netCDF files.
Each time instance is stored in an individual file.
"""
function EN4NetCDFBackend(length, metadata;
                           on_native_grid = false, 
                           cache_inpainted_data = false,
                           inpainting = NearestNeighborInpainting(Inf))

    return EN4NetCDFBackend{on_native_grid, cache_inpainted_data}(1, length, inpainting, metadata)
end

Base.length(backend::EN4NetCDFBackend)  = backend.length
Base.summary(backend::EN4NetCDFBackend) = string("EN4NetCDFBackend(", backend.start, ", ", backend.length, ")")

const EN4FieldTimeSeries{N} = FlavorOfFTS{<:Any, <:Any, <:Any, <:Any, <:EN4NetCDFBackend{N}} where N

new_backend(b::EN4NetCDFBackend{native, cache_data}, start, length) where {native, cache_data} =
    EN4NetCDFBackend{native, cache_data}(start, length, b.inpainting, b.metadata)

on_native_grid(::EN4NetCDFBackend{native}) where native = native
cache_inpainted_data(::EN4NetCDFBackend{native, cache_data}) where {native, cache_data} = cache_data

function set!(fts::EN4FieldTimeSeries) 
    backend = fts.backend
    inpainting = backend.inpainting
    cache_data = cache_inpainted_data(backend)

    for t in time_indices(fts)
        # Set each element of the time-series to the associated file
        metadatum = @inbounds backend.metadata[t] 
        set!(fts[t], metadatum; inpainting, cache_inpainted_data=cache_data)
    end

    fill_halo_regions!(fts)

    return nothing
end

"""
    EN4FieldTimeSeries(metadata::EN4Metadata [, arch_or_grid=CPU() ];
                        time_indices_in_memory = 2,
                        time_indexing = Cyclical(),
                        inpainting = nothing,
                        cache_inpainted_data = true)

Create a field time series object for EN4 data.

Arguments
=========

- `metadata`: `EN4Metadata` containing information about the EN4 dataset.

- `arch_or_grid`: Either a grid to interpolate the EN4 data to, or an `arch`itecture
                  to use for the native EN4 grid. Default: CPU().

Keyword Arguments
=================

- `time_indices_in_memory`: The number of time indices to keep in memory. Default: 2.

- `time_indexing`: The time indexing scheme to use. Default: `Cyclical()`.

- `inpainting`: The inpainting algorithm to use for EN4 interpolation.
                The only option is `NearestNeighborInpainting(maxiter)`, 
                where an average of the valid surrounding values is used `maxiter` times.

- `cache_inpainted_data`: If `true`, the data is cached to disk after inpainting for later retrieving. 
                          Default: `true`.

"""
function EN4FieldTimeSeries(metadata::EN4Metadata, architecture::AbstractArchitecture=CPU(); kw...)
    download_dataset(metadata)
    ftmp = empty_EN4_field(first(metadata); architecture)
    grid = ftmp.grid
    return EN4FieldTimeSeries(metadata, grid; kw...)
end

function EN4FieldTimeSeries(metadata::EN4Metadata, grid::AbstractGrid;
                             time_indices_in_memory = 2,	
                             time_indexing = Cyclical(),
                             inpainting = default_inpainting(metadata),
                             cache_inpainted_data = true)

    # Make sure all the required individual files are downloaded
    download_dataset(metadata)

    inpainting isa Int && (inpainting = NearestNeighborInpainting(inpainting))
    backend = EN4NetCDFBackend(time_indices_in_memory, metadata; on_native_grid, inpainting, cache_inpainted_data)

    times = native_times(metadata)
    loc = LX, LY, LZ = location(metadata)
    boundary_conditions = FieldBoundaryConditions(grid, loc)
    fts = FieldTimeSeries{LX, LY, LZ}(grid, times; backend, time_indexing, boundary_conditions)
    set!(fts)	

    return fts	
end

function EN4FieldTimeSeries(variable_name::Symbol; 
                             dataset = EN4Monthly(),
                             architecture = CPU(),
                             start_date = first_date(dataset, variable_name),
                             end_date = first_date(dataset, variable_name),
                             dir = download_EN4_cache,
                             kw...)

    native_dates = all_dates(dataset, variable_name)
    dates = compute_native_date_range(native_dates, start_date, end_date)                          
    metadata = Metadata(variable_name, dataset, dates, dir)
    return EN4FieldTimeSeries(metadata, architecture; kw...)
end

# Variable names for restorable data
struct Temperature end
struct Salinity end

const oceananigans_fieldnames = Dict(:temperature => Temperature(), 
                                     :salinity    => Salinity())

@inline Base.getindex(fields, i, j, k, ::Temperature) = @inbounds fields.T[i, j, k]
@inline Base.getindex(fields, i, j, k, ::Salinity)    = @inbounds fields.S[i, j, k]

Base.summary(::Temperature) = "temperature"
Base.summary(::Salinity)    = "salinity"

struct EN4Restoring{FTS, G, M, V, N}
    field_time_series :: FTS
    native_grid :: G
    mask :: M
    variable_name :: V
    rate :: N
end

Adapt.adapt_structure(to, p::EN4Restoring) = EN4Restoring(Adapt.adapt(to, p.field_time_series),
                                                            Adapt.adapt(to, p.native_grid),
                                                            Adapt.adapt(to, p.mask),
                                                            Adapt.adapt(to, p.variable_name),
                                                            Adapt.adapt(to, p.rate))

@inline function (p::EN4Restoring)(i, j, k, grid, clock, fields)

    # Figure out all the inputs: time, location, and node
    time = Time(clock.time)
    loc = location(p.field_time_series)

    # Possibly interpolate EN4 data from the EN4 grid to simulation grid.
    # Otherwise, simply extract the pre-interpolated data from p.field_time_series.
    if p.native_grid isa Nothing
        ψ_EN4 = @inbounds p.field_time_series[i, j, k, time]
    else
        ψ_EN4 = interpolate_to_grid(p.field_time_series, i, j, k, p.native_grid, grid, time)
    end
    
    ψ = @inbounds fields[i, j, k, p.variable_name]
    μ = stateindex(p.mask, i, j, k, grid, clock.time, loc)
    r = p.rate

    return r * μ * (ψ_EN4 - ψ)
end

@inline function interpolate_to_grid(fts, i, j, k, native_grid, grid, time)
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
    EN4Restoring(variable_name::Symbol, [ arch_or_grid = CPU(), ];
                  dataset = EN4Monthly(),
                  dates = all_dates(dataset, variable_name),
                  time_indices_in_memory = 2,
                  time_indexing = Cyclical(),
                  mask = 1,
                  rate = 1,
                  dir = download_EN4_cache,
                  inpainting = NearestNeighborInpainting(Inf),
                  cache_inpainted_data = true)

Return a forcing term that restores to values stored in an EN4 field time series.
The restoring is applied as a forcing on the right hand side of the evolution
equations calculated as:

```math
F_ψ = r μ (ψ_{EN4} - ψ)
```

where ``μ`` is the mask, ``r`` is the restoring rate, ``ψ`` is the simulation variable,
and ``ψ_{EN4}`` is the EN4 variable that is linearly interpolated in space and time
from the EN4 dataset of choice to the simulation grid and time.

Arguments
=========

- `variable_name`: The name of the variable to restore. Choices include:
  * `:temperature`,
  * `:salinity`,
  * `:u_velocity`,
  * `:v_velocity`,
  * `:sea_ice_thickness`,
  * `:sea_ice_area_fraction`.

- `arch_or_grid`: Either the architecture of the simulation, or a grid on which the EN4 data
                  is pre-interpolated when loaded. If an `arch`itecture is provided, such as
                  `arch_or_grid = CPU()` or `arch_or_grid = GPU()`, EN4 data are interpolated
                  on-the-fly when the forcing tendency is computed. Default: CPU().

!!! info "Providing `EN4Metadata` instead of `variable_name`"
    Note that `EN4Metadata` may be provided as the first argument instead of `variable_name`.
    In this case the `dataset` and `dates` kwargs (described below) cannot be provided.

Keyword Arguments
=================

- `dataset`: The dataset of the EN4 dataset. Default: `EN4Monthly()`.

- `start_date`: The starting date to use for the EN4 dataset. Default: `first_date(dataset, variable_name)`.

- `end_date`: The ending date to use for the EN4 dataset. Default: `end_date(dataset, variable_name)`.

- `time_indices_in_memory`: The number of time indices to keep in memory. The number is chosen based on 
                            a trade-off between increased performance (more indices in memory) and reduced
                            memory footprint (fewer indices in memory). Default: 2.

- `time_indexing`: The time indexing scheme for the field time series.

- `mask`: The mask value. Can be a function of `(x, y, z, time)`, an array, or a number.

- `rate`: The restoring rate, i.e., the inverse of the restoring timescale (in s⁻¹).

- `dir`: The directory where the native EN4 data is located. If the data does not exist it will
         be automatically downloaded. Default: `download_EN4_cache`.

- `inpainting`: inpainting algorithm, see [`inpaint_mask!`](@ref). Default: `NearestNeighborInpainting(Inf)`.

- `cache_inpainted_data`: If `true`, the data is cached to disk after inpainting for later retrieving. 
                          Default: `true`.
"""
function EN4Restoring(variable_name::Symbol,
                       arch_or_grid = CPU();
                       dataset = EN4Monthly(),
                       start_date = first_date(dataset, variable_name),
                       end_date = last_date(dataset, variable_name),       
                       dir = download_EN4_cache,
                       kw...)

    native_dates = all_dates(dataset, variable_name)
    dates = compute_native_date_range(native_dates, start_date, end_date)                          
    metadata = Metadata(variable_name, dataset, dates, dir)

    return EN4Restoring(metadata, arch_or_grid; kw...)
end

function EN4Restoring(metadata::EN4Metadata,
                       arch_or_grid = CPU();
                       rate,
                       mask = 1,
                       time_indices_in_memory = 2, # Not more than this if we want to use GPU!
                       time_indexing = Cyclical(),
                       inpainting = NearestNeighborInpainting(Inf),
                       cache_inpainted_data = true)

    fts = EN4FieldTimeSeries(metadata, arch_or_grid;
                              time_indices_in_memory, 
                              time_indexing, 
                              inpainting,
                              cache_inpainted_data)

    # Grab the correct Oceananigans field to restore
    variable_name = metadata.name
    field_name = oceananigans_fieldnames[variable_name]

    # If we pass the grid we do not need to interpolate
    # so we can save parameter space by setting the native grid to nothing
    on_native_grid = arch_or_grid isa AbstractArchitecture
    maybe_native_grid = on_native_grid ? fts.grid : nothing

    return EN4Restoring(fts, maybe_native_grid, mask, field_name, rate)
end

function Base.show(io::IO, p::EN4Restoring)
    print(io, "EN4Restoring:", '\n',
              "├── restored variable: ", summary(p.variable_name), '\n',
              "├── restoring dataset: ", summary(p.field_time_series.backend.metadata), '\n',
              "├── restoring rate: ", p.rate, '\n',
              "├── mask: ", summary(p.mask), '\n',
              "└── grid: ", summary(p.native_grid))
end

regularize_forcing(forcing::EN4Restoring, field, field_name, model_field_names) = forcing
