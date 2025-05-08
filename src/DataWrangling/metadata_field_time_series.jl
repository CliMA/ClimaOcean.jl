using Oceananigans.Architectures: AbstractArchitecture
using Oceananigans.Grids: AbstractGrid
using Oceananigans.Fields: interpolate!
using Oceananigans.OutputReaders: Cyclical, AbstractInMemoryBackend, FlavorOfFTS, time_indices

import Oceananigans.OutputReaders: new_backend, update_field_time_series!, FieldTimeSeries

@inline instantiate(T::DataType) = T()
@inline instantiate(T) = T

struct DatasetBackend{N, C, I, M} <: AbstractInMemoryBackend{Int}
    start :: Int
    length :: Int
    inpainting :: I
    metadata :: M

    function DatasetBackend{N, C}(start::Int, length::Int, inpainting, metadata) where {N, C}
        M = typeof(metadata)
        I = typeof(inpainting)
        return new{N, C, I, M}(start, length, inpainting, metadata)
    end
end

Adapt.adapt_structure(to, b::DatasetBackend{N, C}) where {N, C} =
    DatasetBackend{N, C}(b.start, b.length, nothing, nothing)

"""
    DatasetBackend(length, metadata;
                   on_native_grid = false,
                   cache_inpainted_data = false,
                   inpainting = NearestNeighborInpainting(Inf))

Represent a FieldTimeSeries backed by the backend that corresponds to the
dataset with `metadata` (e.g., netCDF). Each time instance is stored in an
individual file.
"""
function DatasetBackend(length, metadata;
                        on_native_grid = false,
                        cache_inpainted_data = false,
                        inpainting = NearestNeighborInpainting(Inf))

    return DatasetBackend{on_native_grid, cache_inpainted_data}(1, length, inpainting, metadata)
end

Base.length(backend::DatasetBackend)  = backend.length
Base.summary(backend::DatasetBackend) = string("DatasetBackend(", backend.start, ", ", backend.length, ")")

new_backend(b::DatasetBackend{native, cache_data}, start, length) where {native, cache_data} =
    DatasetBackend{native, cache_data}(start, length, b.inpainting, b.metadata)

on_native_grid(::DatasetBackend{native}) where native = native
cache_inpainted_data(::DatasetBackend{native, cache_data}) where {native, cache_data} = cache_data

const DatasetFieldTimeSeries{N} = FlavorOfFTS{<:Any, <:Any, <:Any, <:Any, <:DatasetBackend{N}} where N

function set!(fts::DatasetFieldTimeSeries)
    @show backend = fts.backend
    @show inpainting = backend.inpainting
    cache_data = cache_inpainted_data(backend)

    @show time_indices(fts)
    for t in time_indices(fts)
        # Set each element of the time-series to the associated file
        metadatum = @inbounds backend.metadata[t]
        set!(fts[t], metadatum; inpainting, cache_inpainted_data=cache_data)
    end

    fill_halo_regions!(fts)

    return nothing
end

"""
    FieldTimeSeries(metadata::Metadata [, arch_or_grid=CPU() ];
                    time_indices_in_memory = 2,
                    time_indexing = Cyclical(),
                    inpainting = nothing,
                    cache_inpainted_data = true)

Create a FieldTimeSeries from a dataset that corresponds to `metadata`.

Arguments
=========

- `metadata`: `Metadata` containing information about the dataset.

- `arch_or_grid`: Either a grid to interpolate the data to, or an `arch`itecture
                  to use for the native grid. Default: CPU().

Keyword Arguments
=================

- `time_indices_in_memory`: The number of time indices to keep in memory. Default: 2.

- `time_indexing`: The time indexing scheme to use. Default: `Cyclical()`.

- `inpainting`: The inpainting algorithm to use for the interpolation.
                The only option is `NearestNeighborInpainting(maxiter)`,
                where an average of the valid surrounding values is used `maxiter` times.

- `cache_inpainted_data`: If `true`, the data is cached to disk after inpainting for later retrieving.
                          Default: `true`.
"""
function FieldTimeSeries(metadata::Metadata, arch::AbstractArchitecture=CPU(); kw...)
    download_dataset(metadata)
    grid = native_grid(metadata, arch)
    return FieldTimeSeries(metadata, grid; kw...)
end

function FieldTimeSeries(metadata::Metadata, grid::AbstractGrid;
                         time_indices_in_memory = 2,
                         time_indexing = Cyclical(),
                         inpainting = default_inpainting(metadata),
                         cache_inpainted_data = true)

    # Make sure all the required individual files are downloaded
    download_dataset(metadata)

    inpainting isa Int && (inpainting = NearestNeighborInpainting(inpainting))
    backend = DatasetBackend(time_indices_in_memory, metadata; on_native_grid, inpainting, cache_inpainted_data)

    times = native_times(metadata)
    loc = LX, LY, LZ = location(metadata)
    boundary_conditions = FieldBoundaryConditions(grid, loc)
    fts = FieldTimeSeries{LX, LY, LZ}(grid, times; backend, time_indexing, boundary_conditions)
    set!(fts)

    return fts
end

function FieldTimeSeries(variable_name::Symbol;
                         dataset, dir,
                         architecture = CPU(),
                         start_date = first_date(dataset, variable_name),
                         end_date = first_date(dataset, variable_name),
                         kw...)

    native_dates = all_dates(dataset, variable_name)
    dates = compute_native_date_range(native_dates, start_date, end_date)
    metadata = Metadata(variable_name, dataset, dates, dir)
    return FieldTimeSeries(metadata, architecture; kw...)
end
