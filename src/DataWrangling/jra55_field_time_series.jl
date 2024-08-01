using Oceananigans.Units
using Oceananigans.Grids: node, on_architecture
using Oceananigans.Fields: interpolate!, interpolate, location, instantiated_location
using Oceananigans.OutputReaders: Cyclical, TotallyInMemory, AbstractInMemoryBackend, FlavorOfFTS, time_indices
using Oceananigans.Utils: Time

using CUDA: @allowscalar
using Base

using NCDatasets
using JLD2 
using Dates

using ClimaOcean: stateindex
using ClimaOcean.DataWrangling: native_times

import Oceananigans.Fields: set!
import Oceananigans.OutputReaders: new_backend, update_field_time_series!

@inline instantiate(T::DataType) = T()
@inline instantiate(T) = T

struct JRA55NetCDFBackend{N} <: AbstractInMemoryBackend{Int}
    start :: Int
    length :: Int

    JRA55NetCDFBackend{N}(start::Int, length::Int) where N = new{N}(start, length)
end

"""
    JRA55NetCDFBackend(length)

Represents an JRA55 FieldTimeSeries backed by JRA55 native .nc files.
Each time instance is stored in an individual file.
"""
JRA55NetCDFBackend(length; on_native_grid = false) = JRA55NetCDFBackend{on_native_grid}(1, length)

Base.length(backend::JRA55NetCDFBackend)  = backend.length
Base.summary(backend::JRA55NetCDFBackend) = string("JRA55NetCDFBackend(", backend.start, ", ", backend.length, ")")

const JRA55NetCDFFTS{N} = FlavorOfFTS{<:Any, <:Any, <:Any, <:Any, <: JRA55NetCDFBackend{N}} where N

new_backend(::JRA55NetCDFBackend{N}, start, length) where N = JRA55ONetCDFBackend{N}(start, length)
on_native_grid(::JRA55NetCDFBackend{N}) where N = N

function set!(fts::JRA55NetCDFFTS, path::JRA55Metadata=fts.path, name::String=fts.name) 
 
    # Do different things based on the Backend...
    filename = metadata_filename(path)
    ds = filename isa String ? Dataset(filename) : Dataset(filename[1])

    # Note that each file should have the variables
    #   - ds["time"]:     time coordinate 
    #   - ds["lon"]:      longitude at the location of the variable
    #   - ds["lat"]:      latitude at the location of the variable
    #   - ds["lon_bnds"]: bounding longitudes between which variables are averaged
    #   - ds["lat_bnds"]: bounding latitudes between which variables are averaged
    #   - ds[shortname]:  the variable data

    # Nodes at the variable location
    λc = ds["lon"][:]
    φc = ds["lat"][:]
    LX, LY, LZ = location(fts)
    i₁, i₂, j₁, j₂, TX = compute_bounding_indices(nothing, nothing, fts.grid, LX, LY, λc, φc)

    ti = time_indices(fts)

    if on_native_grid(fts.backend)
        ti = collect(ti)
        data = ds[name][i₁:i₂, j₁:j₂, ti]
        copyto!(interior(fts, :, :, 1, :), data)
    end

    close(ds)
    fill_halo_regions!(fts)
    
    return nothing
end

"""
    JRA55_field_time_series(metadata::ECCOMetadata;
                           architecture = CPU(),
                           time_indices_in_memory = 2,
                           time_indexing = Cyclical(),
                           grid = nothing)

Create a field time series object for JRA55 data.

Arguments:
===========
- metadata: An JRA55Metadata object containing information about the JRA55 dataset.

Keyword Arguments:
=====================
- architecture: The architecture to use for computations (default: CPU()).
- time_indices_in_memory: The number of time indices to keep in memory (default: 2).
- time_indexing: The time indexing scheme to use (default: Cyclical()).
- grid: if not a `nothing`, the ECCO data is directly interpolated on the `grid`,
"""
function JRA55_field_time_series(metadata::JRA55Metadata;	
                                 architecture = CPU(),	
                                 time_indices_in_memory = 2,
                                 time_indexing = Cyclical(),
                                 grid = nothing,
                                 latitude = nothing,
                                 longitude = nothing)	

    backend = JRA55NetCDFBackend(time_indices_in_memory; on_native_grid = isnothing(grid))

    # Making sure all the required individual files are downloaded
    download_dataset!(metadata)

    filename = metadata_filename(metadata)
    ds = Dataset(filename)

    # Note that each file should have the variables
    #   - ds["time"]:     time coordinate 
    #   - ds["lon"]:      longitude at the location of the variable
    #   - ds["lat"]:      latitude at the location of the variable
    #   - ds["lon_bnds"]: bounding longitudes between which variables are averaged
    #   - ds["lat_bnds"]: bounding latitudes between which variables are averaged
    #   - ds[shortname]: the variable data

    # Nodes at the variable location
    λc = ds["lon"][:]
    φc = ds["lat"][:]

    # Interfaces for the "native" JRA55 grid
    λn = ds["lon_bnds"][1, :]
    φn = ds["lat_bnds"][1, :]

    # The .nc coordinates lon_bnds and lat_bnds do not include
    # the last interface, so we push them here.
    push!(φn, 90)
    push!(λn, λn[1] + 360)

    # TODO: support loading just part of the JRA55 data.
    # Probably with arguments that take latitude, longitude bounds.
    i₁, i₂, j₁, j₂, TX = compute_bounding_indices(longitude, latitude, grid, LX, LY, λc, φc)

    native_times = ds["time"][time_indices]
    data = ds[shortname][i₁:i₂, j₁:j₂, time_indices_in_memory]
    λr = λn[i₁:i₂+1]
    φr = φn[j₁:j₂+1]
    Nrx, Nry, Nt = size(data)
    close(ds)

    N = (Nrx, Nry)
    H = min.(N, (3, 3))

    JRA55_native_grid = LatitudeLongitudeGrid(architecture, Float32;
                                              halo = H,
                                              size = N,
                                              longitude = λr,
                                              latitude = φr,
                                              topology = (TX, Bounded, Flat))

    boundary_conditions = FieldBoundaryConditions(JRA55_native_grid, (Center, Center, Nothing))

    location  = field_location(metadata)
    shortname = short_name(metadata)

    boundary_conditions = FieldBoundaryConditions(JRA55_native_grid, location)
    times = native_times(metadata)

    fts_grid = isnothing(grid) ? JRA55_native_grid : grid

    path = metadata
    
    fts = FieldTimeSeries{location...}(fts_grid, times;	
                                       backend,	
                                       time_indexing,	
                                       boundary_conditions,	
                                       path,	
                                       name = metadata.name)

    # Let's set the data	
    set!(fts)	

    return fts	
end

JRA55_field_time_series(variable_name::Symbol, version=JRA55RepeatYear(); kw...) = 
    JRA55_field_time_series(Metadata(variable_name, all_dates(version), version); kw...)

"""
    JRA55_field_time_series(variable_name;
                            architecture = CPU(),
                            location = nothing,
                            url = nothing,
                            filename = nothing,
                            shortname = nothing,
                            backend = InMemory(),
                            preprocess_chunk_size = 10,
                            preprocess_architecture = CPU(),
                            time_indices = nothing)

Return a `FieldTimeSeries` containing atmospheric reanalysis data for `variable_name`,
which describes one of the variables in the "repeat year forcing" dataset derived
from the Japanese 55-year atmospheric reanalysis for driving ocean-sea-ice models (JRA55-do).
For more information about the derivation of the repeat year forcing dataset, see

"Stewart et al., JRA55-do-based repeat year forcing datasets for driving ocean–sea-ice models",
Ocean Modelling, 2020, https://doi.org/10.1016/j.ocemod.2019.101557.

The `variable_name`s (and their `shortname`s used in NetCDF files)
available from the JRA55-do are:

    - `:river_freshwater_flux`              ("friver")
    - `:rain_freshwater_flux`               ("prra")
    - `:snow_freshwater_flux`               ("prsn")
    - `:iceberg_freshwater_flux`            ("licalvf")
    - `:specific_humidity`                  ("huss")
    - `:sea_level_pressure`                 ("psl")
    - `:relative_humidity`                  ("rhuss")
    - `:downwelling_longwave_radiation`     ("rlds")
    - `:downwelling_shortwave_radiation`    ("rsds")
    - `:temperature`                        ("ras")
    - `:eastward_velocity`                  ("uas")
    - `:northward_velocity`                 ("vas")

Keyword arguments
=================

    - `architecture`: Architecture for the `FieldTimeSeries`.
                      Default: CPU()

    - `time_indices`: Indices of the timeseries to extract from file. 
                      For example, `time_indices=1:3` returns a 
                      `FieldTimeSeries` with the first three time snapshots
                      of `variable_name`.

    - `url`: The url accessed to download the data for `variable_name`.
             Default: `ClimaOcean.JRA55.urls[variable_name]`.

    - `filename`: The name of the downloaded file.
                  Default: `ClimaOcean.JRA55.filenames[variable_name]`.

    - `shortname`: The "short name" of `variable_name` inside its NetCDF file.
                    Default: `ClimaOcean.JRA55.jra55_short_names[variable_name]`.

    - `interpolated_file`: file holding an Oceananigans compatible version of the JRA55 data.
                            If it does not exist it will be generated.

    - `time_chunks_in_memory`: number of fields held in memory. If `nothing` the whole timeseries is 
                               loaded (not recommended).
"""
