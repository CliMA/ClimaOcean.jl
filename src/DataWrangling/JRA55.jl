module JRA55

using Oceananigans
using Oceananigans.Units
 
using Oceananigans.Architectures: arch_array
using Oceananigans.DistributedComputations
using Oceananigans.DistributedComputations: child_architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: λnodes, φnodes, on_architecture
using Oceananigans.Fields: interpolate!
using Oceananigans.OutputReaders: Cyclical, TotallyInMemory, AbstractInMemoryBackend, FlavorOfFTS, time_indices

using ClimaOcean.OceanSeaIceModels:
    PrescribedAtmosphere,
    TwoBandDownwellingRadiation

using CUDA: @allowscalar

using NCDatasets
using JLD2 
using Dates
using Scratch

import Oceananigans.Fields: set!
import Oceananigans.OutputReaders: new_backend, update_field_time_series!
using Downloads: download

include("jra55_metadata.jl")

download_jra55_cache::String = ""
function __init__()
    global download_jra55_cache = @get_scratch!("JRA55")
end

compute_bounding_nodes(::Nothing, ::Nothing, LH, hnodes) = nothing
compute_bounding_nodes(bounds, ::Nothing, LH, hnodes) = bounds

# TODO: remove the allowscalar
function compute_bounding_nodes(::Nothing, grid, LH, hnodes)
    hg = hnodes(grid, LH())
    h₁ = @allowscalar minimum(hg)
    h₂ = @allowscalar maximum(hg)
    return h₁, h₂
end

function compute_bounding_indices(::Nothing, hc)
    Nh = length(hc)
    return 1, Nh
end

function compute_bounding_indices(bounds, hc)
    h₁, h₂ = bounds
    Nh = length(hc)

    # The following should work. If ᵒ are the extrema of nodes we want to
    # interpolate to, and the following is a sketch of the JRA55 native grid,
    #
    #      1         2         3         4         5
    # |         |         |         |         |         |
    # |    x  ᵒ |    x    |    x    |    x  ᵒ |    x    |
    # |         |         |         |         |         |
    # 1         2         3         4         5         6 
    #
    # then for example, we should find that (iᵢ, i₂) = (1, 5).
    # So we want to reduce the first index by one, and limit them
    # both by the available data. There could be some mismatch due
    # to the use of different coordinate systems (ie whether λ ∈ (0, 360)
    # which we may also need to handle separately.
    i₁ = searchsortedfirst(hc, h₁)
    i₂ = searchsortedfirst(hc, h₂)
    i₁ = max(1, i₁ - 1)
    i₂ = min(Nh, i₂)

    return i₁, i₂
end

infer_longitudinal_topology(::Nothing) = Periodic

function infer_longitudinal_topology(λbounds)
    λ₁, λ₂ = λbounds
    TX = λ₂ - λ₁ ≈ 360 ? Periodic : Bounded
    return TX
end

function compute_bounding_indices(longitude, latitude, grid, LX, LY, λc, φc)
    λbounds = compute_bounding_nodes(longitude, grid, LX, λnodes)
    φbounds = compute_bounding_nodes(latitude, grid, LY, φnodes)

    i₁, i₂ = compute_bounding_indices(λbounds, λc)
    j₁, j₂ = compute_bounding_indices(φbounds, φc)
    TX = infer_longitudinal_topology(λbounds)

    return i₁, i₂, j₁, j₂, TX
end

# Convert dates to range until Oceananigans supports dates natively
function jra55_times(native_times, start_time=native_times[1])

    times = []
    for native_time in native_times
        time = native_time - start_time
        time = Second(time).value
        push!(times, time)
    end

    return times
end

struct JRA55NetCDFBackend <: AbstractInMemoryBackend{Int}
    start :: Int
    length :: Int
end

"""
    JRA55NetCDFBackend(length)

Represents a JRA55 FieldTimeSeries backed by JRA55 native .nc files.
"""
JRA55NetCDFBackend(length) = JRA55NetCDFBackend(1, length)

Base.length(backend::JRA55NetCDFBackend) = backend.length
Base.summary(backend::JRA55NetCDFBackend) = string("JRA55NetCDFBackend(", backend.start, ", ", backend.length, ")")

const JRA55NetCDFFTS = FlavorOfFTS{<:Any, <:Any, <:Any, <:Any, <:JRA55NetCDFBackend}

function set!(fts::JRA55NetCDFFTS, path::String=fts.path, name::String=fts.name) 

    ds = Dataset(path)

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
    ti = collect(ti)
    data = ds[name][i₁:i₂, j₁:j₂, ti]
    close(ds)

    copyto!(interior(fts, :, :, 1, :), data)
    fill_halo_regions!(fts)

    return nothing
end

new_backend(::JRA55NetCDFBackend, start, length) = JRA55NetCDFBackend(start, length)

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
function JRA55_field_time_series(variable_name;
                                 architecture = CPU(),
                                 grid = nothing,
                                 location = nothing,
                                 url = nothing,
                                 dir = download_jra55_cache,
                                 filename = nothing,
                                 shortname = nothing,
                                 latitude = nothing,
                                 longitude = nothing,
                                 backend = InMemory(),
                                 time_indexing = Cyclical(),
                                 preprocess_chunk_size = 10,
                                 preprocess_architecture = CPU(),
                                 time_indices = nothing)

    # OnDisk backends do not support time interpolation!
    # Disallow OnDisk for JRA55 dataset loading 
    if backend isa OnDisk 
        msg = string("We cannot load the JRA55 dataset with an `OnDisk` backend")
        throw(ArgumentError(msg))
    end

    if isnothing(filename) && !(variable_name ∈ JRA55_variable_names)
        variable_strs = Tuple("  - :$name \n" for name in JRA55_variable_names)
        variables_msg = prod(variable_strs)

        msg = string("The variable :$variable_name is not provided by the JRA55-do dataset!", '\n',
                     "The variables provided by the JRA55-do dataset are:", '\n',
                     variables_msg)

        throw(ArgumentError(msg))
    end

    filepath = isnothing(filename) ? joinpath(dir, filenames[variable_name]) : joinpath(dir, filename)

    if !isnothing(filename) && !isfile(filepath) && isnothing(url)
        throw(ArgumentError("A filename was provided without a url, but the file does not exist.\n \
                            If intended, please provide both the filename and url that should be used \n \
                            to download the new file."))
    end

    isnothing(filename)  && (filename  = filenames[variable_name])
    isnothing(shortname) && (shortname = jra55_short_names[variable_name])
    isnothing(url)       && (url       = urls[variable_name])

    # Record some important user decisions
    totally_in_memory = backend isa TotallyInMemory
    on_native_grid = isnothing(grid)
    !on_native_grid && backend isa JRA55NetCDFBackend && error("Can't use custom grid with JRA55NetCDFBackend.")

    jld2_filepath = joinpath(dir, string("JRA55_repeat_year_", variable_name, ".jld2"))
    fts_name = field_time_series_short_names[variable_name]

    # Note, we don't re-use existing jld2 files.
    isfile(filepath) || download(url, filepath)
    isfile(jld2_filepath) && rm(jld2_filepath)

    # Determine default time indices
    if totally_in_memory
        # In this case, the whole time series is in memory.
        # Either the time series is short, or we are doing a limited-area
        # simulation, like in a single column. So, we conservatively
        # set a default `time_indices = 1:2`.
        isnothing(time_indices) && (time_indices = 1:2)
        time_indices_in_memory  = time_indices
        native_fts_architecture = architecture
    else
        # In this case, part or all of the time series will be stored in a file.
        # Note: if the user has provided a grid, we will have to preprocess the
        # .nc JRA55 data into a .jld2 file. In this case, `time_indices` refers
        # to the time_indices that we will preprocess;
        # by default we choose all of them. The architecture is only the
        # architecture used for preprocessing, which typically will be CPU()
        # even if we would like the final FieldTimeSeries on the GPU.
        isnothing(time_indices) && (time_indices = :)

        if backend isa JRA55NetCDFBackend
            time_indices_in_memory = 1:length(backend)
            native_fts_architecture = architecture
        else # then `time_indices_in_memory` refers to preprocessing
            maximum_index = min(preprocess_chunk_size, length(time_indices))
            time_indices_in_memory = 1:maximum_index
            native_fts_architecture = preprocess_architecture
        end
    end

    # Set a default location.
    if isnothing(location)
        LX = LY = Center
    else
        LX, LY = location
    end

    ds = Dataset(filepath)

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

    JRA55_native_grid = LatitudeLongitudeGrid(native_fts_architecture, Float32;
                                              halo = H,
                                              size = N,
                                              longitude = λr,
                                              latitude = φr,
                                              topology = (TX, Bounded, Flat))

    boundary_conditions = FieldBoundaryConditions(JRA55_native_grid, (Center, Center, Nothing))
    times = jra55_times(native_times)
    
    if backend isa JRA55NetCDFBackend
        fts = FieldTimeSeries{Center, Center, Nothing}(JRA55_native_grid, times;
                                                       backend,
                                                       time_indexing,
                                                       boundary_conditions,
                                                       path = filepath,
                                                       name = shortname)

        # Fill the data in a GPU-friendly manner
        copyto!(interior(fts, :, :, 1, :), data)
        fill_halo_regions!(fts)

        return fts
    else
        # Make times into an array for later preprocessing
        if !totally_in_memory
            times = collect(times)
        end

        native_fts = FieldTimeSeries{Center, Center, Nothing}(JRA55_native_grid, times;
                                                              time_indexing,
                                                              boundary_conditions)

        # Fill the data in a GPU-friendly manner
        copyto!(interior(native_fts, :, :, 1, :), data)
        fill_halo_regions!(native_fts)

        if on_native_grid && totally_in_memory
            return native_fts

        elseif totally_in_memory # but not on the native grid!
            boundary_conditions = FieldBoundaryConditions(grid, (LX, LY, Nothing))
            fts = FieldTimeSeries{LX, LY, Nothing}(grid, times; time_indexing, boundary_conditions)
            interpolate!(fts, native_fts)
            return fts
        end
    end

    @info "Pre-processing JRA55 $variable_name data into a JLD2 file..."

    preprocessing_grid = on_native_grid ? JRA55_native_grid : grid

    # Re-open the dataset!
    ds = Dataset(filepath)
    all_datetimes = ds["time"][time_indices]
    all_Nt = length(all_datetimes)

    all_times = jra55_times(all_datetimes)

    on_disk_fts = FieldTimeSeries{LX, LY, Nothing}(preprocessing_grid, all_times;
                                                   boundary_conditions,
                                                   backend = OnDisk(),
                                                   path = jld2_filepath,
                                                   name = fts_name)

    # Save data to disk, one field at a time
    start_clock = time_ns()
    n = 1 # on disk
    m = 0 # in memory

    times_in_memory = all_times[time_indices_in_memory]

    fts = FieldTimeSeries{LX, LY, Nothing}(preprocessing_grid, times_in_memory;
                                           boundary_conditions,
                                           backend = InMemory(),
                                           path = jld2_filepath,
                                           name = fts_name)

    # Re-compute data
    new_data  = ds[shortname][i₁:i₂, j₁:j₂, time_indices_in_memory]

    if !on_native_grid
        copyto!(interior(native_fts, :, :, 1, :), new_data[:, :, :])
        fill_halo_regions!(native_fts)    
        interpolate!(fts, native_fts)
    else
        copyto!(interior(fts, :, :, 1, :), new_data[:, :, :])
    end
                     
    while n <= all_Nt
        print("        ... processing time index $n of $all_Nt \r")

        if time_indices_in_memory isa Colon || n ∈ time_indices_in_memory
            m += 1
        else # load new data
            # Update time_indices
            time_indices_in_memory = time_indices_in_memory .+ preprocess_chunk_size
            n₁ = first(time_indices_in_memory)

            # Clip time_indices if they extend past the end of the dataset
            if last(time_indices_in_memory) > all_Nt
                time_indices_in_memory = UnitRange(n₁, all_Nt)
            end

            # Re-compute times
            new_times = jra55_times(all_times[time_indices_in_memory], all_times[n₁])
            native_fts.times = new_times

            # Re-compute data
            new_data  = ds[shortname][i₁:i₂, j₁:j₂, time_indices_in_memory]
            fts.times = new_times

            if !on_native_grid
                copyto!(interior(native_fts, :, :, 1, :), new_data[:, :, :])
                fill_halo_regions!(native_fts)    
                interpolate!(fts, native_fts)
            else
                copyto!(interior(fts, :, :, 1, :), new_data[:, :, :])
            end

            m = 1 # reset
        end

        set!(on_disk_fts, fts[m], n, fts.times[m])

        n += 1
    end

    elapsed = 1e-9 * (time_ns() - start_clock)
    elapsed_str = prettytime(elapsed)
    @info "    ... done ($elapsed_str)" * repeat(" ", 20)

    close(ds)

    user_fts = FieldTimeSeries(jld2_filepath, fts_name; architecture, backend, time_indexing)
    fill_halo_regions!(user_fts)

    return user_fts
end

const AA = Oceananigans.Architectures.AbstractArchitecture

JRA55_prescribed_atmosphere(time_indices=Colon(); kw...) =
    JRA55_prescribed_atmosphere(CPU(), time_indices; kw...)

JRA55_prescribed_atmosphere(arch::Distributed, time_indices=Colon(); kw...) =
    JRA55_prescribed_atmosphere(child_architecture(arch), time_indices; kw...)

# TODO: allow the user to pass dates
function JRA55_prescribed_atmosphere(architecture::AA, time_indices=Colon();
                                     backend = nothing,
                                     time_indexing = Cyclical(),
                                     reference_height = 10,  # meters
                                     include_rivers_and_icebergs = true, # rivers and icebergs are not needed in single column simulations
                                     other_kw...)

    if isnothing(backend) # apply a default
        Ni = try
            length(time_indices)
        catch
            Inf
        end

        # Manufacture a default for the number of fields to keep InMemory
        Nf = min(24, Ni)
        backend = JRA55NetCDFBackend(Nf)
    end

    kw = (; time_indices, time_indexing, backend, architecture)
    kw = merge(kw, other_kw) 

    ua  = JRA55_field_time_series(:eastward_velocity;               kw...)
    va  = JRA55_field_time_series(:northward_velocity;              kw...)
    Ta  = JRA55_field_time_series(:temperature;                     kw...)
    qa  = JRA55_field_time_series(:specific_humidity;               kw...)
    ra  = JRA55_field_time_series(:relative_humidity;               kw...)
    pa  = JRA55_field_time_series(:sea_level_pressure;              kw...)
    Fra = JRA55_field_time_series(:rain_freshwater_flux;            kw...)
    Fsn = JRA55_field_time_series(:snow_freshwater_flux;            kw...)
    Ql  = JRA55_field_time_series(:downwelling_longwave_radiation;  kw...)
    Qs  = JRA55_field_time_series(:downwelling_shortwave_radiation; kw...)

    freshwater_flux = (rain = Fra,
                       snow = Fsn)

    # Remember that rivers and icebergs are on a different grid and have
    # a different frequency than the rest of the JRA55 data
    if include_rivers_and_icebergs
        Fri = JRA55_field_time_series(:river_freshwater_flux;   kw...)
        Fic = JRA55_field_time_series(:iceberg_freshwater_flux; kw...)
        runoff_flux = (rivers   = Fri,
                       icebergs = Fic)
    else
        runoff_flux = nothing
    end

    times = ua.times

    velocities = (u = ua,
                  v = va)

    tracers = (T = Ta,
               q = qa,
               r = ra)

                       
    pressure = pa

    downwelling_radiation = TwoBandDownwellingRadiation(shortwave=Qs, longwave=Ql)

    FT = eltype(ua)
    reference_height = convert(FT, reference_height)

    atmosphere = PrescribedAtmosphere(times, FT;
                                      velocities,
                                      freshwater_flux,
                                      runoff_flux,
                                      tracers,
                                      downwelling_radiation,
                                      reference_height,
                                      pressure)

    return atmosphere
end

end # module

