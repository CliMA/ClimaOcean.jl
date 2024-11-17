module JRA55

using Oceananigans
using Oceananigans.Units

using Oceananigans: location
using Oceananigans.Architectures: arch_array
using Oceananigans.DistributedComputations
using Oceananigans.DistributedComputations: child_architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: λnodes, φnodes, on_architecture
using Oceananigans.Fields: interpolate!
using Oceananigans.OutputReaders: Cyclical, TotallyInMemory, AbstractInMemoryBackend, FlavorOfFTS, time_indices

using ClimaOcean
using ClimaOcean.DataWrangling: download_progress

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

download_jra55_cache::String = ""
function __init__()
    global download_jra55_cache = @get_scratch!("JRA55")
end

struct JRA55Data end
const JRA55PrescribedAtmosphere = PrescribedAtmosphere{<:Any, <:JRA55Data}

# A list of all variables provided in the JRA55 dataset:
JRA55_variable_names = (:river_freshwater_flux,
                        :rain_freshwater_flux,
                        :snow_freshwater_flux,
                        :iceberg_freshwater_flux,
                        :specific_humidity,
                        :sea_level_pressure,
                        :relative_humidity,
                        :downwelling_longwave_radiation,
                        :downwelling_shortwave_radiation,
                        :temperature,
                        :eastward_velocity,
                        :northward_velocity)

filenames = Dict(
    :river_freshwater_flux           => "RYF.friver.1990_1991.nc",   # Freshwater fluxes from rivers
    :rain_freshwater_flux            => "RYF.prra.1990_1991.nc",     # Freshwater flux from rainfall
    :snow_freshwater_flux            => "RYF.prsn.1990_1991.nc",     # Freshwater flux from snowfall
    :iceberg_freshwater_flux         => "RYF.licalvf.1990_1991.nc",  # Freshwater flux from calving icebergs
    :specific_humidity               => "RYF.huss.1990_1991.nc",     # Surface specific humidity
    :sea_level_pressure              => "RYF.psl.1990_1991.nc",      # Sea level pressure
    :relative_humidity               => "RYF.rhuss.1990_1991.nc",    # Surface relative humidity
    :downwelling_longwave_radiation  => "RYF.rlds.1990_1991.nc",     # Downwelling longwave radiation
    :downwelling_shortwave_radiation => "RYF.rsds.1990_1991.nc",     # Downwelling shortwave radiation
    :temperature                     => "RYF.tas.1990_1991.nc",      # Near-surface air temperature
    :eastward_velocity               => "RYF.uas.1990_1991.nc",      # Eastward near-surface wind
    :northward_velocity              => "RYF.vas.1990_1991.nc",      # Northward near-surface wind
)

jra55_short_names = Dict(
    :river_freshwater_flux           => "friver",   # Freshwater fluxes from rivers
    :rain_freshwater_flux            => "prra",     # Freshwater flux from rainfall
    :snow_freshwater_flux            => "prsn",     # Freshwater flux from snowfall
    :iceberg_freshwater_flux         => "licalvf",  # Freshwater flux from calving icebergs
    :specific_humidity               => "huss",     # Surface specific humidity
    :sea_level_pressure              => "psl",      # Sea level pressure
    :relative_humidity               => "rhuss",    # Surface relative humidity
    :downwelling_longwave_radiation  => "rlds",     # Downwelling longwave radiation
    :downwelling_shortwave_radiation => "rsds",     # Downwelling shortwave radiation
    :temperature                     => "tas",      # Near-surface air temperature
    :eastward_velocity               => "uas",      # Eastward near-surface wind
    :northward_velocity              => "vas",      # Northward near-surface wind
)

field_time_series_short_names = Dict(
    :river_freshwater_flux           => "Fri", # Freshwater fluxes from rivers
    :rain_freshwater_flux            => "Fra", # Freshwater flux from rainfall
    :snow_freshwater_flux            => "Fsn", # Freshwater flux from snowfall
    :iceberg_freshwater_flux         => "Fic", # Freshwater flux from calving icebergs
    :specific_humidity               => "qa",  # Surface specific humidity
    :sea_level_pressure              => "pa",  # Sea level pressure
    :relative_humidity               => "rh",  # Surface relative humidity
    :downwelling_longwave_radiation  => "Ql",  # Downwelling longwave radiation
    :downwelling_shortwave_radiation => "Qs",  # Downwelling shortwave radiation
    :temperature                     => "Ta",  # Near-surface air temperature
    :eastward_velocity               => "ua",  # Eastward near-surface wind
    :northward_velocity              => "va",  # Northward near-surface wind
)

urls = Dict(
    :shortwave_radiation => "https://www.dropbox.com/scl/fi/z6fkvmd9oe3ycmaxta131/" *
                            "RYF.rsds.1990_1991.nc?rlkey=r7q6zcbj6a4fxsq0f8th7c4tc&dl=0",

    :river_freshwater_flux => "https://www.dropbox.com/scl/fi/21ggl4p74k4zvbf04nb67/" * 
                              "RYF.friver.1990_1991.nc?rlkey=ny2qcjkk1cfijmwyqxsfm68fz&dl=0",

    :rain_freshwater_flux => "https://www.dropbox.com/scl/fi/5icl1gbd7f5hvyn656kjq/" *
                             "RYF.prra.1990_1991.nc?rlkey=iifyjm4ppwyd8ztcek4dtx0k8&dl=0",

    :snow_freshwater_flux => "https://www.dropbox.com/scl/fi/1r4ajjzb3643z93ads4x4/" *
                             "RYF.prsn.1990_1991.nc?rlkey=auyqpwn060cvy4w01a2yskfah&dl=0",

    :iceberg_freshwater_flux => "https://www.dropbox.com/scl/fi/44nc5y27ohvif7lkvpyv0/" *
                                "RYF.licalvf.1990_1991.nc?rlkey=w7rqu48y2baw1efmgrnmym0jk&dl=0",

    :specific_humidity => "https://www.dropbox.com/scl/fi/66z6ymfr4ghkynizydc29/" *
                          "RYF.huss.1990_1991.nc?rlkey=107yq04aew8lrmfyorj68v4td&dl=0",

    :sea_level_pressure => "https://www.dropbox.com/scl/fi/0fk332027oru1iiseykgp/" *
                           "RYF.psl.1990_1991.nc?rlkey=4xpr9uah741483aukok6d7ctt&dl=0",

    :relative_humidity => "https://www.dropbox.com/scl/fi/1agwsp0lzvntuyf8bm9la/" *
                          "RYF.rhuss.1990_1991.nc?rlkey=8cd0vs7iy1rw58b9pc9t68gtz&dl=0",

    :downwelling_longwave_radiation  => "https://www.dropbox.com/scl/fi/y6r62szkirrivua5nqq61/" *
                                        "RYF.rlds.1990_1991.nc?rlkey=wt9yq3cyrvs2rbowoirf4nkum&dl=0",

    :downwelling_shortwave_radiation => "https://www.dropbox.com/scl/fi/z6fkvmd9oe3ycmaxta131/" *
                                        "RYF.rsds.1990_1991.nc?rlkey=r7q6zcbj6a4fxsq0f8th7c4tc&dl=0",

    :temperature => "https://www.dropbox.com/scl/fi/fpl0npwi476w635g6lke9/" *
                    "RYF.tas.1990_1991.nc?rlkey=0skb9pe6lgbfbiaoybe7m945s&dl=0",

    :eastward_velocity => "https://www.dropbox.com/scl/fi/86wetpqla2x97isp8092g/" *
                          "RYF.uas.1990_1991.nc?rlkey=rcaf18sh1yz0v9g4hjm1249j0&dl=0",

    :northward_velocity => "https://www.dropbox.com/scl/fi/d38sflo9ddljstd5jwgml/" *
                           "RYF.vas.1990_1991.nc?rlkey=f9y3e57kx8xrb40gbstarf0x6&dl=0",
)

compute_bounding_nodes(::Nothing, ::Nothing, LH, hnodes) = nothing
compute_bounding_nodes(bounds, ::Nothing, LH, hnodes) = bounds

function compute_bounding_nodes(x::Number, ::Nothing, LH, hnodes)
    ϵ = convert(typeof(x), 0.001) # arbitrary?
    return (x - ϵ, x + ϵ)
end

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

function compute_bounding_indices(bounds::Tuple, hc)
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
                            time_indices = nothing,
                            latitude = nothing,
                            longitude = nothing,
                            location = nothing,
                            url = nothing,
                            filename = nothing,
                            shortname = nothing,
                            backend = InMemory(),
                            preprocess_chunk_size = 10,
                            preprocess_architecture = CPU())

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

- `latitude`: Guiding latitude bounds for the resulting grid.
              Used to slice the data when loading into memory.
              Default: nothing, which retains the latitude range of the native grid.

- `longitude`: Guiding longitude bounds for the resulting grid.
              Used to slice the data when loading into memory.
              Default: nothing, which retains the longitude range of the native grid.

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
    @root begin
        isfile(filepath) || download(url, filepath)
        isfile(jld2_filepath) && rm(jld2_filepath)
    end

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
    new_data = ds[shortname][i₁:i₂, j₁:j₂, time_indices_in_memory]

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

JRA55PrescribedAtmosphere(time_indices=Colon(); kw...) =
    JRA55PrescribedAtmosphere(CPU(), time_indices; kw...)

JRA55PrescribedAtmosphere(arch::Distributed, time_indices=Colon(); kw...) =
    JRA55PrescribedAtmosphere(child_architecture(arch), time_indices; kw...)

# TODO: allow the user to pass dates
"""
    JRA55PrescribedAtmosphere(architecture::AA, time_indices=Colon();
                                backend = nothing,
                                time_indexing = Cyclical(),
                                reference_height = 10,  # meters
                                include_rivers_and_icebergs = false,
                                other_kw...)

Return a `PrescribedAtmosphere` representing JRA55 reanalysis data.
"""
function JRA55PrescribedAtmosphere(architecture::AA, time_indices=Colon();
                                     backend = nothing,
                                     time_indexing = Cyclical(),
                                     reference_height = 10,  # meters
                                     include_rivers_and_icebergs = false,
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
    pa  = JRA55_field_time_series(:sea_level_pressure;              kw...)
    Fra = JRA55_field_time_series(:rain_freshwater_flux;            kw...)
    Fsn = JRA55_field_time_series(:snow_freshwater_flux;            kw...)
    Ql  = JRA55_field_time_series(:downwelling_longwave_radiation;  kw...)
    Qs  = JRA55_field_time_series(:downwelling_shortwave_radiation; kw...)

    # In JRA55, rivers and icebergs are on a different grid and stored with
    # a different frequency than the rest of the data. Here, we use the
    # `PrescribedAtmosphere.auxiliary_freshwater_flux` feature to represent them.
    if include_rivers_and_icebergs
        Fri = JRA55_field_time_series(:river_freshwater_flux;   kw...)
        Fic = JRA55_field_time_series(:iceberg_freshwater_flux; kw...)
        auxiliary_freshwater_flux = (rivers=Fri, icebergs=Fic)
    else
        auxiliary_freshwater_flux = nothing
    end

    times = ua.times
    freshwater_flux = (rain=Fra, snow=Fsn)
    velocities = (u=ua, v=va)
    tracers = (T=Ta, q=qa)
    pressure = pa
    downwelling_radiation = TwoBandDownwellingRadiation(shortwave=Qs, longwave=Ql)

    FT = eltype(ua)
    boundary_layer_height = convert(FT, 600)
    reference_height = convert(FT, reference_height)
    thermodynamics_parameters = PrescribedAtmosphereThermodynamicsParameters(FT)
    grid = ua.grid
    metadata = JRA55()

    return PrescribedAtmosphere(grid,
                                metadata,
                                velocities,
                                pressure,
                                tracers,
                                freshwater_flux,
                                auxiliary_freshwater_flux,
                                downwelling_radiation,
                                thermodynamics_parameters,
                                times,
                                reference_height,
                                boundary_layer_height)
end

end # module

