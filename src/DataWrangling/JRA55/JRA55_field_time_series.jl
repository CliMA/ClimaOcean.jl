using ClimaOcean.DataWrangling: all_dates, native_times
using ClimaOcean.DataWrangling: compute_native_date_range
using Oceananigans.Grids: AbstractGrid
using Oceananigans.OutputReaders: PartlyInMemory
using Adapt

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

struct JRA55NetCDFBackend{M} <: AbstractInMemoryBackend{Int}
    start :: Int
    length :: Int
    metadata :: M
end

Adapt.adapt_structure(to, b::JRA55NetCDFBackend) = JRA55NetCDFBackend(b.start, b.length, nothing)

"""
    JRA55NetCDFBackend(length)

Represents a JRA55 FieldTimeSeries backed by JRA55 native netCDF files.
"""
JRA55NetCDFBackend(length, metadata::Metadata) = JRA55NetCDFBackend(1, length, metadata)
JRA55NetCDFBackend(start::Integer, length::Integer) = JRA55NetCDFBackend(start, length, nothing)

# Metadata - agnostic constructor
JRA55NetCDFBackend(length) = JRA55NetCDFBackend(1, length, nothing)

Base.length(backend::JRA55NetCDFBackend) = backend.length
Base.summary(backend::JRA55NetCDFBackend) = string("JRA55NetCDFBackend(", backend.start, ", ", backend.length, ")")

const JRA55NetCDFFTS              = FlavorOfFTS{<:Any, <:Any, <:Any, <:Any, <:JRA55NetCDFBackend}
const JRA55NetCDFFTSRepeatYear    = FlavorOfFTS{<:Any, <:Any, <:Any, <:Any, <:JRA55NetCDFBackend{<:Metadata{<:RepeatYearJRA55}}}
const JRA55NetCDFFTSMultipleYears = FlavorOfFTS{<:Any, <:Any, <:Any, <:Any, <:JRA55NetCDFBackend{<:Metadata{<:MultiYearJRA55}}}

# Note that each file should have the variables
#   - ds["time"]:     time coordinate
#   - ds["lon"]:      longitude at the location of the variable
#   - ds["lat"]:      latitude at the location of the variable
#   - ds["lon_bnds"]: bounding longitudes between which variables are averaged
#   - ds["lat_bnds"]: bounding latitudes between which variables are averaged
#   - ds[shortname]:  the variable data

# Simple case, only one file per variable, no need to deal with multiple files
function set!(fts::JRA55NetCDFFTSRepeatYear, backend=fts.backend)

    metadata = backend.metadata

    filename = metadata_filename(metadata)
    path = joinpath(metadata.dir, filename)
    ds = Dataset(path)

    # Nodes at the variable location

    λc = ds["lon"][:]
    φc = ds["lat"][:]
    LX, LY, LZ = location(fts)
    i₁, i₂, j₁, j₂, TX = compute_bounding_indices(nothing, nothing, fts.grid, LX, LY, λc, φc)

    nn   = time_indices(fts)
    nn   = collect(nn)
    name = dataset_variable_name(fts.backend.metadata)

    if issorted(nn)
        data = ds[name][i₁:i₂, j₁:j₂, nn]
    else
        # The time indices may be cycling past 1; eg ti = [6, 7, 8, 1].
        # However, DiskArrays does not seem to support loading data with unsorted
        # indices. So to handle this, we load the data in chunks, where each chunk's
        # indices are sorted, and then glue the data together.
        m = findfirst(n -> n == 1, nn)
        n1 = nn[1:m-1]
        n2 = nn[m:end]

        data1 = ds[name][i₁:i₂, j₁:j₂, n1]
        data2 = ds[name][i₁:i₂, j₁:j₂, n2]
        data = cat(data1, data2, dims=3)
    end

    close(ds)

    copyto!(interior(fts, :, :, 1, :), data)
    fill_halo_regions!(fts)

    return nothing
end

# Tricky case: multiple files per variable -- one file per year --
# we need to infer the file name from the metadata and split the data loading
function set!(fts::JRA55NetCDFFTSMultipleYears, backend=fts.backend)

    metadata = backend.metadata

    filename   = metadata_filename(metadata)
    filename   = unique(filename)
    name       = dataset_variable_name(metadata)
    start_date = first_date(metadata.dataset, metadata.name)

    for file in filename

        path = joinpath(metadata.dir, file)
        ds = Dataset(path)

        # This can be simplified once we start supporting a
        # datetime `Clock` in Oceananigans
        file_dates = ds["time"][:]
        file_indices = 1:length(file_dates)
        file_times = zeros(length(file_dates))
        for (t, date) in enumerate(file_dates)
            delta = date - start_date
            delta = Second(delta).value
            file_times[t] = delta
        end

        ftsn = time_indices(fts)
        ftsn = collect(ftsn)

        # Intersect the time indices with the file times
        nn   = findall(n -> file_times[n] ∈ fts.times[ftsn], file_indices)
        ftsn = findall(n -> fts.times[n] ∈ file_times[nn], ftsn)

        if !isempty(nn)
            # Nodes at the variable location
            λc = ds["lon"][:]
            φc = ds["lat"][:]
            LX, LY, LZ = location(fts)
            i₁, i₂, j₁, j₂, TX = compute_bounding_indices(nothing, nothing, fts.grid, LX, LY, λc, φc)


            if issorted(nn)
                data = ds[name][i₁:i₂, j₁:j₂, nn]
            else
                # The time indices may be cycling past 1; eg ti = [6, 7, 8, 1].
                # However, DiskArrays does not seem to support loading data with unsorted
                # indices. So to handle this, we load the data in chunks, where each chunk's
                # indices are sorted, and then glue the data together.
                m = findfirst(n -> n == 1, nn)
                n1 = nn[1:m-1]
                n2 = nn[m:end]

                data1 = ds[name][i₁:i₂, j₁:j₂, n1]
                data2 = ds[name][i₁:i₂, j₁:j₂, n2]
                data = cat(data1, data2, dims=3)
            end

            close(ds)

            # We need to set the time index for each file
            # Find start index corresponding to the underlying data
            for n in 1:length(nn)
                copyto!(interior(fts, :, :, 1, ftsn[n]), data[:, :, n])
            end
        end
    end

    fill_halo_regions!(fts)

    return nothing
end

new_backend(b::JRA55NetCDFBackend, start, length) = JRA55NetCDFBackend(start, length, b.metadata)

"""
    JRA55FieldTimeSeries(variable_name, architecture=CPU(), FT=Float32;
                         dataset = RepeatYearJRA55(),
                         dates = all_JRA55_dates(version),
                         latitude = nothing,
                         longitude = nothing,
                         dir = download_JRA55_cache,
                         backend = InMemory(),
                         time_indexing = Cyclical())

Return a `FieldTimeSeries` containing atmospheric reanalysis data for `variable_name`,
which describes one of the variables from the Japanese 55-year atmospheric reanalysis
for driving ocean-sea ice models (JRA55-do). The JRA55-do dataset is described by [tsujino2018jra](@citet).

The `variable_name`s (and their `shortname`s used in the netCDF files) available from the JRA55-do are:
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

- `architecture`: Architecture for the `FieldTimeSeries`. Default: CPU()

- `dataset`: The data dataset; supported datasets are: `RepeatYearJRA55()` and `MultiYearJRA55()`.
             `MultiYearJRA55()` refers to the full length of the JRA55-do dataset; `RepeatYearJRA55()`
             refers to the "repeat-year forcing" dataset derived from JRA55-do. Default: `RepeatYearJRA55()`.

  !!! info "Repeat-year forcing"

      For more information about the derivation of the repeat-year forcing dataset, see [stewart2020jra55](@citet).

  The repeat year in `RepeatYearJRA55()` corresponds to May 1st, 1990 - April 30th, 1991. However, the
  returned dataset has dates that range from January 1st to December 31st. This implies
  that the first 4 months of the `JRA55RepeatYear()` dataset correspond to year 1991 from the JRA55
  reanalysis and the rest 8 months from 1990.

- `start_date`: The starting date to use for the dataset. Default: `first_date(dataset, variable_name)`.

- `end_date`: The ending date to use for the dataset. Default: `end_date(dataset, variable_name)`.

- `dir`: The directory of the data file. Default: `ClimaOcean.JRA55.download_JRA55_cache`.

- `time_indexing`: The time indexing scheme for the field time series. Default: `Cyclical()`.

- `latitude`: Guiding latitude bounds for the resulting grid.
              Used to slice the data when loading into memory.
              Default: nothing, which retains the latitude range of the native grid.

- `longitude`: Guiding longitude bounds for the resulting grid.
               Used to slice the data when loading into memory.
               Default: nothing, which retains the longitude range of the native grid.

- `backend`: Backend for the `FieldTimeSeries`. The two options are:
  * `InMemory()`: the whole time series is loaded into memory.
  * `JRA55NetCDFBackend(total_time_instances_in_memory)`: only a subset of the time series
                                                          is loaded into memory. Default: `InMemory()`.

References
==========

- Tsujino et al. (2018). JRA-55 based surface dataset for driving ocean-sea-ice models (JRA55-do), _Ocean Modelling_, **130(1)**, 79-139.

- Stewart et al. (2020). JRA55-do-based repeat year forcing datasets for driving ocean–sea-ice models, _Ocean Modelling_, **147**, 101557.
"""
function JRA55FieldTimeSeries(variable_name::Symbol, architecture=CPU(), FT=Float32;
                              dataset = RepeatYearJRA55(),
                              start_date = first_date(dataset, variable_name),
                              end_date = last_date(dataset, variable_name),
                              dir = download_JRA55_cache,
                              kw...)

    native_dates = all_dates(dataset, variable_name)
    dates = compute_native_date_range(native_dates, start_date, end_date)
    metadata = Metadata(variable_name; dataset, dates, dir)

    return JRA55FieldTimeSeries(metadata, architecture, FT; kw...)
end

function JRA55FieldTimeSeries(metadata::JRA55Metadata, architecture=CPU(), FT=Float32;
                              latitude = nothing,
                              longitude = nothing,
                              backend = InMemory(),
                              time_indexing = Cyclical())

    # Cannot use `TotallyInMemory` backend with MultiYearJRA55 dataset
    if metadata.dataset isa MultiYearJRA55 && backend isa TotallyInMemory
        msg = string("The `InMemory` backend is not supported for the MultiYearJRA55 dataset.")
        throw(ArgumentError(msg))
    end

    # First thing: we download the dataset!
    download_dataset(metadata)

    # Regularize the backend in case of `JRA55NetCDFBackend`
    if backend isa JRA55NetCDFBackend
        if backend.metadata isa Nothing
            backend = JRA55NetCDFBackend(backend.length, metadata)
        end

        if backend.length > length(metadata)
            backend = JRA55NetCDFBackend(backend.start, length(metadata), metadata)
        end
    end

    # Unpack metadata details
    dataset = metadata.dataset
    name    = metadata.name
    time_indices = JRA55_time_indices(dataset, metadata.dates, name)

    # Change the metadata to reflect the actual time indices
    dates    = all_dates(dataset, name)[time_indices]
    metadata = Metadata(metadata.name; dataset=metadata.dataset, dates, dir=metadata.dir)

    shortname = dataset_variable_name(metadata)
    variable_name = metadata.name

    filepath = metadata_path(metadata) # Might be multiple paths!!!
    filepath = filepath isa AbstractArray ? first(filepath) : filepath

    # OnDisk backends do not support time interpolation!
    # Disallow OnDisk for JRA55 dataset loading
    if ((backend isa InMemory) && !isnothing(backend.length)) || backend isa OnDisk
        msg = string("We cannot load the JRA55 dataset with a $(backend) backend. Use `InMemory()` or `JRA55NetCDFBackend(N)` instead.")
        throw(ArgumentError(msg))
    end

    if !(variable_name ∈ JRA55_variable_names)
        variable_strs = Tuple("  - :$name \n" for name in JRA55_variable_names)
        variables_msg = prod(variable_strs)

        msg = string("The variable :$variable_name is not provided by the JRA55-do dataset!", '\n',
                     "The variables provided by the JRA55-do dataset are:", '\n',
                     variables_msg)

        throw(ArgumentError(msg))
    end

    # Record some important user decisions
    totally_in_memory = backend isa TotallyInMemory

    # Determine default time indices
    if totally_in_memory
        # In this case, the whole time series is in memory.
        # Either the time series is short, or we are doing a limited-area
        # simulation, like in a single column. So, we conservatively
        # set a default `time_indices = 1:2`.
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
        time_indices_in_memory = 1:length(backend)
        native_fts_architecture = architecture
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
    λn = Array(ds["lon_bnds"][1, :])
    φn = Array(ds["lat_bnds"][1, :])

    # The netCDF coordinates lon_bnds and lat_bnds do not include
    # the last interfaces, so we push them here.
    push!(φn, 90)
    push!(λn, λn[1] + 360)

    i₁, i₂, j₁, j₂, TX = compute_bounding_indices(longitude, latitude, nothing, Center, Center, λc, φc)

    Δλ = (λn[2] - λn[1]) / 2
    Δφ = (φn[2] - φn[1]) / 2

    λr = λn[i₁:i₂+1]
    φr = φn[j₁:j₂+1]

    λr = [λr[1] - Δλ, λr..., λr[end] + Δλ]
    φr = [φr[1] - Δφ, φr..., φr[end] + Δφ]

    Nrx = length(λr) - 1
    Nry = length(φr) - 1
    close(ds)

    N = (Nrx, Nry)
    H = min.(N, (3, 3))

    JRA55_native_grid = LatitudeLongitudeGrid(native_fts_architecture, FT;
                                              halo = H,
                                              size = N,
                                              longitude = λr,
                                              latitude = φr,
                                              topology = (TX, Bounded, Flat))

    boundary_conditions = FieldBoundaryConditions(JRA55_native_grid, (Center, Center, Nothing))
    start_time = first_date(metadata.dataset, metadata.name)
    times = native_times(metadata; start_time)

    if backend isa JRA55NetCDFBackend
        fts = FieldTimeSeries{Center, Center, Nothing}(JRA55_native_grid, times;
                                                       backend,
                                                       time_indexing,
                                                       boundary_conditions,
                                                       path = filepath,
                                                       name = shortname)

        set!(fts)
        return fts
    else
        fts = FieldTimeSeries{Center, Center, Nothing}(JRA55_native_grid, times;
                                                       time_indexing,
                                                       backend,
                                                       boundary_conditions)

        # Fill the data in a GPU-friendly manner
        ds = Dataset(filepath)
        data = ds[shortname][i₁:i₂, j₁:j₂, time_indices_in_memory]
        close(ds)

        copyto!(interior(fts, :, :, 1, :), data)
        fill_halo_regions!(fts)

        return fts
    end
end
