module JRA55

using Oceananigans
using Oceananigans.Units
 
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: λnodes, φnodes, on_architecture
using Oceananigans.Fields: interpolate!

using ClimaOcean.OceanSeaIceModels:
    PrescribedAtmosphere,
    TwoStreamDownwellingRadiation

using NCDatasets
using JLD2 

# A list of all variables provided in the JRA55 dataset:
JRA55_variable_names = (:freshwater_river_flux,
                        :rain_freshwater_flux,
                        :snow_freshwater_flux,
                        :freshwater_iceberg_flux,
                        :specific_humidity,
                        :sea_level_pressure,
                        :relative_humidity,
                        :downwelling_longwave_radiation,
                        :downwelling_shortwave_radiation,
                        :temperature,
                        :eastward_velocity,
                        :northward_velocity)

filenames = Dict(
    :freshwater_river_flux           => "RYF.friver.1990_1991.nc",   # Freshwater fluxes from rivers
    :rain_freshwater_flux            => "RYF.prra.1990_1991.nc",     # Freshwater flux from rainfall
    :snow_freshwater_flux            => "RYF.prsn.1990_1991.nc",     # Freshwater flux from snowfall
    :freshwater_iceberg_flux         => "RYF.licalvf.1990_1991.nc",  # Freshwater flux from calving icebergs
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
    :freshwater_river_flux           => "friver",   # Freshwater fluxes from rivers
    :rain_freshwater_flux            => "prra",     # Freshwater flux from rainfall
    :snow_freshwater_flux            => "prsn",     # Freshwater flux from snowfall
    :freshwater_iceberg_flux         => "licalvf",  # Freshwater flux from calving icebergs
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
    :freshwater_river_flux           => "Fri", # Freshwater fluxes from rivers
    :rain_freshwater_flux            => "Fra", # Freshwater flux from rainfall
    :snow_freshwater_flux            => "Fsn", # Freshwater flux from snowfall
    :freshwater_iceberg_flux         => "Fic", # Freshwater flux from calving icebergs
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

    :freshwater_river_flux => "https://www.dropbox.com/scl/fi/21ggl4p74k4zvbf04nb67/" * 
                              "RYF.friver.1990_1991.nc?rlkey=ny2qcjkk1cfijmwyqxsfm68fz&dl=0",

    :rain_freshwater_flux => "https://www.dropbox.com/scl/fi/5icl1gbd7f5hvyn656kjq/" *
                             "RYF.prra.1990_1991.nc?rlkey=iifyjm4ppwyd8ztcek4dtx0k8&dl=0",

    :snow_freshwater_flux => "https://www.dropbox.com/scl/fi/1r4ajjzb3643z93ads4x4/" *
                             "RYF.prsn.1990_1991.nc?rlkey=auyqpwn060cvy4w01a2yskfah&dl=0",

    :freshwater_iceberg_flux => "https://www.dropbox.com/scl/fi/44nc5y27ohvif7lkvpyv0/" *
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

function compute_bounding_indices(grid, LX, LY, λc, φc)

    Nx = length(λc)
    Ny = length(φc)

    if isnothing(grid)
        i₁, i₂ = (1, Nx)
        j₁, j₂ = (1, Ny)
        TX = Periodic
    else # only load the data we need
        # Nodes where we need to find data
        λg = λnodes(grid, LX())
        φg = φnodes(grid, LY())

        λ₁, λ₂ = extrema(λg)
        φ₁, φ₂ = extrema(φg)

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
        
        i₁ = searchsortedfirst(λc, λ₁)
        j₁ = searchsortedfirst(φc, φ₁)

        i₂ = searchsortedfirst(λc, λ₂)
        j₂ = searchsortedfirst(φc, φ₂)

        i₁ = max(1, i₁ - 1)
        j₁ = max(1, j₁ - 1)

        i₂ = min(Nx, i₂)
        j₂ = min(Ny, j₂)

        TX = λ₂ - λ₁ ≈ 360 ? Periodic : Bounded
    end

    return i₁, i₂, j₁, j₂, TX
end



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

    - `:freshwater_river_flux`              ("friver")
    - `:rain_freshwater_flux`               ("prra")
    - `:snow_freshwater_flux`               ("prsn")
    - `:freshwater_iceberg_flux`            ("licalvf")
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
function JRA55_field_time_series(variable_name; #, grid=nothing;
                                 architecture = CPU(),
                                 location = nothing,
                                 url = nothing,
                                 filename = nothing,
                                 shortname = nothing,
                                 backend = InMemory(),
                                 preprocess_chunk_size = 10,
                                 preprocess_architecture = CPU(),
                                 time_indices = nothing)

    if isnothing(filename) && !(variable_name ∈ JRA55_variable_names)
        variable_strs = Tuple("  - :$name \n" for name in JRA55_variable_names)
        variables_msg = prod(variable_strs)

        msg = string("The variable :$variable_name is not provided by the JRA55-do dataset!", '\n',
                     "The variables provided by the JRA55-do dataset are:", '\n',
                     variables_msg)

        throw(ArgumentError(msg))
    end

    isnothing(shortname) && (shortname = jra55_short_names[variable_name])

    !isnothing(filename) && !isfile(filename) && isnothing(url) &&
        throw(ArgumentError("A filename was provided without a url, but the file does not exist.\n \
                            If intended, please provide both the filename and url that should be used \n \
                            to download the new file."))

    isnothing(filename) && (filename = filenames[variable_name])
    isnothing(url) && (url = urls[variable_name])

    # Decision tree:
    #   1. jld2 file exists?
    #       - yes -> load and return FieldTimeSeries
    #                check time_indices and all that?
    #       - no -> download .nc data if not available

    jld2_filename = string("JRA55_repeat_year_", variable_name, ".jld2")
    fts_name = field_time_series_short_names[variable_name]
    totally_in_memory = backend isa InMemory{Colon}

    if isfile(jld2_filename)
        isnothing(time_indices) && (time_indices = Colon())

        # Infer the `times` before loading data
        temporary_fts = FieldTimeSeries(jld2_filename, fts_name; backend=OnDisk())

        #try
            times = temporary_fts.times[time_indices]
            fts = FieldTimeSeries(jld2_filename, fts_name; backend, architecture, times)
            return fts
        #catch 
        #    if !totally_in_memory # will need to overwrite
        #        msg = string("Cannot use backend=$backend with time_indices=$time_indices", '\n',
        #                     " and the existing $jld2_filename, which does not", '\n',
        #                     " have enough `times`. Delete $jld2_filename in order", '\n',
        #                     " to re-generate it.")
        #        error(msg)
        #    end
        #end
    end

    isfile(filename) || download(url, filename)

    # Extract variable data
    if totally_in_memory
        # Set a sensible default
        isnothing(time_indices) && (time_indices = 1:1)
        time_indices_in_memory = time_indices
        native_fts_architecture = architecture
    else
        isnothing(time_indices) && (time_indices = :)
        time_indices_in_memory = 1:preprocess_chunk_size
        native_fts_architecture = preprocess_architecture
    end

    # Get location
    if isnothing(location)
        LX = LY = Center
    else
        LX, LY = location
    end

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
    # i₁, i₂, j₁, j₂, TX = compute_bounding_indices(grid, LX, LY, λc, φc)
    Nx = length(λc)
    Ny = length(φc)
    i₁, i₂ = (1, Nx)
    j₁, j₂ = (1, Ny)
    TX = Periodic

    times = ds["time"][time_indices_in_memory]
    data = ds[shortname][i₁:i₂, j₁:j₂, time_indices_in_memory]
    λr = λn[i₁:i₂+1]
    φr = φn[j₁:j₂+1]
    Nrx, Nry, Nt = size(data)
    close(ds)

    JRA55_native_grid = LatitudeLongitudeGrid(native_fts_architecture, Float32;
                                              halo = (1, 1),
                                              size = (Nrx, Nry),
                                              longitude = λr,
                                              latitude = φr,
                                              topology = (TX, Bounded, Flat))

    boundary_conditions = FieldBoundaryConditions(JRA55_native_grid, (Center, Center, Nothing))

    # Hack together the `times` for the JRA55 dataset we are currently using.
    # We might want to use the acutal dates instead though.
    # So the following code might need to change.
    Δt = 3hours # just what it is
    Nt = length(times)
    start_time = 0 # Note: the forcing start at Jan 1 of the repeat year.
    stop_time = Δt * (Nt - 1)
    times = start_time:Δt:stop_time

    # Make times into an array for later preprocessing
    !totally_in_memory && (times = collect(times))

    native_fts = FieldTimeSeries{Center, Center, Nothing}(JRA55_native_grid, times; boundary_conditions)

    # Fill the data in a GPU-friendly manner
    copyto!(interior(native_fts, :, :, 1, :), data)
    fill_halo_regions!(native_fts)

    #=
    if isnothing(grid)
        fts = native_fts
    else # make a new FieldTimeSeries and interpolate native data onto it.
        boundary_conditions = FieldBoundaryConditions(grid, (LX, LY, Nothing))
        fts = FieldTimeSeries{LX, LY, Nothing}(grid, times; boundary_conditions)
        interpolate!(fts, native_fts)
    end
    =#

    if totally_in_memory
        return native_fts
    else # we're gonna save to disk!
        @info "Pre-processing JRA55 data into a JLD2 file to be used with FieldTimeSeries..."

        on_disk_fts = FieldTimeSeries{LX, LY, Nothing}(native_fts.grid;
                                                       backend = OnDisk(),
                                                       path = jld2_filename,
                                                       name = fts_name)

        # Re-open the dataset!
        ds = Dataset(filename)
        all_datetimes = ds["time"][time_indices]
        all_Nt = length(all_datetimes)
        chunk = last(preprocess_chunk_size)

        Δt = 3hours # just what it is
        start_time = 0 # Note: the forcing starts at Jan 1 of the repeat year.
        stop_time = Δt * (all_Nt - 1)
        all_times = start_time:Δt:stop_time

        # Save data to disk, one field at a time
        start_clock = time_ns()
        n = 1 # on disk
        m = 0 # in memory
        while n <= all_Nt
            print("        ... processing time index $n of $all_Nt \r")

            if time_indices_in_memory isa Colon || n ∈ time_indices_in_memory
                m += 1
            else
                # Update time_indices
                time_indices_in_memory = time_indices_in_memory .+ preprocess_chunk_size

                # Clip time_indices if they extend past the end of the dataset
                if last(time_indices_in_memory) > all_Nt
                    n₁ = first(time_indices_in_memory)
                    time_indices_in_memory = UnitRange(n₁, all_Nt)
                end

                # Re-load .nc times and data
                # new_times = ds["time"][time_indices_in_memory]
                # native_fts.times .= new_times

                new_data = ds[shortname][i₁:i₂, j₁:j₂, time_indices_in_memory]
                copyto!(interior(native_fts, :, :, 1, :), new_data[:, :, :])
                fill_halo_regions!(native_fts)

                m = 1 # reset
            end

            set!(on_disk_fts, native_fts[m], n, all_times[n])
            n += 1
        end

        elapsed = 1e-9 * (time_ns() - start_clock)
        elapsed_str = prettytime(elapsed)
        @info "    ... done ($elapsed_str)" * repeat(" ", 20)

        close(ds)

        grid = on_architecture(architecture, JRA55_native_grid)

        backend_fts = FieldTimeSeries{LX, LY, Nothing}(grid, all_times;
                                                       backend,
                                                       boundary_conditions,
                                                       path = jld2_filename,
                                                       name = fts_name)

        return backend_fts
    end
end

const AA = Oceananigans.Architectures.AbstractArchitecture

JRA55_prescribed_atmosphere(time_indices=Colon(); kw...) =
    JRA55_prescribed_atmosphere(CPU(), time_indices; kw...)

# TODO: allow the user to pass dates
<<<<<<< HEAD
function JRA55_prescribed_atmosphere(architecture::AA, time_indices=Colon();
                                     backend = InMemory(24), # 3 days of data
                                     reference_height = 2,  # meters
                                     other_kw...)

    ua  = JRA55_field_time_series(:eastward_velocity;               time_indices, backend, architecture, other_kw...)
    va  = JRA55_field_time_series(:northward_velocity;              time_indices, backend, architecture, other_kw...)
    Ta  = JRA55_field_time_series(:temperature;                     time_indices, backend, architecture, other_kw...)
    qa  = JRA55_field_time_series(:specific_humidity;               time_indices, backend, architecture, other_kw...)
    pa  = JRA55_field_time_series(:sea_level_pressure;              time_indices, backend, architecture, other_kw...)
    Fra = JRA55_field_time_series(:rain_freshwater_flux;            time_indices, backend, architecture, other_kw...)
    Fsn = JRA55_field_time_series(:snow_freshwater_flux;            time_indices, backend, architecture, other_kw...)
    Ql  = JRA55_field_time_series(:downwelling_longwave_radiation;  time_indices, backend, architecture, other_kw...)
    Qs  = JRA55_field_time_series(:downwelling_shortwave_radiation; time_indices, backend, architecture, other_kw...)

    # NOTE: these have a different frequency than 3 hours so some changes are needed to 
    # JRA55_field_time_series to support them.
    # Fv_JRA55  = JRA55_field_time_series(:freshwater_river_flux,           grid; time_indices, architecture)
    # Fi_JRA55  = JRA55_field_time_series(:freshwater_iceberg_flux,         grid; time_indices, architecture)

    times = ua.times

    velocities = (u = ua,
                  v = va)

    tracers = (T = Ta,
               q = qa)

    freshwater_flux = (rain     = Fra,
                       snow     = Fsn)
                       # rivers   = Fv_JRA55,
                       # icebergs = Fi_JRA55)
                       
    pressure = pa

    downwelling_radiation = TwoStreamDownwellingRadiation(shortwave=Qs, longwave=Ql)

    atmosphere = PrescribedAtmosphere(times, eltype(ua);
                                      velocities,
                                      freshwater_flux,
                                      tracers,
                                      downwelling_radiation,
                                      reference_height,
                                      pressure)

    return atmosphere
end

end # module

