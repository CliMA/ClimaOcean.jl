module JRA55

using Oceananigans
using Oceananigans.Units

using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: λnodes, φnodes
using Oceananigans.Fields: interpolate!

using ClimaOcean.OceanSeaIceModels:
    PrescribedAtmosphere,
    TwoStreamDownwellingRadiation

using NCDatasets

# A list of all variables provided in the JRA55 dataset:
jra55_variable_names = (:freshwater_river_flux,
                        :freshwater_rain_flux,
                        :freshwater_snow_flux,
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
    :freshwater_rain_flux            => "RYF.prra.1990_1991.nc",     # Freshwater flux from rainfall
    :freshwater_snow_flux            => "RYF.prsn.1990_1991.nc",     # Freshwater flux from snowfall
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

shortnames = Dict(
    :freshwater_river_flux           => "friver",   # Freshwater fluxes from rivers
    :freshwater_rain_flux            => "prra",     # Freshwater flux from rainfall
    :freshwater_snow_flux            => "prsn",     # Freshwater flux from snowfall
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

urls = Dict(
    :shortwave_radiation => "https://www.dropbox.com/scl/fi/z6fkvmd9oe3ycmaxta131/" *
                            "RYF.rsds.1990_1991.nc?rlkey=r7q6zcbj6a4fxsq0f8th7c4tc&dl=0",

    :freshwater_river_flux => "https://www.dropbox.com/scl/fi/21ggl4p74k4zvbf04nb67/" * 
                              "RYF.friver.1990_1991.nc?rlkey=ny2qcjkk1cfijmwyqxsfm68fz&dl=0",

    :freshwater_rain_flux => "https://www.dropbox.com/scl/fi/5icl1gbd7f5hvyn656kjq/" *
                             "RYF.prra.1990_1991.nc?rlkey=iifyjm4ppwyd8ztcek4dtx0k8&dl=0",

    :freshwater_snow_flux => "https://www.dropbox.com/scl/fi/1r4ajjzb3643z93ads4x4/" *
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

"""
    jra55_field_time_series(variable_name;
                            architecture = CPU(),
                            time_indices = :,    
                            url = urls[name],
                            filename = filenames[variable_name],
                            shortname = shortnames[variable_name])

Return a `FieldTimeSeries` containing atmospheric reanalysis data for `variable_name`,
which describes one of the variables in the "repeat year forcing" dataset derived
from the Japanese 55-year atmospheric reanalysis for driving ocean-sea-ice models (JRA55-do).
For more information about the derivation of the repeat year forcing dataset, see

"Stewart et al., JRA55-do-based repeat year forcing datasets for driving ocean–sea-ice models",
Ocean Modelling, 2020, https://doi.org/10.1016/j.ocemod.2019.101557.

The `variable_name`s (and their `shortname`s used in NetCDF files)
available from the JRA55-do are:

    - `:freshwater_river_flux`              ("friver")
    - `:freshwater_rain_flux`               ("prra")
    - `:freshwater_snow_flux`               ("prsn")
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
                    Default: `ClimaOcean.JRA55.shortnames[variable_name]`.
"""
function jra55_field_time_series(variable_name, grid=nothing;
                                 architecture = CPU(),
                                 location = nothing,
                                 time_indices = :,    
                                 url = nothing,
                                 filename = nothing,
                                 shortname = nothing)

    if isnothing(filename) && !(variable_name ∈ jra55_variable_names)
        variable_strs = Tuple("  - :$name \n" for name in jra55_variable_names)
        variables_msg = prod(variable_strs)

        msg = string("The variable :$variable_name is not provided by the JRA55-do dataset!", '\n',
                     "The variables provided by the JRA55-do dataset are:", '\n',
                     variables_msg)

        throw(ArgumentError(msg))
    end

    isnothing(url)       && (url       = urls[variable_name])
    isnothing(filename)  && (filename  = filenames[variable_name])
    isnothing(shortname) && (shortname = shortnames[variable_name])

    isfile(filename) || download(url, filename)

    # Get location
    if isnothing(location)
        LX = LY = Center
    else
        LX, LY = location
    end

    ds = Dataset(filename)

    # Note that each file should have the variables
    # ds["time"]:     time coordinate 
    # ds["lon"]:      longitude at the location of the variable
    # ds["lat"]:      latitude at the location of the variable
    # ds["lon_bnds"]: bounding longitudes between which variables are averaged
    # ds["lat_bnds"]: bounding latitudes between which variables are averaged
    # ds[shortname]: the variable data

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

    times = ds["time"][time_indices]

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

    # Extract variable data
    data = ds[shortname][i₁:i₂, j₁:j₂, time_indices]
    λr = λn[i₁:i₂+1]
    φr = φn[j₁:j₂+1]
    Nrx, Nry, Nt = size(data)

    close(ds)

    # Hack together the `times` for the JRA55 dataset we are currently using.
    # We might want to use the acutal dates instead though.
    # So the following code might need to change.
    Δt = 3hours # just what it is
    Nt = length(times)
    start_time = 0 # Note: the forcing start at Jan 1 of the repeat year.
    stop_time = Δt * (Nt - 1)
    times = start_time:Δt:stop_time

    jra55_native_grid = LatitudeLongitudeGrid(architecture,
                                              size = (Nrx, Nry);
                                              longitude = λr,
                                              latitude = φr,
                                              topology = (TX, Bounded, Flat))

    boundary_conditions = FieldBoundaryConditions(jra55_native_grid, (Center, Center, Nothing))
    native_fts = FieldTimeSeries{Center, Center, Nothing}(jra55_native_grid, times; boundary_conditions)

    # Fill the data
    copyto!(interior(native_fts, :, :, 1, :), data[:, :, :])

    # Fill halo regions so we can interpolate to finer grids
    Nt = length(times)
    fill_halo_regions!(native_fts)

    if isnothing(grid)
        return native_fts
    else # make a new FieldTimeSeries and interpolate native data onto it.

        boundary_conditions = FieldBoundaryConditions(grid, (LX, LY, Nothing))
        fts = FieldTimeSeries{LX, LY, Nothing}(grid, times; boundary_conditions)

        interpolate!(fts, native_fts)

        return fts
    end
end

# TODO: allow the user to pass dates
function jra55_prescribed_atmosphere(grid, time_indices=:)
    architecture = Oceananigans.architecture(grid)

    u_jra55   = jra55_field_time_series(:eastward_velocity,               grid; time_indices, architecture, location=(Face, Center))
    v_jra55   = jra55_field_time_series(:northward_velocity,              grid; time_indices, architecture, location=(Center, Face)) 
    T_jra55   = jra55_field_time_series(:temperature,                     grid; time_indices, architecture)
    q_jra55   = jra55_field_time_series(:specific_humidity,               grid; time_indices, architecture)
    Fr_jra55  = jra55_field_time_series(:freshwater_rain_flux,            grid; time_indices, architecture)
    Fs_jra55  = jra55_field_time_series(:freshwater_snow_flux,            grid; time_indices, architecture)
    Fv_jra55  = jra55_field_time_series(:freshwater_river_flux,           grid; time_indices, architecture)
    Fi_jra55  = jra55_field_time_series(:freshwater_iceberg_flux,         grid; time_indices, architecture)
    Qlw_jra55 = jra55_field_time_series(:downwelling_longwave_radiation,  grid; time_indices, architecture)
    Qsw_jra55 = jra55_field_time_series(:downwelling_shortwave_radiation, grid; time_indices, architecture)

    times = u_jra55.times

    velocities = (u = u_jra55,
                  v = v_jra55)

    tracers = (T = T_jra55,
               q = q_jra55)

    freshwater_flux = (rain     = Fr_jra55,
                       snow     = Fs_jra55,
                       rivers   = Fv_jra55,
                       icebergs = Fi_jra55)

    downwelling_radiation = TwoStreamDownwellingRadiation(shortwave=Qsw_jra55, longwave=Qsw_jra55)
    atmosphere = PrescribedAtmosphere(times; velocities, freshwater_flux, tracers, downwelling_radiation)

    return atmosphere
end
 

end # module

