module JRA55

using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: fill_halo_regions!
using NCDatasets

# A list of all variables provided in the JRA55 dataset:
jra55_short_names = (:freshwater_river_flux,
                     :freshwater_rain_flux,
                     :freshwater_snow_flux,
                     :freshwater_iceberg_flux,
                     :specific_humidity,
                     :sea_level_pressure,
                     :relative_humidity,
                     :downwelling_longwave_radiation,
                     :downwelling_shortwave_radiation,
                     :atmospheric_temperature,
                     :atmospheric_eastward_velocity,
                     :atmospheric_northward_velocity)

file_names = Dict(
    :freshwater_river_flux           => "RYF.friver.1990_1991.nc",   # Freshwater fluxes from rivers
    :freshwater_rain_flux            => "RYF.prra.1990_1991.nc",     # Freshwater flux from rainfall
    :freshwater_snow_flux            => "RYF.prsn.1990_1991.nc",     # Freshwater flux from snowfall
    :freshwater_iceberg_flux         => "RYF.licalvf.1990_1991.nc",  # Freshwater flux from calving icebergs
    :specific_humidity               => "RYF.huss.1990_1991.nc",     # Surface specific humidity
    :sea_level_pressure              => "RYF.psl.1990_1991.nc",      # Sea level pressure
    :relative_humidity               => "RYF.rhuss.1990_1991.nc",    # Surface relative humidity
    :downwelling_longwave_radiation  => "RYF.rlds.1990_1991.nc",     # Downwelling longwave radiation
    :downwelling_shortwave_radiation => "RYF.rsds.1990_1991.nc",     # Downwelling shortwave radiation
    :atmospheric_temperature         => "RYF.tas.1990_1991.nc",      # Near-surface air temperature
    :atmospheric_eastward_velocity   => "RYF.uas.1990_1991.nc",      # Eastward near-surface wind
    :atmospheric_northward_velocity  => "RYF.vas.1990_1991.nc",      # Northward near-surface wind
)

short_names = Dict(
    :freshwater_river_flux           => "friver",   # Freshwater fluxes from rivers
    :freshwater_rain_flux            => "prra",     # Freshwater flux from rainfall
    :freshwater_snow_flux            => "prsn",     # Freshwater flux from snowfall
    :freshwater_iceberg_flux         => "licalvf",  # Freshwater flux from calving icebergs
    :specific_humidity               => "huss",     # Surface specific humidity
    :sea_level_pressure              => "psl",      # Sea level pressure
    :relative_humidity               => "rhuss",    # Surface relative humidity
    :downwelling_longwave_radiation  => "rlds",     # Downwelling longwave radiation
    :downwelling_shortwave_radiation => "rsds",     # Downwelling shortwave radiation
    :atmospheric_temperature         => "tas",      # Near-surface air temperature
    :atmospheric_eastward_velocity   => "uas",      # Eastward near-surface wind
    :atmospheric_northward_velocity  => "vas",      # Northward near-surface wind
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

    :atmospheric_temperature => "https://www.dropbox.com/scl/fi/fpl0npwi476w635g6lke9/" *
                                "RYF.tas.1990_1991.nc?rlkey=0skb9pe6lgbfbiaoybe7m945s&dl=0",

    :atmospheric_eastward_velocity => "https://www.dropbox.com/scl/fi/86wetpqla2x97isp8092g/" *
                                      "RYF.uas.1990_1991.nc?rlkey=rcaf18sh1yz0v9g4hjm1249j0&dl=0",

    :atmospheric_northward_velocity => "https://www.dropbox.com/scl/fi/d38sflo9ddljstd5jwgml/" *
                                       "RYF.vas.1990_1991.nc?rlkey=f9y3e57kx8xrb40gbstarf0x6&dl=0",
)

"""
    jra55_field_time_series(variable_name;
                            architecture = CPU(),
                            time_indices = :,    
                            url = urls[name],
                            filename = file_names[variable_name],
                            short_name = short_names[variable_name])

Return a `FieldTimeSeries` containing atmospheric reanalysis data for `variable_name`,
which describes one of the variables in the "repeat year forcing" dataset derived
from the Japanese 55-year atmospheric reanalysis for driving ocean-sea-ice models (JRA55-do).
For more information about the derivation of the repeat year forcing dataset, see

"Stewart et al., JRA55-do-based repeat year forcing datasets for driving ocean–sea-ice models",
Ocean Modelling, 2020, https://doi.org/10.1016/j.ocemod.2019.101557.

The `variable_name`s (and their `short_name`s used in NetCDF files)
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
    - `:atmospheric_temperature`            ("ras")
    - `:atmospheric_eastward_velocity`      ("uas")
    - `:atmospheric_northward_velocity`     ("vas")

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

    - `short_name`: The "short name" of `variable_name` inside its NetCDF file.
                    Default: `ClimaOcean.JRA55.short_names[variable_name]`.
"""
function jra55_field_time_series(variable_name;
                                 architecture = CPU(),
                                 time_indices = :,    
                                 url = urls[variable_name],
                                 filename = file_names[variable_name],
                                 short_name = short_names[variable_name])

    isfile(filename) || download(url, filename)

    ds = Dataset(filename)

    # Note that each file should have the variables
    # ds["time"]:     time coordinate 
    # ds["lon"]:      longitude at the location of the variable
    # ds["lat"]:      latitude at the location of the variable
    # ds["lon_bnds"]: bounding longitudes between which variables are averaged
    # ds["lat_bnds"]: bounding latitudes between which variables are averaged
    # ds[short_name]: the variable data

    # Extract variable data
    data = ds[short_name][:, :, time_indices]

    # Make the JRA55 grid
    λ = ds["lon_bnds"][1, :]
    φ = ds["lat_bnds"][1, :]
    times = ds["time"][time_indices]

    close(ds)

    # The .nc coordinates lon_bnds and lat_bnds do not include
    # the last interface, so we push them here.
    push!(φ, 90)
    push!(λ, λ[1] + 360)

    Nx = length(λ) - 1
    Ny = length(φ) - 1

    grid = LatitudeLongitudeGrid(architecture,
                                 size = (Nx, Ny);
                                 longitude = λ,
                                 latitude = φ,
                                 topology = (Periodic, Bounded, Flat))

    # Hack together the `times` for the JRA55 dataset we are currently using.
    # We might want to use the acutal dates instead though.
    # So the following code maybe should change.
    Δt = 3hours # just what it is
    Nt = length(times)
    start_time = 0 # Note: the forcing start at Jan 1 of the repeat year.
    stop_time = Δt * (Nt - 1)
    times = start_time:Δt:stop_time

    boundary_conditions = FieldBoundaryConditions(grid, (Center, Center, Nothing))
    fts = FieldTimeSeries{Center, Center, Nothing}(grid, times; boundary_conditions)

    # Fill the data
    interior(fts, :, :, 1, :) .= data[:, :, :]

    # Fill halo regions so we can interpolate to finer grids
    Nt = length(times)
    chunk_size = 100
    if Nt <= chunk_size # one chunk will do
        fill_halo_regions!(fts)
    else # need multiple chunks
        start = 1
        while start < Nt
            stop = min(Nt, start + chunk_size - 1)
            fts_chunk = Tuple(fts[n] for n = start:stop)
            fill_halo_regions!(fts_chunk)
            start += chunk_size
        end
    end

    return fts
end

end # module

