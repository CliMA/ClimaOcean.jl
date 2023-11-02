module JRA55

using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: fill_halo_regions!
using NCDatasets

filenames = Dict(
    :shortwave_radiation             => "RYF.rsds.1990_1991.nc",
    :freshwater_river_flux           => "RYF.friver.1990_1991.nc"    # Freshwater fluxes from rivers
    :freshwater_rain_flux            => "RYF.prra.1990_1991.nc"      # Freshwater flux from rainfall
    :freshwater_snow_flux            => "RYF.prsn.1990_1991.nc"      # Freshwater flux from snowfall
    :freshwater_iceberg_flux         => "RYF.licalvf.1990_1991.nc"   # Freshwater flux from calving icebergs
    :specific_humidity               => "RYF.huss.1990_1991.nc"      # Surface specific humidity
    :sea_level_pressure              => "RYF.psl.1990_1991.nc"       # Sea level pressure
    :relative_humidity               => "RYF.rhuss.1990_1991.nc"     # Surface relative humidity
    :downwelling_longwave_radiation  => "RYF.rlds.1990_1991.nc"      # Downwelling longwave radiation
    :downwelling_shortwave_radiation => "RYF.rsds.1990_1991.nc"      # Downwelling shortwave radiation
    :atmospheric_temperature         => "RYF.tas.1990_1991.nc"       # Near-surface air temperature
    :atmospheric_eastward_velocity   => "RYF.uas.1990_1991.nc"       # Eastward near-surface wind
    :atmospheric_westward_velocity   => "RYF.vas.1990_1991.nc"       # Northward near-surface wind
)

variable_names = Dict(
    :shortwave_radiation             => "rsds",
    :freshwater_river_flux           => "friver",   # Freshwater fluxes from rivers
    :freshwater_rain_flux            => "prra",     # Freshwater flux from rainfall
    :freshwater_snow_flux            => "prsn",     # Freshwater flux from snowfall
    :freshwater_iceberg_flux         => "licalvf",  # Freshwater flux from calving icebergs
    :specific_humidity               => "huss",     # Surface specific humidity
    :sea_level_pressure              => "psl",      # Sea level pressure
    :relative_humidity               => "rhuss",    # Surface relative humidity
    :downwelling_longwave_radiation  => "rlds",     # Downwelling longwave radiation
    :downwelling_shortwave_radiation => "rsds"      # Downwelling shortwave radiation
    :atmospheric_temperature         => "tas",      # Near-surface air temperature
    :atmospheric_eastward_velocity   => "uas",      # Eastward near-surface wind
    :atmospheric_westward_velocity   => "vas",      # Northward near-surface wind
)

urls = Dict(
    :shortwave_radiation => "https://www.dropbox.com/scl/fi/z6fkvmd9oe3ycmaxta131/" *
                            "RYF.rsds.1990_1991.nc?rlkey=r7q6zcbj6a4fxsq0f8th7c4tc&dl=0"

    :freshwater_river_flux => "https://www.dropbox.com/scl/fi/21ggl4p74k4zvbf04nb67/" * 
                              "RYF.friver.1990_1991.nc?rlkey=ny2qcjkk1cfijmwyqxsfm68fz&dl=0"

    :freshwater_rain_flux => "https://www.dropbox.com/scl/fi/5icl1gbd7f5hvyn656kjq/" *
                             "RYF.prra.1990_1991.nc?rlkey=iifyjm4ppwyd8ztcek4dtx0k8&dl=0"

    :freshwater_snow_flux => "https://www.dropbox.com/scl/fi/1r4ajjzb3643z93ads4x4/" *
                             "RYF.prsn.1990_1991.nc?rlkey=auyqpwn060cvy4w01a2yskfah&dl=0"

    :freshwater_iceberg_flux => "https://www.dropbox.com/scl/fi/44nc5y27ohvif7lkvpyv0/" *
                                "RYF.licalvf.1990_1991.nc?rlkey=w7rqu48y2baw1efmgrnmym0jk&dl=0"

    :specific_humidity => "https://www.dropbox.com/scl/fi/66z6ymfr4ghkynizydc29/" *
                          "RYF.huss.1990_1991.nc?rlkey=107yq04aew8lrmfyorj68v4td&dl=0"

    :sea_level_pressure => "https://www.dropbox.com/scl/fi/0fk332027oru1iiseykgp/" *
                           "RYF.psl.1990_1991.nc?rlkey=4xpr9uah741483aukok6d7ctt&dl=0"

    :relative_humidity => "https://www.dropbox.com/scl/fi/1agwsp0lzvntuyf8bm9la/" *
                          "RYF.rhuss.1990_1991.nc?rlkey=8cd0vs7iy1rw58b9pc9t68gtz&dl=0"

    :downwelling_longwave_radiation  => "https://www.dropbox.com/scl/fi/y6r62szkirrivua5nqq61/" *
                                        "RYF.rlds.1990_1991.nc?rlkey=wt9yq3cyrvs2rbowoirf4nkum&dl=0"

    :downwelling_shortwave_radiation => "https://www.dropbox.com/scl/fi/z6fkvmd9oe3ycmaxta131/" *
                                        "RYF.rsds.1990_1991.nc?rlkey=r7q6zcbj6a4fxsq0f8th7c4tc&dl=0"

    :atmospheric_temperature => "https://www.dropbox.com/scl/fi/fpl0npwi476w635g6lke9/" *
                                "RYF.tas.1990_1991.nc?rlkey=0skb9pe6lgbfbiaoybe7m945s&dl=0"

    :atmospheric_eastward_velocity => "https://www.dropbox.com/scl/fi/86wetpqla2x97isp8092g/" *
                                      "RYF.uas.1990_1991.nc?rlkey=rcaf18sh1yz0v9g4hjm1249j0&dl=0"

    :atmospheric_westward_velocity => "https://www.dropbox.com/scl/fi/d38sflo9ddljstd5jwgml/" *
                                      "RYF.vas.1990_1991.nc?rlkey=f9y3e57kx8xrb40gbstarf0x6&dl=0"
)

"""
    jra55_field_time_series(name, architecture=CPU();
                            time_indices = :,    
                            url = urls[name],
                            variable_name = variable_names[name])

Return a FieldTimeSeries representing JRA55 data at the interface between the
atmosphere and ocean.
"""
function jra55_field_time_series(name, arch=CPU();
                                 time_indices = :,    
                                 url = urls[name],
                                 filename = filenames[name],
                                 variable_name = variable_names[name])

    isfile(filename) || download(url, filename)

    ds = Dataset(filename)

    # Extract variable data
    data = ds[variable_name][:, :, time_indices]

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

    grid = LatitudeLongitudeGrid(arch,
                                 size = (Nx, Ny);
                                 longitude = λ,
                                 latitude = φ,
                                 topology = (Periodic, Bounded, Flat))

    # Hack together the `times` for the JRA55 dataset we are currently using.
    # In the future, when we have the "real" (total) JRA55 dataset, we'll have to
    # use the built-in dates.
    Δt = 3hours
    Nt = length(times)
    start_time = 0
    stop_time = Δt * (Nt - 1)
    times = start_time:Δt:stop_time

    boundary_conditions = FieldBoundaryConditions(grid, (Center, Center, Nothing))
    fts = FieldTimeSeries{Center, Center, Nothing}(grid, times; boundary_conditions)

    # Fill the data
    interior(fts, :, :, 1, :) .= data[:, :, :]

    # Fill halo regions so we can interpolate to finer grids
    Nt = length(times)
    fill_halo_regions!(fts)

    return fts
end

end # module

