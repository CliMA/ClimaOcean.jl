using Oceananigans.Fields: interpolate!

using ClimaOcean.JRA55: jra55_field_time_series

using ClimaOcean.OceanSeaIceModels: 
    PrescribedAtmosphere,
    TwoStreamDownwellingRadiation

function prescribed_jra55_atmosphere(grid, time_indices=:;
                                     reference_height = 2) # meters

    architecture = Oceananigans.architecture(grid)

    u_jra55   = jra55_field_time_series(:eastward_velocity,               grid; time_indices, architecture)
    v_jra55   = jra55_field_time_series(:northward_velocity,              grid; time_indices, architecture)
    T_jra55   = jra55_field_time_series(:temperature,                     grid; time_indices, architecture)
    p_jra55   = jra55_field_time_series(:surface_pressure,                grid; time_indices, architecture)
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

    atmosphere = PrescribedAtmosphere(times; velocities,
                                      freshwater_flux,
                                      tracers,
                                      downwelling_radiation,
                                      reference_height,
                                      pressure = p_jra55)

    return atmosphere
end
