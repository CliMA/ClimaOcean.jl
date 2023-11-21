#####
##### Setup JRA55 atmosphere
#####

time_indices = 1:10

u_jra55_native = jra55_field_time_series(:eastward_velocity;  time_indices, architecture=arch)
v_jra55_native = jra55_field_time_series(:northward_velocity; time_indices, architecture=arch)

T_jra55_native = jra55_field_time_series(:temperature;  time_indices, architecture=arch)
q_jra55_native = jra55_field_time_series(:relative_humidity;  time_indices, architecture=arch)

Fr_jra55_native = jra55_field_time_series(:freshwater_rain_flux;  time_indices, architecture=arch)
Fs_jra55_native = jra55_field_time_series(:freshwater_snow_flux;  time_indices, architecture=arch)
#Fv_jra55_native = jra55_field_time_series(:freshwater_river_flux;  time_indices, architecture=arch)
#Fi_jra55_native = jra55_field_time_series(:freshwater_iceberg_flux;  time_indices, architecture=arch)
                                                          
times = u_jra55_native.times

u_bcs = FieldBoundaryConditions(grid, (Face, Center, Nothing))
v_bcs = FieldBoundaryConditions(grid, (Center, Face, Nothing))
c_bcs = FieldBoundaryConditions(grid, (Center, Center, Nothing))

u_jra55 = FieldTimeSeries{Face, Center, Nothing}(grid, times; boundary_conditions=u_bcs)
v_jra55 = FieldTimeSeries{Center, Face, Nothing}(grid, times; boundary_conditions=v_bcs)

T_jra55  = FieldTimeSeries{Center, Center, Nothing}(grid, times; boundary_conditions=c_bcs)
q_jra55  = FieldTimeSeries{Center, Center, Nothing}(grid, times; boundary_conditions=c_bcs)
Fr_jra55 = FieldTimeSeries{Center, Center, Nothing}(grid, times; boundary_conditions=c_bcs)
Fs_jra55 = FieldTimeSeries{Center, Center, Nothing}(grid, times; boundary_conditions=c_bcs)

interpolate!(u_jra55,  u_jra55_native)
interpolate!(v_jra55,  v_jra55_native)
interpolate!(T_jra55,  T_jra55_native)
interpolate!(q_jra55,  q_jra55_native)
interpolate!(Fr_jra55, Fr_jra55_native)
interpolate!(Fs_jra55, Fs_jra55_native)

velocities = (u = u_jra55,
              v = v_jra55)

tracers = (T = T_jra55,
           q = q_jra55)

freshwater_fluxes = (rain = Fr_jra55,
                     snow = Fs_jra55)

atmosphere = PrescribedAtmosphere(velocities, freshwater_fluxes, tracers, times)

#####
##### Radiation
#####

Qlw_jra55_native = jra55_field_time_series(:downwelling_longwave_radiation;  time_indices, architecture=arch)
Qsw_jra55_native = jra55_field_time_series(:downwelling_shortwave_radiation; time_indices, architecture=arch)

Qlw_jra55 = FieldTimeSeries{Center, Center, Nothing}(grid, times; boundary_conditions=c_bcs)
Qsw_jra55 = FieldTimeSeries{Center, Center, Nothing}(grid, times; boundary_conditions=c_bcs)

interpolate!(Qlw_jra55, Qlw_jra55_native)
interpolate!(Qsw_jra55, Qsw_jra55_native)

radiation = Radiation(downwelling_shortwave_radiation = Qsw_jra55,
                       downwelling_longwave_radiation = Qlw_jra55)

