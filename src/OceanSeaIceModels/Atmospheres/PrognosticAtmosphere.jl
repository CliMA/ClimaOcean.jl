struct PrognosticAtmosphere{FT, M, T, U, P, C, F, R, TP}
    model :: M
    velocities :: U
    pressure :: P
    tracers :: C
    freshwater_flux :: F
    downwelling_radiation :: R
    thermodynamics_parameters :: TP
    reference_height :: FT
    boundary_layer_height :: FT
end

function PrognosticAtmosphere(arch, grid; model = nothing)
    
    # Extract surface values from the atmospheric model
    # and move them to the correct architecture if needed
    velocities = surface_winds(model, arch)
    pressure = surface_pressure(model, arch)
    tracers = surface_tracers(model, arch)
    freshwater_flux = precipitation(model, arch)

    return PrognosticAtmosphere(model, 
                                velocities, 
                                pressure, 
                                tracers, 
                                freshwater_flux, 
                                downwelling_radiation,
                                thermodynamics_parameters,
                                reference_height,
                                boundary_layer_height)
end

update_model_field_time_series!(atmos::PrognosticAtmosphere) = update_surface_fields!(atmos)

# Out-source the time_step! to the prognostic atmosphere model
time_step!(atmos::PrognosticAtmosphere) = time_step!(atmos.model)

function interpolate_atmospheric_state!(atmosphere::PrognosticAtmosphere, surface_atmosphere_state, grid, clock)

    # Here we need to interpolate atmospheric state to the ocean grid?


end


