struct PrognosticAtmosphere{FT, S, TP}
    simulation :: S
    thermodynamics_parameters :: TP
    reference_height :: FT
    boundary_layer_height :: FT
end

update_model_field_time_series!(atmos::PrognosticAtmosphere) = update_surface_fields!(atmos)

# Out-source the time_step! to the prognostic atmosphere model
time_step!(atmos::PrognosticAtmosphere) = time_step!(atmos.simulation)

# To be extended by the various atmospheric models
# The idea is that the atmospheric model provides the surface data and interpolates it into the
# surface_atmosphere_state that represents surface atmospheric data on the ocean/sea-ice grid
interpolate_atmospheric_state!(atmosphere::PrognosticAtmosphere, surface_atmosphere_state, grid, clock) = nothing

# To be extended by the various atmospheric models.
# Convert the fluxes from the ocean/sea-ice grid to the atmospheric model grid
regrid_fluxes_to_atmospheric_model!(atmosphere::PrognosticAtmosphere, net_tracer_fluxes, centered_velocity_fluxes) = nothing