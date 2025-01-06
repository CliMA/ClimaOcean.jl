struct PrognosticAtmosphere{FT, S, TP}
    simulation :: S
    thermodynamics_parameters :: TP # TODO: Make sure these are consistent with `simulation`
    reference_height :: FT
    boundary_layer_height :: FT
end

function PrognosticAtmosphere(FT=Float64; 
                              simulation=nothing,
                              reference_height = 10, # meters
                              boundary_layer_height = 600, # meters,
                              thermodynamics_parameters = PrescribedAtmosphereThermodynamicsParameters(FT))
    
    return PrognosticAtmosphere(simulation, thermodynamics_parameters, reference_height, boundary_layer_height)
end

# Nothing for the moment, to be fixed by the various atmospheric models
update_model_field_time_series!(::PrognosticAtmosphere, time) = nothing

# To be extended by the various atmospheric models
# The idea is that the atmospheric model provides the surface data and interpolates it into the
# surface_atmosphere_state that represents surface atmospheric data on the ocean/sea-ice grid
interpolate_atmospheric_state!(surface_atmosphere_state, 
                               auxiliary_freshwater_flux, 
                               ::PrognosticAtmosphere, 
                               grid, clock) = nothing

# To be extended by the various atmospheric models.
# Convert the fluxes from the ocean/sea-ice grid to the atmospheric model grid
regrid_fluxes_to_atmospheric_model!(::PrognosticAtmosphere, net_tracer_fluxes, centered_velocity_fluxes) = nothing

Base.show(io::IO, atmos::PrognosticAtmosphere) = 
    print(io, "Prognostic Atmosphere: ", '\n',
              "├── simulation: ", summary(atmos.simulation), '\n',
              "├── thermodynamics_parameters: ", summary(atmos.thermodynamics_parameters), '\n',
              "├── reference_height: ", atmos.reference_height, '\n',
              "└── boundary_layer_height: ", atmos.boundary_layer_height)