#####
##### Functions extended by sea-ice and ocean models
#####

reference_density(::Nothing) = 0
heat_capacity(::Nothing) = 0

#####
##### Functions extended by sea-ice models
#####

sea_ice_thickness(::Nothing) = ZeroField()
sea_ice_concentration(::Nothing) = ZeroField()

#####
##### Functions extended by atmosphere models
#####

function thermodynamics_parameters end
function surface_layer_height end
function boundary_layer_height end

#####
##### Functions extended by all component models
#####

function interpolate_state! end
function compute_net_fluxes! end

# Fallbacks for a ``Nothing`` component model
compute_net_fluxes!(coupled_model, ::Nothing) = nothing
interpolate_state!(exchanger, grid, ::Nothing, coupled_model) = nothing
