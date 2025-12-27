#####
##### Functions extended by ocean models
#####

reference_density(::Nothing) = 0
heat_capacity(::Nothing) = 0
ocean_temperature(ocean) = ZeroField()
ocean_salinity(ocean) = ZeroField()
ocean_surface_salinity(ocean) = ZeroField()
ocean_surface_velocities(ocean) = ZeroField(), ZeroField()

#####
##### Functions extended by sea-ice models
#####

sea_ice_thickness(sea_ice) = ZeroField()
sea_ice_concentration(sea_ice) = ZeroField()
sea_ice_top_temperature(sea_ice) = ZeroField()
sea_ice_salinity(sea_ice) = ZeroField()
sea_ice_ocean_stress(sea_ice) = nothing
sea_ice_velocities(sea_ice) = ZeroField(), ZeroField()
liquidus(sea_ice) = LinearLiquidus()

function default_sea_ice end
function default_ai_temperature end

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
function update_net_fluxes! end

# Fallbacks for a  generic component model
update_net_fluxes!(coupled_model, component) = nothing
interpolate_state!(exchanger, grid, component, coupled_model) = nothing
