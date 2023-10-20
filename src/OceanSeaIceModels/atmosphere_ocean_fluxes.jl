# If there is no atmosphere, do not compute fluxes! (this model has the ocean component which 
# will handle the top boundary_conditions, for example if we want to impose a value BC)
compute_air_sea_flux!(coupled_model::NoAtmosphereModel) = nothing

function compute_air_sea_flux!(coupled_model)
    ocean   = coupled_model.ocean
    forcing = coupled_model.atmospheric_forcing

    (; T, S) = ocean.model.tracers
    (; u, v) = ocean.model.velocities

    grid   = ocean.model.grid
    clock  = ocean.model.clock
    fields = prognostic_fields(ocean.model)

    Qˢ = T.boundary_conditions.top.condition
    Fˢ = S.boundary_conditions.top.condition
    τˣ = u.boundary_conditions.top.condition
    τʸ = v.boundary_conditions.top.condition

    ε  = coupled_model.ocean_emissivity
    ρₒ = coupled_model.ocean_density
    cₒ = coupled_model.ocean_heat_capacity
    I₀ = coupled_model.solar_insolation

    ice_thickness = coupled_model.ice.model.ice_thickness

    launch!(ocean, :xy, _calculate_air_sea_fluxes!, Qˢ, Fˢ, τˣ, τʸ, ρₒ, cₒ, ε, I₀, 
            grid, clock, fields, forcing, ice_thickness)

    return nothing
end


