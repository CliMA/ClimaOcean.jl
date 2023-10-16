
#####
##### No ice-ocean fluxes in this model!!
#####

compute_ice_ocean_salinity_flux!(::OnlyOceanModel) = nothing
ice_ocean_latent_heat!(::OnlyOceanModel) = nothing

#####
##### Air-sea fluxes
#####

function time_step!(coupled_model::OnlyOceanModel, Δt; callbacks=nothing)
    compute_air_sea_flux!(coupled_model)
    time_step!(ocean)
    tick!(coupled_model.clock, Δt)
    return nothing
end

function compute_air_sea_fluxes!(coupled_model::OnlyOceanModel) 
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

    launch!(ocean, :xy, _calculate_air_sea_fluxes!, Qˢ, Fˢ, τˣ, τʸ, ρₒ, cₒ, ε, Iₒ, 
            grid, clock, fields, forcing, nothing)

    return nothing
end
