const OceanOnlyModel = OceanSeaIceModel{Nothing}

OceanOnlyModel(ocean; kw...) = OceanSeaIceModel(nothing, ocean; kw...)

#####
##### No ice-ocean fluxes in this model!!
#####

function time_step!(coupled_model::OceanOnlyModel, Δt; callbacks=[], compute_tendencies=true)
    ocean = coupled_model.ocean

    # Be paranoid and update state at iteration 0
    coupled_model.clock.iteration == 0 && update_state!(coupled_model, callbacks)

    time_step!(ocean)

    tick!(coupled_model.clock, ocean.Δt) # An Ocean-only model advances with the ocean time-step!
    update_state!(coupled_model, callbacks; compute_tendencies)
    
    return nothing
end

#=
compute_ice_ocean_salinity_flux!(::OceanOnlyModel) = nothing
ice_ocean_latent_heat!(::OceanOnlyModel) = nothing

#####
##### Air-sea fluxes
#####

function time_step!(coupled_model::OceanOnlyModel, Δt; callbacks=nothing)
    time_step!(ocean)
    tick!(coupled_model.clock, Δt)
    return nothing
end

function update_state!(coupled_model::OceanOnlyModel; callbacks=nothing)
    compute_air_sea_flux!(coupled_model)
    return nothing
end

function compute_air_sea_fluxes!(coupled_model::OceanOnlyModel) 
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
    ρₒ = coupled_model.ocean_reference_density
    cₒ = coupled_model.ocean_heat_capacity
    I₀ = coupled_model.solar_insolation

    launch!(ocean, :xy, _calculate_air_sea_fluxes!, Qˢ, Fˢ, τˣ, τʸ, ρₒ, cₒ, ε, Iₒ, 
            grid, clock, fields, forcing, nothing)

    return nothing
end
=#
