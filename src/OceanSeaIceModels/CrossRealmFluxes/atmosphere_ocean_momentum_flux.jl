using Oceananigans.Operators: ℑxᶜᵃᵃ, ℑyᵃᶜᵃ

#####
##### Relative velocities
#####

@inline function x_atmosphere_ocean_momentum_flux(i, j, grid, clock, cᴰ, ρₒ, ρₐ, Uₒ, Uₐ)
    x_transfer_velocity_transfer_velocity = Δu_transfer_velocity(i, j, grid, time, Uₒ.u, Uₒ.v, Uₐ.u, Uₐ.v)
    return ρₐ / ρₒ * cᴰ * x_transfer_velocity_transfer_velocity
end 

@inline function y_atmosphere_ocean_momentum_flux(i, j, grid, clock, cᴰ, ρₒ, ρₐ, Uₒ, Uₐ)
    time = Time(clock.time)
    y_transfer_velocity_transfer_velocity = Δv_transfer_velocity(i, j, grid, time, Uₒ.u, Uₒ.v, Uₐ.u, Uₐ.v)
    return ρₐ / ρₒ * cᴰ * y_transfer_velocity_transfer_velocity
end

    
