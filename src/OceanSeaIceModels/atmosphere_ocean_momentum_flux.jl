using Oceananigans.Operators: ℑxᶜᵃᵃ, ℑyᵃᶜᵃ

#####
##### Relative velocities
#####

# Compute the square of a FieldTimeSeries at `time`
@inline Δϕt²(i, j, k, grid, ϕ1t, ϕ2, time) = @inbounds (ϕ1t[i, j, k, time] - ϕ2[i, j, k])^2

const ThreeDArray = AbstractArray{<:Any, 3}
const FourDArray = AbstractArray{<:Any, 4}

@inline transfer_velocityᶠᶜᶜ(i, j, grid, time, Uₒ, Uₐ) = transfer_velocityᶠᶜᶜ(i, j, grid, time, Uₒ.u, Uₒ.v, Uₐ.u, Uₐ.v)
@inline transfer_velocityᶜᶠᶜ(i, j, grid, time, Uₒ, Uₐ) = transfer_velocityᶜᶠᶜ(i, j, grid, time, Uₒ.u, Uₒ.v, Uₐ.u, Uₐ.v)

@inline function transfer_velocityᶠᶜᶜ(i, j, grid, time::Time, uₒ, vₒ,
                       uₐ::FourDArray, vₐ::FourDArray)

    Δu = @inbounds uₐ[i, j, 1, time] - uₒ[i, j, 1]
    Δv² = ℑyᵃᶜᵃ(i, j, 1, grid, Δϕt², vₐ, vₒ, time)
    return sqrt(Δu^2 + Δv²)
end

@inline function transfer_velocityᶜᶠᶜ(i, j, grid, time::Time, uₒ, vₒ,
                       uₐ::FourDArray, vₐ::FourDArray)

    Δu² = ℑxᶜᵃᵃ(i, j, 1, grid, Δϕt², uₐ, uₒ, time)
    Δv = @inbounds vₐ[i, j, 1, time] - vₒ[i, j, 1]
    return sqrt(Δu² + Δv^2)
end

@inline function Δu_transfer_velocity(i, j, grid, time::Time, uₒ, vₒ,
                       uₐ::FourDArray, vₐ::FourDArray)
    transfer_velocity = transfer_velocityᶠᶜᶜ(i, j, grid, time, uₒ, vₒ, uₐ, vₐ)
    Δu = @inbounds uₐ[i, j, 1, time] - uₒ[i, j, 1]
    return Δu * transfer_velocity
end

@inline function Δv_transfer_velocity(i, j, grid, time::Time, uₒ, vₒ,
                       uₐ::FourDArray, vₐ::FourDArray)
    transfer_velocity = transfer_velocityᶠᶜᶜ(i, j, grid, time, uₒ, vₒ, uₐ, vₐ)
    Δv = @inbounds vₐ[i, j, 1, time] - vₒ[i, j, 1]
    return Δv * transfer_velocity
end

@inline function x_atmosphere_ocean_momentum_flux(i, j, grid, clock, cᴰ, ρₒ, ρₐ, Uₒ, Uₐ)
    x_transfer_velocity_transfer_velocity = Δu_transfer_velocity(i, j, grid, time, Uₒ.u, Uₒ.v, Uₐ.u, Uₐ.v)
    return ρₐ / ρₒ * cᴰ * x_transfer_velocity_transfer_velocity
end 

@inline function y_atmosphere_ocean_momentum_flux(i, j, grid, clock, cᴰ, ρₒ, ρₐ, Uₒ, Uₐ)
    time = Time(clock.time)
    y_transfer_velocity_transfer_velocity = Δv_transfer_velocity(i, j, grid, time, Uₒ.u, Uₒ.v, Uₐ.u, Uₐ.v)
    return ρₐ / ρₒ * cᴰ * y_transfer_velocity_transfer_velocity
end

