struct AtmosphereOceanMomentumFlux{F, CD}
    formulation :: F
    drag_coefficient :: CD
end

#####
##### Relative velocities
#####

# Compute the square of a FieldTimeSeries at `time`
@inline Δϕt²(i, j, k, grid, ϕ1t, ϕ2, time) = @inbounds (ϕ1t[i, j, k, time] - ϕ2[i, j, k])^2

const ThreeDArray = AbstractArray{<:Any, 3}
const FourDArray = AbstractArray{<:Any, 4}

@inline function Δu_mod_ΔU(i, j, grid, time::Time,
                           uₒ::ThreeDArray,
                           vₒ::ThreeDArray,
                           uₐ::FourDArray,
                           vₐ::FourDArray)

    Δu = @inbounds uₐ[i, j, 1, time] - uₒ[i, j, 1]
    Δv² = ℑyᶠᶜᶜ(i, j, 1, grid, Δϕt², vₐ, vₒ, time)
    return Δu * sqrt(Δu^2 + Δv²)
end

@inline function Δv_mod_ΔU(i, j, grid, time::Time,
                           uₒ::ThreeDArray,
                           vₒ::ThreeDArray,
                           uₐ::FourDArray,
                           vₐ::FourDArray)

    Δu² = ℑxᶜᶠᶜ(i, j, 1, grid, Δϕt², uₐ, uₒ, time)
    Δv = @inbounds vₐ[i, j, 1, time] - vₒ[i, j, 1]
    return Δv * sqrt(Δu² + Δv^2)
end

@inline function x_atmopsphere_ocean_momentum_flux(i, j, grid, clock, cᴰ, ρₒ, ρₐ, Uₒ, Uₐ)
    time = Time(clock.time)
    x_ΔU_ΔU = Δu_mod_ΔU(i, j, grid, time, Uₒ.u, Uₒ.v, Uₐ.u, Uₐ.v)
    return ρₐ / ρₒ * cᴰ * x_ΔU_ΔU
end 

@inline function y_atmopsphere_ocean_momentum_flux(i, j, grid, clock, cᴰ, ρₒ, ρₐ, Uₒ, Uₐ)
    time = Time(clock.time)
    y_ΔU_ΔU = Δv_mod_ΔU(i, j, grid, time, Uₒ.u, Uₒ.v, Uₐ.u, Uₐ.v)
    return ρₐ / ρₒ * cᴰ * y_ΔU_ΔU
end

