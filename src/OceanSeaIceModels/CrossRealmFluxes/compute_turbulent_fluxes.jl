using Printf

import Thermodynamics as AtmosphericThermodynamics

using Thermodynamics: PhasePartition
using KernelAbstractions.Extras.LoopInfo: @unroll

@inline update_turbulent_flux_fields!(::Nothing, args...) = nothing

@inline function update_turbulent_flux_fields!(fields, i, j, grid, fluxes)
    Qv = fields.latent_heat
    Qc = fields.sensible_heat
    Fv = fields.water_vapor
    τx = fields.x_momentum
    τy = fields.y_momentum
    kᴺ = size(grid, 3) # index of the top ocean cell
    inactive = inactive_node(i, j, kᴺ, grid, c, c, c)
    @inbounds begin
        # +0: cooling, -0: heating
        Qv[i, j, 1] = ifelse(inactive, 0, fluxes.latent_heat)
        Qc[i, j, 1] = ifelse(inactive, 0, fluxes.sensible_heat)
        Fv[i, j, 1] = ifelse(inactive, 0, fluxes.water_vapor)
        τx[i, j, 1] = ifelse(inactive, 0, fluxes.x_momentum)
        τy[i, j, 1] = ifelse(inactive, 0, fluxes.y_momentum)
    end
    return nothing
end

@inline compute_turbulent_fluxes(turbulent_fluxes, atmos_state, ocean_state) =
    compute_turbulent_fluxes(turbulent_fluxes.roughness_lengths, turbulent_fluxes, atmos_state, ocean_state)

#####
##### Struct that represents a 3-tuple of momentum, heat, and water vapor
#####

struct SimilarityScales{U, T, Q}
    momentum :: U
    temperature :: T
    water_vapor :: Q
end

# Convenience default with water_vapor component = nothing
SimilarityScales(momentum, temperature) = SimilarityScales(momentum, temperature, nothing)

#####
##### Interface into SurfaceFluxes.jl
#####

# This is the case that SurfaceFluxes.jl can do
const NothingVaporRoughnessLength = SimilarityScales{<:Number, <:Number, Nothing}

@inline function compute_turbulent_fluxes(roughness_lengths::NothingVaporRoughnessLength,
                                          turbulent_fluxes,
                                          atmos_state,
                                          ocean_state)

    # Constant roughness lengths
    ℓu = roughness_lengths.momentum
    ℓθ = roughness_lengths.temperature

    # Solve for the surface fluxes with initial roughness length guess
    Uᵍ = zero(zm) # gustiness
    β = one(zm)   # surface "resistance"
    values = SurfaceFluxes.ValuesOnly(atmos_state, ocean_state, ℓu, ℓθ, Uᵍ, β)
    conditions = SurfaceFluxes.surface_conditions(turbulent_fluxes, values)

    fluxes = (;
        sensible_heat = conditions.shf,
        latent_heat = conditions.lhf,
        water_vapor = conditions.evaporation,
        x_momentum = conditions.ρτxz,
        y_momentum = conditions.ρτyz,
    )

    return fluxes
end

#####
##### Fixed-point iteration for roughness length
#####

const ConstantRoughnessLength = SimilarityScales{<:Number, <:Number, <:Number}

struct SimilarityFunction{FT, C}
    a :: FT
    b :: FT
    c :: C
end

@inline function (ψ::SimilarityFunction)(Ri)
    a = ψ.a
    b = ψ.b
    c = ψ.c

    Ri⁻ = min(zero(Ri), Ri)
    ϕ⁻¹ = (1 - b * Ri⁻)^c
    ψ_unstable = log((1 + ϕ⁻¹)^2 * (1 + ϕ⁻¹^2) / 8) - (4 * atan(ϕ⁻¹) + π) / 2

    ψ_stable = - a * Ri

    return ifelse(Ri < 0, ψ_unstable, ψ_stable)
end

struct OneQuarter end
struct OneHalf end

import Base: ^
@inline ^(x, ::OneQuarter) = sqrt(sqrt(x))
@inline ^(x, ::OneHalf) = sqrt(x)

function businger_similarity_functions(FT=Float64)
    au = convert(FT, 4.7)
    bu = convert(FT, 15)
    cu = OneQuarter()
    ψu = SimilarityFunction(au, bu, cu)

    ah = convert(FT, 6.35)
    bh = convert(FT, 9)
    ch = OneHalf()
    ψh = SimilarityFunction(ah, bh, ch)

    ψq = ψh

    return SimilarityScales(ψu, ψh, ψq)
end

@inline function bulk_factor(ψ, h, ℓ, Ri)
    L★ = h / Ri
    χ⁻¹ = log(h / ℓ) - ψ(Ri) + ψ(ℓ / L★)
    return 1 / χ⁻¹
end

@inline function buoyancy_scale(θ★, q★, 𝒬, parameters)
    ℂ = parameters.thermodynamics_parameters

    𝒯₀ = AtmosphericThermodynamics.virtual_temperature(ℂ, 𝒬)
    θ₀ = AtmosphericThermodynamics.air_temperature(ℂ, 𝒬)
    q₀ = AtmosphericThermodynamics.vapor_specific_humidity(ℂ, 𝒬)

    ε = AtmosphericThermodynamics.Parameters.molmass_ratio(parameters)
    δ = ε - 1
    g = SurfaceFluxes.Parameters.grav(parameters)

    b★ = g / 𝒯₀ * (θ★ * (1 + δ * q₀) + δ * θ₀ * q★)

    return b★
end


@inline function state_differences(ℂ, 𝒰₁, 𝒰₀)
    z₁ = 𝒰₁.z
    z₀ = 𝒰₀.z
    Δh = z₁ - z₀

    U₁ = 𝒰₁.u
    U₀ = 𝒰₀.u

    @inbounds begin
        Δu = U₁[1] - U₀[1]
        Δv = U₁[2] - U₀[2]
    end

    # Thermodynamic state
    𝒬₁ = 𝒰₁.ts
    𝒬₀ = 𝒰₀.ts

    θ₁ = AtmosphericThermodynamics.air_temperature(ℂ, 𝒬₁)
    θ₀ = AtmosphericThermodynamics.air_temperature(ℂ, 𝒬₀)
    Δθ = θ₁ - θ₀

    q₁ = AtmosphericThermodynamics.vapor_specific_humidity(ℂ, 𝒬₁)
    q₀ = AtmosphericThermodynamics.vapor_specific_humidity(ℂ, 𝒬₀)
    Δq = q₁ - q₀

    return Δh, Δu, Δv, Δθ, Δq
end

@inline function compute_turbulent_fluxes(roughness_lengths::ConstantRoughnessLength,
                                          turbulent_fluxes,
                                          atmos_state,
                                          ocean_state)

    # Prescribed difference between two states
    ℂₐ = thermodynamics_params(turbulent_fluxes)
    Δh, Δu, Δv, Δθ, Δq = state_differences(ℂₐ, atmos_state, ocean_state)
    differences = (; u=Δu, v=Δv, θ=Δθ, q=Δq, h=Δh)

    # Solve for the characteristic scales u★, θ★, q★, and thus for fluxes.
    Γ₀ = Γ★ = SimilarityScales(1e-3, 1e-3, 1e-3)

    @unroll for iter = 1:10
        Γ★ = refine_characteristic_scales(Γ★,
                                          roughness_lengths, 
                                          turbulent_fluxes.similarity_functions,
                                          differences,
                                          ocean_state,
                                          turbulent_fluxes)

        @debug begin
            u★ = Γ★.momentum
            θ★ = Γ★.temperature
            q★ = Γ★.water_vapor
            # u★² = Cᴰ * (Δu² + Δv²)
            Cᴰ = u★^2 / (Δu^2 + Δv^2)
            @sprintf("Iter: %d, Cᴰ: %.4e, u★: %.4e, θ★: %.4e, q★: %.4e", iter, Cᴰ, u★ , θ★, q★)
        end
    end

    u★ = Γ★.momentum
    θ★ = Γ★.temperature
    q★ = Γ★.water_vapor

    # u★² ≡ sqrt(τx² + τy²)
    τx = u★^2 * Δu / sqrt(Δu^2 + Δv^2)
    τy = u★^2 * Δv / sqrt(Δu^2 + Δv^2)

    𝒬ₐ = atmos_state.ts
    ρₐ = AtmosphericThermodynamics.air_density(ℂₐ, 𝒬ₐ)
    cₚ = AtmosphericThermodynamics.cp_m(ℂₐ, 𝒬ₐ) # moist heat capacity
    ℰv = AtmosphericThermodynamics.latent_heat_vapor(ℂₐ, 𝒬ₐ)

    fluxes = (;
        water_vapor   = ρₐ * u★ * q★,
        sensible_heat = ρₐ * cₚ * u★ * θ★,
        latent_heat   = ρₐ * u★ * q★ * ℰv,
        x_momentum    = ρₐ * τx,
        y_momentum    = ρₐ * τy,
    )

    return fluxes
end

@inline function refine_characteristic_scales(estimated_characteristic_scales,
                                              roughness_lengths,
                                              similarity_functions,
                                              differences,
                                              ocean_state,
                                              similarity_parameters)

    # "initial" scales because we will recompute them
    u★ = estimated_characteristic_scales.momentum
    θ★ = estimated_characteristic_scales.temperature
    q★ = estimated_characteristic_scales.water_vapor

    # Extract roughness lengths
    ℓu = roughness_lengths.momentum
    ℓθ = roughness_lengths.temperature
    ℓq = roughness_lengths.water_vapor

    # Compute flux Richardson number
    h = differences.h
    ϰ = similarity_parameters.von_karman_constant

    𝒬ₒ = ocean_state.ts # thermodyanmic state
    b★ = buoyancy_scale(θ★, q★, 𝒬ₒ, similarity_parameters)
    Riₕ = - ϰ * h * b★ / u★^2

    # Compute similarity functions
    fu = similarity_functions.momentum
    fθ = similarity_functions.temperature
    fq = similarity_functions.water_vapor

    χu = bulk_factor(fu, h, ℓu, Riₕ)
    χθ = bulk_factor(fθ, h, ℓθ, Riₕ)
    χq = bulk_factor(fq, h, ℓq, Riₕ)

    Δu = differences.u
    Δv = differences.v
    Δθ = differences.θ
    Δq = differences.q

    u★ = ϰ * χu * sqrt(Δu^2 + Δv^2)
    θ★ = ϰ * χθ * Δθ
    q★ = ϰ * χq * Δq

    return SimilarityScales(u★, θ★, q★)
end

struct GravityWaveRoughnessLength{FT}
    gravitational_acceleration :: FT
    gravity_wave_parameter :: FT
    laminar_parameter :: FT
end

