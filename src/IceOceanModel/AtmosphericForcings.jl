module AtmosphericForcings

export PrescribedAtmosphere, PrescribedFluxes

using Adapt
using Oceananigans
using Oceananigans.Utils
using Oceananigans.BoundaryConditions: getbc
using KernelAbstractions: @kernel, @index

import IceOceanModel: compute_air_sea_fluxes!

abstract type AbstractAtmospericForcing end

# We generally have 2 types of atmospheric forcing: Prescribed fluxes and
# Prescribed atmospheric state (to treat with bulk formulae)
# This implementation also allows to have in future a prognostic atmospheric model

# Prescribed fluxes can be arrays, fields, of functions. 
# When functions, the signature should be 
# `f(i, j, grid, clock, fields)` where `fields` are the ocean model's prognostic fields
# in case of OnyOceanModel and the coupled model's prognostic fields in case of an `IceOceanModel`
# Parameters can be implemented using callable structs that subtype `Function`
struct PrescribedFluxes{T, S, U, V} <: AbstractAtmospericForcing
    heat_flux          :: T # heat flux
    freshwater_flux    :: S # freshwater flux
    zonal_stress       :: U # zonal stress
    meriodional_stress :: V # meriodional stress
end

Adapt.adapt_structure(to, f::PrescribedFluxes) = 
    PrescribedFluxes(Adapt.adapt(to, f.heat_flux),
                     Adapt.adapt(to, f.freshwater_flux),
                     Adapt.adapt(to, f.zonal_stress),
                     Adapt.adapt(to, f.meriodional_stress))

# Here we leverage a `getflux` function similar to the `getbc` from Oceananigans.jl to extract the fluxes,
# In this way we allow prescribed fluxes as well as relaxation fluxes
@kernel function _calculate_air_sea_fluxes!(Qˢ, Fˢ, τˣ, τʸ, ρₒ, cₒ, ε, grid, clock, fields, ice_thickness, solar_insolation, f::PrescribedFluxes)
    i, j = @index(Global, NTuple)
    @inbounds begin
        I₀ = solar_insolation[i, j, 1]
        Qˢ[i, j] = getflux(ice_thickness, f.heat_flux,          i, j, grid, clock, fields) + ε * I₀ / (ρₒ * cₒ)
        Fˢ[i, j] = getflux(ice_thickness, f.freshwater_fluxes,  i, j, grid, clock, fields)
        τˣ[i, j] = getflux(ice_thickness, f.zonal_stress,       i, j, grid, clock, fields)
        τʸ[i, j] = getflux(ice_thickness, f.meriodional_stress, i, j, grid, clock, fields)
    end
end

struct PrescribedAtmosphere{R, H, P, W, T, Q, D, C, G} <: AbstractAtmospericForcing
    adiabatic_lapse_rate :: R    # -
    atmosphere_state_height :: H # m
    reference_height :: H        # m 
    surface_pressure :: P        # Pa
    atmosphere_velocities :: W   # (m/s, m/s)
    air_temperature :: T         # deg ᵒC
    air_humidity :: Q            # kg/m³
    air_density :: D             # kg/m³
    cloud_cover_feedback :: C    # - 
    gamma_air :: C               # -
end

# The atmospheric state (T, u, v, q, ρ and p) can all be Values, Arrays, Functions, Fields or FieldTimeSerieses
function PrescribedAtmosphere(; 
                     adiabatic_lapse_rate    = 0.08,  
                     atmosphere_state_height = 10,       # m 
                     reference_height        = 10,       # m
                     surface_pressure        = 1e5,      # Pa
                     atmosphere_velocities,              # (m/s, m/s)
                     air_temperature,                    # deg ᵒC
                     air_humidity            = 0.01,     # kg/m³
                     air_density             = 1.25,     # kg/m³
                     cloud_cover_coeff       = 0.8,
                     gamma_air               = 0.01)

    return PrescribedAtmosphere(adiabatic_lapse_rate, 
                                atmosphere_state_height,
                                reference_height,
                                surface_pressure,
                                atmosphere_velocities,
                                air_temperature,
                                air_humidity,
                                air_density,
                                cloud_cover_coeff,
                                gamma_air)
end

Adapt.adapt_structure(to, f::PrescribedAtmosphere) = 
    PrescribedAtmosphere(Adapt.adapt(to, f.adiabatic_lapse_rate),  
                         Adapt.adapt(to, f.atmosphere_state_height),
                         Adapt.adapt(to, f.reference_height),
                         Adapt.adapt(to, f.surface_pressure),
                         Adapt.adapt(to, f.atmosphere_velocities),
                         Adapt.adapt(to, f.air_temperature),
                         Adapt.adapt(to, f.air_humidity),
                         Adapt.adapt(to, f.air_density),
                         Adapt.adapt(to, f.cloud_cover_feedback),
                         Adapt.adapt(to, f.gamma_air))

@inline clausius_clapeyron(FT, Tₛ) = convert(FT, 611.2) * exp(convert(FT, 17.67) * Tₛ / (Tₛ + convert(FT, 243.5)))

# Follows MITgcm
@kernel function _calculate_air_sea_fluxes!(Qˢ, Fˢ, τˣ, τʸ, ε, ρₒ, cₒ, grid, clock, fields, ice_thickness, f::PrescribedAtmosphere)
    
    hᵀ   = f.atmosphere_state_height
    α    = f.adiabatic_lapse_rate
    uˢ, vˢ = f.atmosphere_velocities

    Tₐ = getflux(f.air_temperature,  i, j, grid, clock, fields)
    uₐ = getflux(uˢ,                 i, j, grid, clock, fields)
    vₐ = getflux(vˢ,                 i, j, grid, clock, fields)
    qₐ = getflux(f.air_humidity,     i, j, grid, clock, fields)
    ρₐ = getflux(f.air_density,      i, j, grid, clock, fields)
    p₀ = getflux(f.surface_pressure, i, j, grid, clock, fields)
    
    h  = getflux(ice_thickness, i, j, grid, clock, fields)
    I₀ = solar_insolation[i, j, 1]

    s = sqrt(uₐ^2 + vₐ^2) # speed m / s
    γ = f.gamma_air
    
    # Physical constants
    σ = convert(FT, σᴮ) # W/m²/K⁴ Stefan-Boltzmann constant
    ℒ = convert(FT, ℒₑ) # J/kg Latent heat of evaporation
    
    FT = eltype(grid)

    # latent heat of evaporation
    ℒ = convert(FT, ℒₑ)

    Tₛ = fields.T[i, j, grid.Nz]
    T₀ = Tₐ*(1 - γ * qₐ)

    # sea-air temperature difference
    ΔT = T₀ - Tₛ + α*hᵀ

    # saturation vapour pressure (Pa)
    eₛ = clausius_clapeyron(FT, Tₛ)

    # saturation air humidity (kg/m³)
    qₛ = convert(FT, 0.622) * eₛ / (p₀ - eₛ)

    # air excess humidity
    Δq = qₐ - qₛ

    # Turbulent heat transfer coefficients
    Cᵀ, Cᵁ, Cq = turbulent_heat_transfer_coefficients(FT, f, T₀, qₐ, uₛ, ΔT, Δq)

    # sensible heat flux (W/m²)
    H = ρₐ * uₛ * Cᵀ * Cᵁ * ΔT

    # latent heat flux (W/m²)
    L = ρₐ * uₛ * Cq * Cᵁ * Δq * ℒ

    # net longwave radiation (W/m²)
    Rₙ = ε * σ * (Tₛ + convert(FT, 273.15))^4 * (1 - f.cloud_cover_feedback)
    
    @inbounds begin
        Qˢ[i, j, 1] = ifelse(ice_thickness(H + L + Rₙ + I₀ * ε) / (ρₒ * cₒ)
        # Fˢ[i, j, 1] = L / ℒ
        τˣ[i, j, 1] = ρₐ * uₛ * Cᵁ * uₐ / ρₒ 
        τʸ[i, j, 1] = ρₐ * uₛ * Cᵁ * vₐ / ρₒ 
    end

    return nothing
end

# Follows MITgcm 
@inline function turbulent_heat_transfer_coefficients(FT, f, T₀, qₐ, uₛ, ΔT, Δq)
    hᵀ = f.atmosphere_state_height
    zᴿ = f.reference_height
    λ  = log(hᵀ / zᴿ)

    Cᵀ = Cᵁ = Cq = convert(FT, 0.41) / log(zᴿ * 2)
    u★ = Cᵁ * uₛ
    T★ = Cᵀ * ΔT
    q★ = Cq * Δq

    @unroll for iter in 1:5
        G = Γ(FT, u★, T★, q★, T₀, qₐ, f)
        χ = sqrt(1 - 16 * G)

        ψˢ = ifelse(G > 0, -5G, 2 * log((1 + χ^2) / 2))
        ψᵐ = ifelse(G > 0, -5G, 2 * log((1 + χ) / 2) + ψˢ / 2 - 2 * atan(χ) + convert(FT, π/2))

        Cᵁ = Cᵁ / (1 + Cᵁ * (λ - ψᵐ) / convert(FT, 0.41))
        Cᵀ = Cᵀ / (1 + Cᵀ * (λ - ψˢ) / convert(FT, 0.41))
        u★ = Cᵁ * uₛ
        T★ = Cᵀ * ΔT
        q★ = Cq * Δq    
    end

    return Cᵀ, Cᵁ, Cq
end

@inline Γ(FT, u★, T★, q★, T₀, qₐ, f) = convert(FT, 0.41) * convert(FT, 9.80655) * f.reference_height / u★^2 * 
                                       (T★ / T₀ - q★ / (1/f.gamma_air - qₐ))

end