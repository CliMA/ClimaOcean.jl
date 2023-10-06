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
@kernel function _calculate_prescribed_fluxes!(Qˢ, Fˢ, τˣ, τʸ, fields, f::PrescribedFluxes)
    i, j = @index(Global, NTuple)
    @inbounds begin
        Qˢ[i, j] = getflux(f.heat_flux,          i, j, grid, clock, fields)
        Fˢ[i, j] = getflux(f.freshwater_fluxes,  i, j, grid, clock, fields)
        τˣ[i, j] = getflux(f.zonal_stress,       i, j, grid, clock, fields)
        τʸ[i, j] = getflux(f.meriodional_stress, i, j, grid, clock, fields)
    end
end


struct PrescribedAtmosphere{R, H, P, W, T, Q, D, C, G} <: AbstractAtmospericForcing
    adiabatic_lapse_rate :: R    # -
    atmosphere_state_height :: H # m
    reference_height :: H        # m 
    surface_pressure :: P        # Pa
    atmosphere_velocity :: W     # (m/s, m/s)
    air_temperature :: T         # deg ᵒC
    air_humidity :: Q            # kg/m³
    air_density :: D             # kg/m³
    cloud_cover_feedback :: C    # - 
    gamma_air :: C               # -
end

const PrescribedAtmosphereModel = IceOceanModel{<:Any, <:Any, PrescribedAtmosphere}

# To put in ClimaOcean.jl integrating with ClimaSeaIce.jl (To modify)
#=
using Adapt 
using KernelAbstractions.Extras.LoopInfo: @unroll

using Oceananigans.Utils: Time

struct AirSeaFlux{A, R, H, P, W, T, Q, D, C, E} <: Function
    adiabatic_lapse_rate :: R    # -
    atmosphere_state_height :: H # m
    reference_height :: H        # m 
    surface_pressure :: P        # Pa
    wind_speed :: W              # m/s
    air_temperature :: T         # deg ᵒC
    air_humidity :: Q            # kg/m³
    air_density :: D             # kg/m³
    cloud_cover_feedback :: C    # - 
    gamma_air :: C               # -
    ocean_emissivity :: E        # -

    AirSeaFlux{A}(adiabatic_lapse_rate::R, 
                  atmosphere_state_height::H,
                  reference_height::H,
                  surface_pressure::P,
                  wind_speed::W,
                  air_temperature::T,
                  air_humidity::Q,
                  air_density::D,
                  cloud_cover_coeff::C,
                  gamma_air::C,
                  ocean_emissivity::E) where {A, R, H, P, W, T, Q, D, C, E} = 
                  new{A, R, H, P, W, T, Q, D, C, E}(adiabatic_lapse_rate, 
                                                    atmosphere_state_height,
                                                    reference_height,
                                                    surface_pressure,
                                                    wind_speed,
                                                    air_temperature,
                                                    air_humidity,
                                                    air_density,
                                                    cloud_cover_coeff,
                                                    gamma_air,
                                                    ocean_emissivity)
end

struct HeatFlux end
struct WindStress end

function AirSeaFlux(; 
                     adiabatic_lapse_rate    = 0.08,  
                     atmosphere_state_height = 10,       # m 
                     reference_height        = 10,       # m
                     surface_pressure        = 1e5,      # Pa
                     wind_speed,                         # m/s
                     air_temperature,                    # deg ᵒC
                     air_humidity            = 0.01,     # kg/m³
                     air_density             = 1.25,     # kg/m³
                     cloud_cover_coeff       = 0.8,
                     gamma_air               = 0.01,
                     ocean_emissivity        = 0.9, 
                     flux_or_stress::A       = HeatFlux()) where A

    return AirSeaFlux{A}(adiabatic_lapse_rate, 
                         atmosphere_state_height,
                         reference_height,
                         surface_pressure,
                         wind_speed,
                         air_temperature,
                         air_humidity,
                         air_density,
                         cloud_cover_coeff,
                         gamma_air,
                         ocean_emissivity)
end

Adapt.adapt_structure(to, f::AirSeaFlux{A}) where A = 
                AirSeaFlux{A}(Adapt.adapt(to, f.adiabatic_lapse_rate),  
                              Adapt.adapt(to, f.atmosphere_state_height),
                              Adapt.adapt(to, f.reference_height),
                              Adapt.adapt(to, f.surface_pressure),
                              Adapt.adapt(to, f.wind_speed),
                              Adapt.adapt(to, f.air_temperature),
                              Adapt.adapt(to, f.air_humidity),
                              Adapt.adapt(to, f.air_density),
                              Adapt.adapt(to, f.cloud_cover_feedback),
                              Adapt.adapt(to, f.gamma_air),
                              Adapt.adapt(to, f.ocean_emissivity))

const AirSeaHeatFlux   = AirSeaFlux{HeatFlux}
const AirSeaWindStress = AirSeaFlux{WindStress}

function (f::AirSeaHeatFlux)(i, j, grid, clock, fields)
    
    hᵀ = f.atmosphere_state_height
    α  = f.adiabatic_lapse_rate
    Tₐ = f.air_temperature[i, j, 1, Time(clock.time)]
    uₛ = f.wind_speed[i, j, 1, Time(clock.time)]
    qₐ = f.air_humidity
    ρₐ = f.air_density
    γ  = f.gamma_air
    p₀ = f.surface_pressure

    FT = eltype(grid)

    # latent heat of evaporation
    ℒ = convert(FT, 2.5e6)

    Tₛ = fields.T[i, j, grid.Nz]
    T₀ = Tₐ*(1 - γ * qₐ)

    # sea-air temperature difference
    ΔT = T₀ - Tₛ + α*hᵀ

    # saturation vapour pressure (Pa)
    eₛ = convert(FT, 611.2) * exp(convert(FT, 17.67) * Tₛ / (Tₛ + convert(FT, 243.5)))

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
    Rₙ = f.ocean_emissivity * convert(FT, 5.67e-8) * (Tₛ + convert(FT, 273.15))^4 * (1 - f.cloud_cover_feedback)

    return H + L + Rₙ
end

function (f::AirSeaWindStress)(i, j, grid, clock, fields)
    
    uₛ = f.wind_speed[i, j, 1, Time(clock.time)]
    ρₐ = f.air_density
    FT = eltype(grid)
    
    # Drag coefficient 
    Cᴰ = ifelse(uₛ > 25, convert(FT, 2.3e-3), convert(FT, 1.3e-3))

    return ρₐ * Cᴰ * uₛ^2 
end

# Follows MITgcm 
@inline function turbulent_heat_transfer_coefficients(FT, f, T₀, qₐ, uₛ, ΔT, Δq)
    hᵀ = f.atmosphere_state_height
    zᴿ = f.reference_height
    λ  = log(hᵀ / zᴿ)
    γ  = f.gamma_air

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

AirSeaHeatFluxBoundaryCondition(; kwargs...) = 
    FluxBoundaryCondition(AirSeaFlux(; flux_or_stress = HeatFlux(), kwargs...), discrete_form = true)

AirSeaWindStressBoundaryCondition(; kwargs...) = 
    FluxBoundaryCondition(AirSeaFlux(; flux_or_stress = WindStress(), kwargs...), discrete_form = true)


=#


end