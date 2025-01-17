using CUDA: @allowscalar

import ClimaSeaIce
import Thermodynamics as AtmosphericThermodynamics  
using Thermodynamics: Liquid, Ice

#####
##### Interface properties
#####

struct InterfaceProperties{R, Q, T}
    radiation :: R
    specific_humidity_formulation :: Q
    temperature_formulation :: T
end

#####
##### Interface specific humidity formulations
#####

# TODO: allow different saturation models
# struct ClasiusClapyeronSaturation end
struct SpecificHumidityFormulation{Φ, X}
    # saturation :: S 
    phase :: Φ
    water_mole_fraction :: X
end

"""
    SpecificHumidityFormulation(phase [, water_mole_fraction=1])

Return the formulation for computing specific humidity at an interface.
"""
SpecificHumidityFormulation(phase) = SpecificHumidityFormulation(phase, nothing)

@inline compute_water_mole_fraction(::Nothing, salinity) = 1
@inline compute_water_mole_fraction(x_H₂O::Number, salinity) = x_H₂O

@inline function saturation_specific_humidity(formulation::SpecificHumidityFormulation, ℂₐ, ρₛ, Tₛ, Sₛ=zero(Tₛ))
    x_H₂O = compute_water_mole_fraction(formulation.water_mole_fraction, Sₛ)
    phase = formulation.phase

    CT = eltype(ℂₐ)
    p★ = Thermodynamics.saturation_vapor_pressure(ℂₐ, convert(CT, Tₛ), phase)
    q★ = Thermodynamics.q_vap_saturation_from_density(ℂₐ, convert(CT, Tₛ), convert(CT, ρₛ), p★)

    # Compute saturation specific humidity according to Raoult's law
    return q★ * x_H₂O
end

struct SalinityConstituent{FT}
    molar_mass :: FT
    mass_fraction :: FT
end

struct WaterMoleFraction{FT, C}
    water_molar_mass :: FT
    salinity_constituents :: C
end

function WaterMoleFraction(FT=Float64)
    water_molar_mass = convert(FT, 18.02)

    # TODO: find reference for these
    salinity_constituents = (
        chloride  = SalinityConstituent{FT}(35.45, 0.56),
        sodium    = SalinityConstituent{FT}(22.99, 0.31),
        sulfate   = SalinityConstituent{FT}(96.06, 0.08),
        magnesium = SalinityConstituent{FT}(24.31, 0.05),
    )

    return SeawaterComposition(water_molar_mass, salinity_constituents)
end

@inline function compute_water_mole_fraction(wmf::WaterMoleFraction, S)
    # TODO: express the concept of "ocean_salinity_units"?
    s = S / 1000 # convert g/kg to concentration

    # Molecular weights
    μ_H₂O = wmf.water_molar_mass

    # Salinity constituents: Cl⁻, Na, SO₄, Mg
    μ_Cl  = wmf.salinity_constituents.chloride.molar_mass
    μ_Na  = wmf.salinity_constituents.sodium.molar_mass
    μ_SO₄ = wmf.salinity_constituents.sulfate.molar_mass
    μ_Mg  = wmf.salinity_constituents.magnesium.molar_mass

    # Salinity constituent fractions
    ϵ_Cl  = wmf.salinity_constituents.chloride.mass_fraction
    ϵ_Na  = wmf.salinity_constituents.sodium.mass_fraction
    ϵ_SO₄ = wmf.salinity_constituents.sulfate.mass_fraction
    ϵ_Mg  = wmf.salinity_constituents.magnesium.mass_fraction

    α = μ_H₂O * (ϵ_Cl/μ_Cl + ϵ_Na/μ_Na  + ϵ_SO₄/μ_SO₄ + ϵ_Mg/μ_Mg)

    return (1 - s) / (1 - s + α * s)
end

####
#### Interface temperature formulations
####

"""
    struct BulkTemperature

A type to represent the interface temperature used in fixed-point iteration for interface
fluxes following similarity theory. The interface temperature is not calculated but instead
provided by either the ocean or the sea ice model.
"""
struct BulkTemperature end

# Do nothing (just copy the temperature)
@inline compute_interface_temperature(::BulkTemperature, Ψₛ, args...) = Ψₛ.T

####
#### Skin interface temperature calculated as a flux balance
####

""" 
    struct SkinTemperature     

A type to represent the interface temperature used in the flux calculation.
The interface temperature is calculated from the flux balance at the interface.
In particular, the interface temperature ``Tₛ`` is the root of:
 
```math
F(Tₛ) - Jᵀ = 0
```

where ``Jᵀ`` are the fluxes at the top of the interface (turbulent + radiative), and
``F`` is the internal diffusive flux dependent on the interface temperature itself.

Note that all fluxes positive upwards.
"""
struct SkinTemperature{I}
    internal_flux :: I
end

struct DiffusiveFlux{Z, K}
    δ :: Z # Boundary layer thickness, as a first guess we will use half the grid spacing
    κ :: K # diffusivity in m² s⁻¹
end

# A default constructor for SkinTemperature
function SkinTemperature(FT::DataType=Float64; κ=1e-2, δ=1.0) 
    internal_flux = DiffusiveFlux(FT; κ, δ)
    return SkinTemperature(internal_flux)
end

DiffusiveFlux(FT; κ = 1e-2, δ = 1.0) = DiffusiveFlux(convert(FT, δ), convert(FT, κ))

# The flux balance is solved by computing
# 
#            κ 
# Jᵃ(Tₛⁿ) + --- (Tₛⁿ⁺¹ - Tᵢ) = 0
#            δ
#
# where Jᵃ is the external flux impinging on the surface from above and
# Jᵢ = - κ (Tₛ - Tᵢ) / δ is the "internal flux" coming up from below.
# We have indicated that Jᵃ may depend on the surface temperature from the previous
# iterate. We thus find that 
#
# Tₛⁿ⁺¹ = Tᵢ + δ * Jᵃ(Tₛⁿ) / κ
#
# Note that we could also use the fact that Jᵃ(T) = σ * ϵ * T^4 + ⋯
# to expand Jᵃ around Tⁿ⁺¹,
#
# Jᵃ(Tⁿ⁺¹) ≈ Jᵃ(Tⁿ) + (Tⁿ⁺¹ - Tⁿ) * ∂T_Jᵃ(Tⁿ)
#          ≈ Jᵃ(Tⁿ) + 4 * (Tⁿ⁺¹ - Tⁿ) σ * ϵ * Tⁿ^3 / (ρ c)
#
# which produces the alternative, semi-implicit flux balance
# 
#                                      κ 
# Jᵃ(Tₛⁿ) - 4 α Tₛⁿ⁴ + 4 α Tₛⁿ Tₛⁿ³ + --- (Tₛⁿ⁺¹ - Tᵢ) = 0
#                                      δ
#
# with α = σ ϵ / (ρ c) such that
#
# Tₛⁿ⁺¹ (κ / δ + 4 α Tₛⁿ³) = κ * Tᵢ / δ - Jᵃ + 4 α Tₛⁿ⁴)
#
# or
#
# Tₛⁿ⁺¹ = = (Tᵢ - δ / κ * (Jᵃ - 4 α Tₛⁿ⁴)) / (1 + 4 δ σ ϵ Tₛⁿ³ / ρ c κ) 
#
# corresponding to a linearization of the outgoing longwave radiation term.
@inline function flux_balance_temperature(F::DiffusiveFlux, Qₐ, Ψₛ, ℙₛ, Ψᵢ, ℙᵢ)
    ρ = ℙ.reference_density
    c = ℙ.heat_capacity
    Jᵀ = Qₐ / (ρ * c)
    return Ψᵢ.T + Jᵀ * F.δ / F.κ
end

# Q + k / h * (Tˢ - Tᵢ) = 0
# ⟹  Tₛ = Tᵢ - Q * h / k
@inline function flux_balance_temperature(F::ClimaSeaIce.ConductiveFlux, Qₐ, Ψₛ, ℙₛ, Ψᵢ, ℙᵢ)
    k = F.conductivity
    h = Ψᵢ.h

    # Bottom temperature at the melting temperature
    Tᵢ = ClimaSeaIce.SeaIceThermodynamics.melting_temperature(ℙᵢ.liquidus, Ψᵢ.S)
    Tᵢ = convert_to_kelvin(ℙᵢ.temperature_units, Tᵢ)
    Tₛ⁻ = Ψₛ.T

    #=
    @show Tᵢ Tₛ⁻
    @show Qₐ
    @show h
    @show k
    @show Qₐ * h / k
    =#

    k *= 100
    Tₛ = Tᵢ - Qₐ * h / k

    #=
    σ = ℙₛ.radiation.σ
    ϵ = ℙₛ.radiation.ϵ
    α = σ * ϵ
    Tₛ = (Tᵢ - h / k * (Qₐ + 4α * Tₛ⁻^4)) / (1 + 4α * h * Tₛ⁻^3 / k)
    Tₛ = ifelse(isnan(Tₛ), Tₛ⁻, Tₛ)
    =#

    # @show Tₛ

    # Under heating fluxes, cap surface temperature by melting temperature
    Tₘ = ℙᵢ.liquidus.freshwater_melting_temperature
    Tₘ = convert_to_kelvin(ℙᵢ.temperature_units, Tₘ)

    # Don't let it go below 0?
    Tₛ = max(zero(Tₛ), Tₛ)

    return min(Tₛ, Tₘ)
end

@inline function compute_interface_temperature(st::SkinTemperature,
                                               interface_state,
                                               atmosphere_state,
                                               interior_state,
                                               downwelling_radiation,
                                               interface_properties,
                                               atmosphere_properties,
                                               interior_properties)

    ℂₐ = atmosphere_properties.thermodynamics_parameters
    𝒬ₐ = atmosphere_state.𝒬
    ρₐ = AtmosphericThermodynamics.air_density(ℂₐ, 𝒬ₐ)
    cₐ = AtmosphericThermodynamics.cp_m(ℂₐ, 𝒬ₐ) # moist heat capacity
    ℰv = AtmosphericThermodynamics.latent_heat_vapor(ℂₐ, 𝒬ₐ)

    # upwelling radiation is calculated explicitly 
    Tₛ⁻ = interface_state.T # approximate interface temperature from previous iteration
    σ = interface_properties.radiation.σ
    ϵ = interface_properties.radiation.ϵ
    α = interface_properties.radiation.α

    Qu = upwelling_radiation(Tₛ⁻, σ, ϵ)
    Qd = net_downwelling_radiation(downwelling_radiation, α, ϵ)
    Qr = Qd + Qu # Net radiation (positive out of the ocean)

    u★ = interface_state.u★
    θ★ = interface_state.θ★
    q★ = interface_state.q★
 
    # Turbulent heat fluxes, sensible + latent (positive out of the ocean)
    Qc = - ρₐ * cₐ * u★ * θ★
    Qv = - ρₐ * ℰv * u★ * q★

    # Net heat flux
    Qa = Qr + Qc + Qv

    Tₛ = flux_balance_temperature(st.internal_flux, Qa,
                                  interface_state,
                                  interface_properties,
                                  interior_state,
                                  interior_properties)

    return Tₛ
end

######
###### Interface state
######

struct InterfaceState{FT}
    u★ :: FT # friction velocity
    θ★ :: FT # flux characteristic temperature
    q★ :: FT # flux characteristic specific humidity
    u :: FT  # interface x-velocity
    v :: FT  # interface y-velocity
    T :: FT  # interface temperature
    S :: FT  # interface salinity
    q :: FT  # interface specific humidity
end

Base.eltype(::InterfaceState{FT}) where FT = FT

zero_interface_state(FT) = InterfaceState(zero(FT),
                                          zero(FT),
                                          zero(FT),
                                          zero(FT),
                                          zero(FT),
                                          convert(FT, 273.15),
                                          zero(FT),
                                          zero(FT))

