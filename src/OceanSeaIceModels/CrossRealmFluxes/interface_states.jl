using CUDA: @allowscalar

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

# The flux balance could be solved either
# 
#   Tᵇ - Tₛⁿ⁺¹
# κ ---------- = Jᵀ (all fluxes positive upwards)
#       δ
#
# Where the LHS is the internal diffusive flux inside the ocean (within the boundary layer of thickness δ) 
# and the RHS are the atmospheric and radiative fluxes are provided explicitly, or
# 
#   Tᵇ - Tₛⁿ⁺¹    σ ϵ Tₛⁿ⁺¹Tₛⁿ³
# κ ---------- - ------------ = Jᵀ (all fluxes positive upwards)
#       δ           ρₒ cpₒ
#
# Where the LHS is the internal diffusive flux inside the ocean (within the boundary layer of thickness δ) 
# plus the (semi-implicit) outgoing longwave radiation and the RHS are the remaining atmospheric and radiative fluxes
# provided explicitly. Here we implement the fully explicit version, the linearized version is an optimization
# that can be explored in the future.
@inline flux_balance_temperature(F::DiffusiveFlux, Tᵇ, Jᵀ) = Tᵇ - Jᵀ / F.κ * F.δ

# the flaw here is that the ocean emissivity and albedo are fixed, but they might be a function of the
# interface temperature, so we might need to pass the radiation and the albedo and emissivity as arguments.
@inline function compute_interface_temperature(st::SkinTemperature, Tₛ, ℂ, 𝒬₀, 𝒬₁, Σ★, ρᵇ, cᵇ, Qd, σ, α, ϵ)
    ρₐ = AtmosphericThermodynamics.air_density(ℂ, 𝒬₁)
    cₚ = AtmosphericThermodynamics.cp_m(ℂ, 𝒬₁) # moist heat capacity
    ℰv = AtmosphericThermodynamics.latent_heat_vapor(ℂ, 𝒬₁)

    # upwelling radiation is calculated explicitly 
    Qu = upwelling_radiation(Tₛ, σ, ϵ)
    Qn = Qd + Qu # Net radiation (positive out of the ocean)

    u★ = Σ★.momentum
    T★ = Σ★.temperature
    q★ = Σ★.water_vapor
 
    # Turbulent heat fluxes, sensible + latent (positive out of the ocean)
    Qt = - ρₐ * u★ * (cₚ * T★ + q★ * ℰv)

    # Net temperature flux (positive upwards)
    Jᵀ = (Qt + Qn) / (ρᵇ * cᵇ)

    Tₒ = AtmosphericThermodynamics.air_temperature(ℂ, 𝒬₀)
    Tₛ = flux_balance_temperature(st.internal_flux, Tₒ, Jᵀ) # new interface temperature

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

