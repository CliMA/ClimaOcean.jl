using CUDA: @allowscalar
using Printf

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
    ρ = ℙᵢ.reference_density
    c = ℙᵢ.heat_capacity
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
    σ = ℙₛ.radiation.σ
    ϵ = ℙₛ.radiation.ϵ
    α = σ * ϵ
    Tₛ = (Tᵢ - h / k * (Qₐ + 4α * Tₛ⁻^4)) / (1 + 4α * h * Tₛ⁻^3 / k)
    =#

    Tₛ = Tᵢ - Qₐ * h / k

    # Under heating fluxes, cap surface temperature by melting temperature
    Tₘ = ℙᵢ.liquidus.freshwater_melting_temperature
    Tₘ = convert_to_kelvin(ℙᵢ.temperature_units, Tₘ)

    # Fix a NaN
    Tₛ = ifelse(isnan(Tₛ), Tₛ⁻, Tₛ)
    Tₛ = min(Tₛ, Tₘ)

    #=
    # Don't let it go below some minimum number?
    # FT = typeof(Tₛ⁻)
    # min_Tₛ = convert(FT, 230)
    # Tₛ = max(min_Tₛ, Tₛ)
    # 

    # Tₛ = (Tₛ + 9Tₛ⁻) / 10
    =#

    if Tₛ < -5
        @show Tₛ, Tₘ, Qₐ, h, k, Tᵢ
    end

    return Tₛ
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

    # TODO: this depends on the phase of the interface
    #ℰv = 0 #AtmosphericThermodynamics.latent_heat_vapor(ℂₐ, 𝒬ₐ)
    ℰs = AtmosphericThermodynamics.latent_heat_sublim(ℂₐ, 𝒬ₐ)

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
    Qc = - ρₐ * cₐ * u★ * θ★ # = - ρₐ cₐ u★ Ch / sqrt(Cd) * (θₐ - Tₛ)
    Qv = - ρₐ * ℰs * u★ * q★

    # Net heat flux
    Qa = Qr + Qc + Qv

    if Qa > 1000
        @show Qa, Qr, Qc, Qv, u★, θ★, q★, Tₛ⁻
    end

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
    melting :: Bool
end

InterfaceState(u★, θ★, q★, u, v, T, S, q) =
    InterfaceState(u★, θ★, q★, u, v, T, S, q, false)

Base.eltype(::InterfaceState{FT}) where FT = FT

function Base.show(io::IO, is::InterfaceState)
    print(io, "InterfaceState(",
          "u★=", prettysummary(is.u★), " ",
          "θ★=", prettysummary(is.θ★), " ",
          "q★=", prettysummary(is.q★), " ",
          "u=", prettysummary(is.u), " ",
          "v=", prettysummary(is.v), " ",
          "T=", prettysummary(is.T), " ",
          "S=", prettysummary(is.S), " ",
          "q=", prettysummary(is.q), ")")
end

@inline zero_interface_state(FT) = InterfaceState(zero(FT),
                                                  zero(FT),
                                                  zero(FT),
                                                  zero(FT),
                                                  zero(FT),
                                                  convert(FT, 273.15),
                                                  zero(FT),
                                                  zero(FT))

# Iterating condition for the characteristic scales solvers
@inline function iterating(Ψⁿ, Ψ⁻, iteration, maxiter, tolerance)
    hasnt_started = iteration == 0
    reached_maxiter = iteration ≥ maxiter
    drift = abs(Ψⁿ.u★ - Ψ⁻.u★) + abs(Ψⁿ.θ★ - Ψ⁻.θ★) + abs(Ψⁿ.q★ - Ψ⁻.q★)
    converged = drift < tolerance
    return !(converged | reached_maxiter) | hasnt_started
end

@inline function compute_interface_state(flux_formulation,
                                         initial_interface_state,
                                         atmosphere_state,
                                         interior_state,
                                         downwelling_radiation,
                                         interface_properties,
                                         atmosphere_properties,
                                         interior_properties)

    Ψₐ = atmosphere_state
    Ψᵢ = interior_state
    Ψₛⁿ = Ψₛ⁻ = initial_interface_state
    maxiter = flux_formulation.solver_maxiter
    tolerance = flux_formulation.solver_tolerance
    iteration = 0

    while iterating(Ψₛⁿ, Ψₛ⁻, iteration, maxiter, tolerance)
        Ψₛ⁻ = Ψₛⁿ
        Ψₛⁿ = iterate_interface_state(flux_formulation,
                                      Ψₛ⁻, Ψₐ, Ψᵢ,
                                      downwelling_radiation,
                                      interface_properties,
                                      atmosphere_properties,
                                      interior_properties)
        iteration += 1
    end

    return Ψₛⁿ
end

"""
    iterate_interface_state(flux_formulation, Ψₛⁿ⁻¹, Ψₐ, Ψᵢ, Qᵣ, ℙₛ, ℙₐ, ℙᵢ)

Return the nth iterate of the interface state `Ψₛⁿ` computed according to the
`flux_formulation`, given the interface state at the previous iterate `Ψₛⁿ⁻¹`,
as well as the atmosphere state `Ψₐ`, the interior state `Ψᵢ`,
downwelling radiation `Qᵣ`, and the interface, atmosphere,
and interior properties `ℙₛ`, `ℙₐ`, and `ℙᵢ`.
"""
@inline function iterate_interface_state(flux_formulation,
                                         approximate_interface_state,
                                         atmosphere_state,
                                         interior_state,
                                         downwelling_radiation,
                                         interface_properties,
                                         atmosphere_properties,
                                         interior_properties)
        
    Tₛ = compute_interface_temperature(interface_properties.temperature_formulation,
                                       approximate_interface_state,
                                       atmosphere_state,
                                       interior_state,
                                       downwelling_radiation,
                                       interface_properties,
                                       atmosphere_properties,
                                       interior_properties)

    # Thermodynamic state
    FT = eltype(approximate_interface_state)
    ℂₐ = atmosphere_properties.thermodynamics_parameters
    𝒬ₐ = atmosphere_state.𝒬
    ρₐ = 𝒬ₐ.ρ

    # Recompute the saturation specific humidity at the interface based on the new temperature
    q_formulation = interface_properties.specific_humidity_formulation
    Sₛ = approximate_interface_state.S
    qₛ = saturation_specific_humidity(q_formulation, ℂₐ, ρₐ, Tₛ, Sₛ)

    # Compute the specific humidity increment
    qₐ = AtmosphericThermodynamics.vapor_specific_humidity(ℂₐ, 𝒬ₐ)
    Δq = qₐ - qₛ

    # Temperature increment including the ``lapse rate'' `α = g / cₚ`
    zₐ = atmosphere_state.z
    zₛ = zero(FT)
    Δh = zₐ - zₛ
    Tₐ = AtmosphericThermodynamics.air_temperature(ℂₐ, 𝒬ₐ)
    g  = flux_formulation.gravitational_acceleration
    cₐ = AtmosphericThermodynamics.cp_m(ℂₐ, 𝒬ₐ)
    θₐ = Tₐ + g * Δh / cₐ
    Δθ = θₐ - Tₛ

    u★, θ★, q★ = iterate_interface_fluxes(flux_formulation,
                                          Tₛ, qₛ, Δθ, Δq, Δh,
                                          approximate_interface_state,
                                          atmosphere_state,
                                          atmosphere_properties)

    u = approximate_interface_state.u
    v = approximate_interface_state.v
    S = approximate_interface_state.S

    return InterfaceState(u★, θ★, q★, u, v, Tₛ, S, convert(FT, qₛ))
end

