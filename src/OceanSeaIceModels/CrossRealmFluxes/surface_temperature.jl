using CUDA: @allowscalar
import Thermodynamics as AtmosphericThermodynamics  

####
#### Utilities
####

@inline upwelling_radiation(θ₀, ::Nothing) = zero(θ₀)
@inline upwelling_radiation(θ₀, r) = r.σ * r.ϵ * θ₀^4

# For any surface temperture type that does not depend on the grid
regularize_surface_temperature(surface_temperature_type, grid) = surface_temperature_type

####
#### Prescribed surface temperature (the easiest case)
####

struct PrescribedSurfaceTemperature end

# Do nothing (just copy the temperature)
@inline retrieve_temperature(::PrescribedSurfaceTemperature, θ₀, args...) = θ₀

####
#### Diagnostic surface temperature calculated as a flux balance
####

struct DiagnosticSurfaceTemperature{I}
    internal_flux :: I
end

struct DiffusiveFlux{Z, K}
    δ :: Z # Boundary layer thickness, as a first guess we will use half the grid spacing
    κ :: K # diffusivity in m²s⁻¹
end

# A default constructor for DiagnosticSurfaceTemperature
function DiagnosticSurfaceTemperature(; κ = 0.1) 
    internal_flux = DiffusiveFlux(; κ)
    return DiagnosticSurfaceTemperature(internal_flux)
end

DiffusiveFlux(; κ = 1e-2) = DiffusiveFlux(nothing, κ)

function DiffusiveFlux(grid; κ = 1e-2) 
    Δz = @allowscalar Δzᶜᶜᶜ(1, 1, grid.Nz, grid)
    return DiffusiveFlux(Δz, κ)
end

regularize_surface_temperature(T::DiagnosticSurfaceTemperature{<:DiffusiveFlux}, grid) =
    DiagnosticSurfaceTemperature(DiffusiveFlux(grid; κ = T.internal_flux.κ))

# The balance solved is
# 
#   Θo - Θₛ
# κ ------- = Qn (all fluxes positive upwards)
#     δ
#
# Where the LHS is the internal diffusive flux inside the ocean and the RHS are the 
# atmospheric and radiative fluxes (provided explicitly and iterated upon).
@inline flux_balance_temperature(F::DiffusiveFlux, θo, Qn) = θo - Qn / F.κ * F.δ

# Change 𝒬₀ as a function of incoming and outgoing fluxes. The flaw here is that
# the ocean emissivity and albedo are fixed, but they might be a function of the 
# surface temperature, so we might need to pass actually the radiation and the 
# albedo and emissivity as arguments.
@inline function retrieve_temperature(st::DiagnosticSurfaceTemperature, θ₀, ℂ, 𝒬₀, 
                                      ρₐ, cₚ, ℰv, Σ★, ρₒ, cpₒ, g, 
                                      prescribed_heat_fluxes, 
                                      radiation_properties)

    Rd = prescribed_heat_fluxes # net downwelling radiation (positive out of the ocean)
    
    # upwelling radiation is calculated explicitly 
    # TODO: we could calculate it semi-implicitly as ϵσTⁿ⁺¹Tⁿ³
    Ru = upwelling_radiation(θ₀, radiation_properties) 
    Rn = Rd + Ru

    u★ = Σ★.momentum
    θ★ = Σ★.temperature
    q★ = Σ★.water_vapor
 
    # sensible heat flux + latent heat flux (positive out of the ocean)
    Qs = - ρₐ * u★ * (cₚ * θ★ + q★ * ℰv)

    # Net heat flux (positive out of the ocean)
    Qn = (Qs + Rn) / ρₒ / cpₒ

    θo = AtmosphericThermodynamics.air_temperature(ℂ, 𝒬₀)

    return flux_balance_temperature(st.internal_flux, θo, Qn)
end