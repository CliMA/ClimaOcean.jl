using CUDA: @allowscalar
import Thermodynamics as AtmosphericThermodynamics  

####
#### Utilities
####

@inline upwelling_radiation(θ₀, ::Nothing) = zero(θ₀)
@inline upwelling_radiation(θ₀, r) = r.σ * r.ϵ * θ₀^4

# For any surface temperture type that does not depend on the grid
regularize_surface_temperature_type(surface_temperature_type, grid) = surface_temperature_type

####
#### Prescribed surface temperature (the easiest case)
####

struct PrescribedSurfaceTemperature end

# Do nothing (just copy the temperature)
@inline compute_surface_temperature(::PrescribedSurfaceTemperature, θ₀, args...) = θ₀

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
function DiagnosticSurfaceTemperature(; κ = 0.1, δ = nothing) 
    internal_flux = DiffusiveFlux(; κ, δ)
    return DiagnosticSurfaceTemperature(internal_flux)
end

DiffusiveFlux(; κ = 1e-2, δ = nothing) = DiffusiveFlux(δ, κ)

function DiffusiveFlux(grid; κ = 0.1, δ = nothing)
    if isnothing(δ)
        δ = @allowscalar Δzᶜᶜᶜ(1, 1, grid.Nz, grid)
    end
    return DiffusiveFlux(δ, κ)
end

regularize_surface_temperature_type(T::DiagnosticSurfaceTemperature{<:DiffusiveFlux}, grid) =
    DiagnosticSurfaceTemperature(DiffusiveFlux(grid; κ = T.internal_flux.κ, δ = T.internal_flux.δ))

# The flux balance could be solved either
# 
#   Θₒ - Θₛⁿ⁺¹
# κ ---------- = Jᵀ (all fluxes positive upwards)
#       δ
#
# Where the LHS is the internal diffusive flux inside the ocean (within the boundary layer of thickness δ) 
# and the RHS are the atmospheric and radiative fluxes are provided explicitly, or
# 
#   Θₒ - Θₛⁿ⁺¹    σϵ θₛⁿ⁺¹θₛⁿ³
# κ ---------- - ------------ = Jᵀ (all fluxes positive upwards)
#       δ           ρₒ cpₒ
#
# Where the LHS is the internal diffusive flux inside the ocean (within the boundary layer of thickness δ) 
# plus the (semi-implicit) outgoing longwave flux and the RHS are the remaining atmospheric and radiative fluxes
# provided explicitly.
@inline flux_balance_temperature(F::DiffusiveFlux, θₒ, Jᵀ) = (θₒ - Jᵀ / F.κ * F.δ)

# he flaw here is that the ocean emissivity and albedo are fixed, but they might be a function of the 
# surface temperature, so we might need to pass the radiation and the albedo and emissivity as arguments.
@inline function compute_surface_temperature(st::DiagnosticSurfaceTemperature, θₛ, ℂ, 𝒬₀, 
                                            ρₐ, cₚ, ℰv, Σ★, ρₒ, cpₒ, g, 
                                            prescribed_heat_fluxes, 
                                            radiation_properties)

    Rd = prescribed_heat_fluxes # net downwelling radiation (positive out of the ocean)
    
    # upwelling radiation is calculated explicitly 
    Ru = upwelling_radiation(θₛ, radiation_properties) 
    Rn = Rd + Ru # Net radiation (positive out of the ocean)
    
    u★ = Σ★.momentum
    θ★ = Σ★.temperature
    q★ = Σ★.water_vapor
 
    # Turbulent heat fluxes, sensible + latent (positive out of the ocean)
    Qt = - ρₐ * u★ * (cₚ * θ★ + q★ * ℰv)

    # Net temperature flux (positive out of the ocean)
    Jᵀ = (Qt + Rn) / ρₒ / cpₒ 

    θₒ = AtmosphericThermodynamics.air_temperature(ℂ, 𝒬₀)

    return flux_balance_temperature(st.internal_flux, θₒ, Jᵀ)
end