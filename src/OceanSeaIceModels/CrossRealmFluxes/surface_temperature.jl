using CUDA: @allowscalar
import Thermodynamics as AtmosphericThermodynamics  

####
#### Prescribed surface temperature
####

struct PrescribedSurfaceTemperature end

regularize_surface_temperature(::PrescribedSurfaceTemperature, grid) = PrescribedSurfaceTemperature()

# Do nothing (just copy the temperature)
@inline retrieve_temperature(::PrescribedSurfaceTemperature, θ₀, args...) = θ₀

####
#### Diagnostic surface temperature calculated as a flux balance
####

struct DiagnosticSurfaceTemperature{I}
    internal_flux :: I
end

struct DiffusiveFlux{Z, K}
    Δz :: Z
    κ  :: K # diffusivity in m²s⁻¹
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

@inline flux_balance_temperature(F::DiffusiveFlux, θo, Qn) = θo + Qn / F.κ * F.Δz

# Change 𝒬₀ as a function of incoming and outgoing fluxes. The flaw here is that
# the ocean emissivity and albedo are fixed, but they might be a function of the 
# surface temperature, so we might need to pass actually the radiation and the 
# albedo and emissivity as arguments.
@inline function retrieve_temperature(st::DiagnosticSurfaceTemperature, θ₀, ℂ, 𝒬₀, 
                                      ρₐ, cₚ, ℰv, Σ★, ρₒ, cpₒ,
                                      prescribed_heat_fluxes, 
                                      radiation_properties)

    Rd = prescribed_heat_fluxes # net downwelling radiation (positive out of the ocean)
    
    # upwelling radiation is calculated explicitly 
    # TODO: we could calculate it semi-implicitly as ϵσTⁿ⁺¹Tⁿ³
    Ru = upwelling_radiation(θ₀, radiation_properties) 
    Rn = Ru + Rd # net radiation

    u★ = Σ★.momentum
    θ★ = Σ★.temperature
    q★ = Σ★.water_vapor

    # sensible heat flux + latent heat flux
    Qs = - ρₐ * u★ * (cₚ * θ★ + q★ * ℰv)

    # Net heat flux (positive out of the ocean)
    Qn = (Qs + Rn) / ρₒ / cpₒ
    θo = AtmosphericThermodynamics.air_temperature(ℂ, 𝒬₀)

    # surface temperature calculated as a balance of fluxes
    return flux_balance_temperature(st.internal_flux, θo, Qn)
end

@inline upwelling_radiation(θ₀, ::Nothing) = zero(θ₀)
@inline upwelling_radiation(θ₀, r) = r.σ * r.ϵ * θ₀^4