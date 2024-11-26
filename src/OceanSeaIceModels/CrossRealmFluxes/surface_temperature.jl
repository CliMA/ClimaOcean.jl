import Thermodynamics as AtmosphericThermodynamics  

struct PrescribedSurfaceTemperature end

regularize_surface_temperature(::PrescribedSurfaceTemperature, grid) = PrescribedSurfaceTemperature()

struct DiagnosticSurfaceTemperature{I}
    internal_flux :: I
end

struct DiffusiveFlux{Δz, K}
    Δz :: Z
    κ  :: K # diffusivity in W m⁻² K⁻¹
end

function DiffusiveFlux(grid; κ = 0.2) 
    Δz = @allowscalar Δzᶜᶜᶜ(1, 1, grid.Nz, grid)
    return DiffusiveFlux(Δz, κ)
end

regularize_surface_temperature(T::DiagnosticSurfaceTemperature{<:DiffusiveFlux}, grid) =
    DiagnosticSurfaceTemperature(DiffusiveFlux(grid; κ = T.internal_flux.κ))

@inline flux_balance_temperature(F::DiffusiveFlux, θo, Qn) = θo + Qn / F.κ * F.Δz

# Do nothing (just copy the temperature)
@inline retrieve_temperature(::PrescribedSurfaceTemperature, θ₀, args...) = θ₀

# Change 𝒬₀ as a function of incoming and outgoing fluxes. The flaw here is that
# the ocean emissivity and albedo are fixed, but they might be a function of the 
# surface temperature, so we might need to pass actually the radiation and the 
# albedo and emissivity as arguments.
@inline function retrieve_temperature(st::DiagnosticSurfaceTemperature, θ₀, ℂ, 𝒬₀, 
                                      ρₐ, cₚ, ℰv, Σ★, 
                                      prescribed_heat_fluxes, 
                                      radiation_properties)

    Rd = prescribed_heat_fluxes # net downwelling radiation 
    ϵ  = radiation_properties.ocean_emissivity
    σ  = radiation_properties.stefan_boltzmann_constant

    Ru = ϵ * σ * θ₀^4 # upwelling radiation (calculated explicitly)
    Rn = Ru - Rd # net radiation

    u★ = Σ★.momentum
    θ★ = Σ★.temperature
    q★ = Σ★.water_vapor

    # sensible heat flux + latent heat flux
    Qs = - ρₐ * u★ * (cₚ * θ★ + q★ * ℰv)

    # Net heat flux (positive out of the ocean)
    Qn = Qs + Rn
    θo = AtmosphericThermodynamics.air_temperature(ℂ, 𝒬₀)

    # surface temperature calculated as a balance of fluxes
    return flux_balance_temperature(st.internal_flux, θo, Qn)
end
