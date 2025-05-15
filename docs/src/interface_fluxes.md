# Interface fluxes

`ClimaOcean`'s `OceanSeaIceModel` has essentially two goals:
    1. Manage time-stepping multiple component models forward simulatneously,
    2. Compute and/or pass fluxes between the component models.

This page describes how OceanSeaIceModel computes fluxes at the interfaces between
component models.

## Computing the neutral drag coefficient

```@example interface_fluxes
using Oceananigans
using ClimaOcean

# Atmosphere velocities
Nx = Ny = 200
uₐ = range(0.5, stop=40, length=Nx) # winds at 10 m, m/s

# Ocean state parameters
T₀ = 20   # Surface temperature, ᵒC
S₀ = 35   # Surface salinity

x = y = (0, 1)
z = (-1, 0)
atmos_grid = RectilinearGrid(size=(Nx, Ny); x, y, topology=(Periodic, Periodic, Flat))
ocean_grid = RectilinearGrid(size=(Nx, Ny, 1); x, y, z, topology=(Periodic, Periodic, Bounded))

# Build the atmosphere
atmosphere = PrescribedAtmosphere(atmos_grid, surface_layer_height=10)
interior(atmosphere.tracers.T) .= 273.15 + T₀ # K
interior(atmosphere.velocities.u, :, :, 1, 1) .= uₐ # m/s

ocean = ocean_simulation(ocean_grid, momentum_advection=nothing,
                         tracer_advection=nothing, closure=nothing)
set!(ocean.model, T=T₀, S=S₀)

atmosphere_ocean_flux_formulation = SimilarityTheoryFluxes(stability_functions=nothing)
interfaces = ComponentInterfaces(atmosphere, ocean; atmosphere_ocean_flux_formulation)
model = OceanSeaIceModel(ocean; atmosphere)
```

The drag coefficient is defined as

```math
C^D \equiv \frac{u_\star^2}{| \Delta \bm{u} |^2}
\qquad \text{where} \qquad
\Delta \bm{u} \equiv \left ( u_a - u_o \right ) \bm{\hat x} + \left ( v_a - v_o \right ) \bm{\hat y} \, .
```

is the difference between the atmosphere velocity at the
prescribed `surface_layer_height` and the ocean surface velocity,
and``u_\star`` is the friction velocity. 

```@example interface_fluxes
using CairoMakie
u★ = model.interfaces.atmosphere_ocean_interface.fluxes.friction_velocity
u★ = interior(u★, :, 1, 1)

c₁ = 0.0027
c₂ = 0.000142
c₃ = 0.0000764
u★_LY = @. sqrt(c₁ * uₐ + c₂ * uₐ^2 + c₃ * uₐ^3)

fig = Figure(size=(800, 400))
axu = Axis(fig[1, 1], xlabel="uₐ (m s⁻¹) at 10 m", ylabel="u★ (m s⁻¹)")
lines!(axu, uₐ, u★, label="SimilarityTheoryFluxes")
lines!(axu, uₐ, u★_LY, label="Polynomial fit \n from Large and Yeager (2009)")
axislegend(axu, position=:lt)

axd = Axis(fig[1, 2], xlabel="uₐ (m s⁻¹) at 10 m", ylabel="1000 × Cᴰ")
Cᴰ = @. (u★ / uₐ)^2
Cᴰ_LY = @. (u★_LY / uₐ)^2
lines!(axd, uₐ, 1000 .* Cᴰ, label="SimilarityTheoryFluxes")
lines!(axd, uₐ, 1000 .* Cᴰ_LY, label="Polynomial fit \n from Large and Yeager (2009)")
axislegend(axd, position=:rt)

fig
```
