# # Nonhydrostatic nighttime LES
#
# This script assembles a two-dimensional large-eddy simulation that resolves the
# nighttime marine boundary layer with a prescribed, dry atmosphere and no sea ice.
# The goal is to keep the setup minimal so we can diagnose coupling issues before
# adding more physics or richer initial conditions.

using ClimaOcean
using Oceananigans
using Oceananigans.Units
using Printf

# ## Ocean LES configuration
#
# We keep a 2D (x–z) grid with 4 m horizontal and 2 m vertical spacing, spanning
# 256 m × 128 m.  The Coriolis parameter corresponds to 33° N.

Nx, Ny, Nz = 64, 1, 64
x = y = (0, 256)
z = (-128, 0)
arch = CPU()
grid = RectilinearGrid(arch; size = (Nx, Ny, Nz), x, y, z, topology = (Periodic, Periodic, Bounded), halo=(5, 1, 5))
coriolis = FPlane(latitude=33)

ocean = nonhydrostatic_ocean_simulation(grid; coriolis)
Tᵢ(x, y, z) = 25 + 0.01 * rand()
Sᵢ(x, y, z) = 35 + 0.01 * rand()
set!(ocean.model, T=Tᵢ, S=Sᵢ) # warm ocean, uniform salinity

# ## Prescribed nighttime atmosphere
#
# We prescribe a clear-sky, dry atmosphere: zero humidity, no downwelling radiation,
# and a weak 5 m s⁻¹ eastward wind that cools the ocean via sensible heat fluxes.

atmosphere_grid = RectilinearGrid(arch; size=(Nx, Ny), x, y, topology = (Periodic, Periodic, Flat))
atmosphere_times = range(0, 8hours, length=3)
atmosphere = PrescribedAtmosphere(atmosphere_grid, collect(atmosphere_times))

parent(atmosphere.tracers.T) .= 290 # Kelvin
parent(atmosphere.tracers.q) .= 0   # dry
parent(atmosphere.velocities.u) .= 5 # m s⁻¹
parent(atmosphere.velocities.v) .= 0
parent(atmosphere.downwelling_radiation.shortwave) .= 0
parent(atmosphere.downwelling_radiation.longwave) .= 0

# ## Coupled ocean–atmosphere model

coupled_model = OceanSeaIceModel(ocean; atmosphere)
#simulation = Simulation(coupled_model; Δt=1second, stop_time=8hours)
simulation = Simulation(coupled_model; Δt=1second, stop_iteration=100)

heat_flux = coupled_model.interfaces.net_fluxes.ocean_surface.Q
τx = coupled_model.interfaces.net_fluxes.ocean_surface.u

function progress(sim)
    ocean = sim.model.ocean
    u, v, w = ocean.model.velocities
    T = ocean.model.tracers.T
    Tmax = maximum(interior(T))
    Tmin = minimum(interior(T))
    umax = maximum(abs, u)
    vmax = maximum(abs, v)
    wmax = maximum(abs, w)
    Qmax = maximum(sim.model.interfaces.net_fluxes.ocean_surface.Q)
    Qmin = minimum(sim.model.interfaces.net_fluxes.ocean_surface.Q)
    τxmax = maximum(sim.model.interfaces.net_fluxes.ocean_surface.u)
    τxmin = minimum(sim.model.interfaces.net_fluxes.ocean_surface.u)

    @info @sprintf("Time: %s, iteration: %d, Δt: %s, max|u|: %.2e, max|v|: %.2e, max|w|: %.2e, max(Q): %.2e, min(Q): %.2e, max(τx): %.2e, min(τx): %.2e",
                   prettytime(sim), iteration(sim), prettytime(sim.Δt), umax, vmax, wmax, Qmax, Qmin, τxmax, τxmin)
end

add_callback!(simulation, progress, IterationInterval(10))

run!(simulation)