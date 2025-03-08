
using ClimaOcean
using Oceananigans
using OrthogonalSphericalShellGrids
using Oceananigans.Units
using Printf

import Oceananigans.BuoyancyFormulations: ρ′ 
import ClimaOcean.OceanSeaIceModels: reference_density, heat_capacity

struct ZeroBuoyancy{R}
    reference_density :: R
end

@inline ρ′(i, j, k, grid, ::ZeroBuoyancy, T, S) = zero(grid)
@inline reference_density(r::ZeroBuoyancy) = r.reference_density
@inline heat_capacity(r::ZeroBuoyancy) = 3990.0

# all passive tracers
equation_of_state = ZeroBuoyancy(1026.0)

# A barotropic ocean
grid = TripolarGrid(size = (700, 400, 1), halo = (7, 7, 7), z = (-10, 0))
bottom_height = regrid_bathymetry(grid)
grid  = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)
ocean = ocean_simulation(grid; closure=nothing, bottom_drag_coefficient=0, tracer_advection=nothing, equation_of_state)

# A prescribed tidal forcing
Φ = FieldTimeSeries("tidal_potential_jra55.jld2", "Φ")
atmos = PrescribedAtmosphere(Φ.grid, Φ.times; tidal_potential=Φ)
set!(ocean.model, T=first(atmos.tracers.T), S=35)

# neutral radiation
radiation=Radiation(ocean_albedo=1, ocean_emissivity=0)

# No fluxes
struct ZeroFluxes{N, B}
    solver_stop_criteria :: N
    bulk_velocity :: B
end

import ClimaOcean.OceanSeaIceModels.InterfaceComputations: compute_interface_state, ComponentInterfaces
using ClimaOcean.OceanSeaIceModels.InterfaceComputations: zero_interface_state

@inline compute_interface_state(::ZeroFluxes, args...) = zero_interface_state(Float64)

interfaces = ComponentInterfaces(atmos, ocean, nothing;
                                 atmosphere_ocean_flux_formulation=ZeroFluxes(nothing, ClimaOcean.OceanSeaIceModels.InterfaceComputations.WindVelocity()),
                                 radiation)

barotropic_earth_model = OceanSeaIceModel(ocean, nothing; atmosphere=atmos, radiation, interfaces)

barotropic_earth = Simulation(barotropic_earth_model, Δt=20minutes, stop_time=2days)

wall_time = Ref(0.0)

function progress(barotropic_earth)
    clock   = barotropic_earth.model.clock
    u, v, w = barotropic_earth.model.ocean.model.velocities

    maxu = maximum(abs, u)
    maxv = maximum(abs, v)
    maxw = maximum(abs, w)
    
    @info @sprintf("Iteration: %d, time: %s, wall_time: %s, max(|u|, |v|, |w|): %.2e %.2e %.2e\n",
                   clock.iteration, prettytime(clock.time), prettytime(1e-9 * (time_ns() - wall_time[])), maxu, maxv, maxw)

    wall_time[] = time_ns()
end

add_callback!(barotropic_earth, progress, IterationInterval(100))

run!(barotropic_earth)