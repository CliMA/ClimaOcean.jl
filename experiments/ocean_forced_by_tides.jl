
using ClimaOcean
using Oceananigans
using Oceananigans.BoundaryConditions: FieldBoundaryConditions
using Oceananigans.Grids
using OrthogonalSphericalShellGrids
using Oceananigans.Units
using Printf

Φ = FieldTimeSeries("tidal_potential_jra55.jld2", "Φ")
Φ = FieldTimeSeries("tidal_potential_jra55.jld2", "Φ"; boundary_conditions=FieldBoundaryConditions(Φ.grid, (Center, Center, Nothing)))
Oceananigans.BoundaryConditions.fill_halo_regions!(Φ)

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

z = MutableVerticalDiscretization([-100, 0])

# A barotropic ocean
grid = TripolarGrid(; size = (700, 400, 1), halo = (7, 7, 7), z)
bottom_height = regrid_bathymetry(grid)
grid  = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)
ocean = ocean_simulation(grid; closure=nothing, tracer_advection=nothing, equation_of_state)

# A prescribed tidal forcing
atmos = PrescribedAtmosphere(Φ.grid, Φ.times; tidal_potential=Φ)
set!(ocean.model, T=first(atmos.tracers.T) + 273.15, S=35)

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
barotropic_earth = Simulation(barotropic_earth_model, Δt=20minutes, stop_time=360days) 

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

ocean.output_writers[:velocities] = JLD2OutputWriter(ocean.model, ocean.model.velocities;
                                                     filename="barotropic_earth_velocities.jld2",
                                                     overwrite_existing=true,
                                                     schedule=IterationInterval(40))

ocean.output_writers[:free_surface] = JLD2OutputWriter(ocean.model, (; η = ocean.model.free_surface.η);
                                                       filename="barotropic_earth_free_surface.jld2",
                                                       overwrite_existing=true,
                                                       schedule=IterationInterval(40))

add_callback!(barotropic_earth, progress, IterationInterval(100))

run!(barotropic_earth)

u = FieldTimeSeries("barotropic_earth_velocities.jld2", "u"; backend=InMemory(100))
v = FieldTimeSeries("barotropic_earth_velocities.jld2", "v"; backend=InMemory(100))
# η = FieldTimeSeries("barotropic_earth_free_surface.jld2", "η")
# w = FieldTimeSeries("barotropic_earth_velocities.jld2", "w")

using GLMakie, JLD2

file = jldopen("barotropic_earth_free_surface.jld2")
iters = keys(file["timeseries/t"])

fig = Figure()
axu = Axis(fig[1, 1], title = "u-velocity")
axv = Axis(fig[1, 2], title = "v-velocity")
axη = Axis(fig[2, 1], title = "free surface")
#axw = Axis(fig[1, 3], title = "w-velocity")

iter = Observable(1)

un = @lift(interior(u[$iter], :, :, 1))
vn = @lift(interior(v[$iter], :, :, 1))
ηn = @lift(file["timeseries/η/" * iters[$iter]][:, :, 1])
# wn = @lift(interior(w[$iter], :, :, 2))

heatmap!(axu, un, colorrange=(-0.1, 0.1))
heatmap!(axv, vn, colorrange=(-0.1, 0.1))
heatmap!(axη, ηn, colorrange=(-0.7, 0.7))

record(fig, "tides.mp4", 1:length(u.times), framerate = 8) do i
    iter[] = i
    @info "doing iter $i"
end