using GLMakie
using Oceananigans
using Oceananigans: architecture
using Oceananigans.Fields: interpolate!
using ClimaOcean
using ClimaOcean.ECCO2
using ClimaOcean.InitialConditions: three_dimensional_regrid!, adjust_tracers!
using Oceananigans.ImmersedBoundaries: mask_immersed_field!
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: CATKEVerticalDiffusivity
using Oceananigans.Coriolis: ActiveCellEnstrophyConserving
using Oceananigans.Units
using Printf

#####
##### Regional Mediterranean grid 
#####

# A stretched vertical grid with a Δz of 1.5 meters in the first 50 meters
z = stretched_vertical_faces(depth = 5000, 
                             surface_layer_Δz = 2.5, 
                             stretching = PowerLawStretching(1.070), 
                             surface_layer_height = 50)

Nx = 4 * 42 # 1 / 4th of a degree
Ny = 4 * 15 # 1 / 4th of a degree
Nz = length(z) - 1

@info "grid size: ($Nx, $Ny, $Nz)"

grid = LatitudeLongitudeGrid(CPU();
                             size = (Nx, Ny, Nz),
                             latitude = (30, 45),
                             longitude = (0, 42),
                             z,
                             halo = (7, 7, 7))

h = regrid_bathymetry(grid, height_above_water=1)

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(h))

# Correct oceananigans
import Oceananigans.Advection: nothing_to_default

nothing_to_default(user_value; default) = isnothing(user_value) ? default : user_value

# Construct the model and run it, will it run or we have to diffuse?
model = HydrostaticFreeSurfaceModel(; grid,
                                      momentum_advection = WENOVectorInvariant(),
                                      tracer_advection = WENO(grid; order = 7),
                                      free_surface = SplitExplicitFreeSurface(; cfl = 0.75, grid),
                                      closure  = CATKEVerticalDiffusivity(),
                                      buoyancy = SeawaterBuoyancy(),
                                      tracers  = (:T, :S, :e),
                                      coriolis = HydrostaticSphericalCoriolis(scheme = ActiveCellEnstrophyConserving()))

initialize!(model, T = :ecco2_temperature, S = :ecco2_salinity, e = 1e-6)

fig = Figure()
ax  = Axis(fig[1, 1])
heatmap!(ax, interior(model.tracers.T, :, :, Nz), colorrange = (10, 20), colormap = :thermal)
ax  = Axis(fig[1, 2])
heatmap!(ax, interior(model.tracers.S, :, :, Nz), colorrange = (35, 40), colormap = :haline)

simulation = Simulation(model, Δt = 20, stop_iteration = 100)

function progress(sim) 
    u, v, w = sim.model.velocities  
    T, S = sim.model.tracers

    @info @sprintf("Time: %s, Iteration %d, Δt %s, max(vel): (%.2e, %.2e, %.2e), max(T, S): %.2f, %.2f\n",
                   prettytime(sim.model.clock.time),
                   sim.model.clock.iteration,
                   prettytime(sim.Δt),
                   maximum(abs, u), maximum(abs, v), maximum(abs, w),
                   maximum(abs, T), maximum(abs, S))
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

# warm up for 100 iterations!
run!(simulation)

simulation.stop_time = 10*365days

wizard = TimeStepWizard(; cfl = 0.2, max_Δt = 10minutes, max_change = 1.1)

simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

simulation.output_writers[:surface_fields] = JLD2OutputWriter(model, merge(model.velocities, model.tracers);
                                                              indices = (:, :, Nz),
                                                              schedule = TimeInterval(1day),
                                                              overwrite_existing = true,
                                                              filename = "med_surface_field")

run!(simulation)



