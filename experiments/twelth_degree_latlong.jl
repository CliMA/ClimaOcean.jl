using Oceananigans
using Oceananigans: architecture
using ClimaOcean
using ClimaOcean.ECCO2
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: CATKEVerticalDiffusivity
using Oceananigans.Coriolis: ActiveCellEnstrophyConserving
using Oceananigans.Units
using ClimaOcean.OceanSimulations
using ClimaOcean.VerticalGrids: exponential_z_faces
using ClimaOcean.JRA55
using Printf

#####
##### Global Ocean at 1/12th of a degree
#####

# 100 vertical levels
z_faces = exponential_z_faces(100, 6000)

Nx = 360
Ny = 150
Nz = length(z_faces) - 1

arch = Distributed(CPU(), partition = Partition(8))

@show grid = load_balanced_regional_grid(arch; 
                                         size = (Nx, Ny, Nz), 
                                         z = z_faces, 
                                         latitude  = (-75, 75),
                                         longitude = (0, 360),
                                         halo = (7, 7, 7))
          
#####
##### The Ocean component
#####                             

simulation = ocean_simulation(grid)

model = simulation.model

# Initializing the model
set!(model, 
     T = ECCO2Metadata(:temperature), 
     S = ECCO2Metadata(:salinity),
     e = 1e-6)

#####
##### The atmosphere
#####

backend    = JRA55NetCDFBackend(10) # InMemory(8)
atmosphere = JRA55_prescribed_atmosphere(arch; backend)
radiation  = Radiation()

coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)

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

# Simulation warm up!
#
# We have regridded from a coarse solution (1/4er of a degree) to a
# fine grid (1/15th of a degree). Also, the bathymetry has little mismatches 
# that might crash our simulation. We warm up the simulation with a little 
# time step for few iterations to allow the solution to adjust to the new_grid
# bathymetry
simulation.Δt = 10
simulation.stop_iteration = 1000
run!(simulation)

# Run the real simulation
#
# Now that the solution has adjusted to the bathymetry we can ramp up the time
# step size. We use a `TimeStepWizard` to automatically adapt to a cfl of 0.2
wizard = TimeStepWizard(; cfl = 0.2, max_Δt = 10minutes, max_change = 1.1)

simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# Let's reset the maximum number of iterations
simulation.stop_iteration = Inf

simulation.output_writers[:surface_fields] = JLD2OutputWriter(model, merge(model.velocities, model.tracers);
                                                              indices = (:, :, Nz),
                                                              schedule = TimeInterval(1day),
                                                              overwrite_existing = true,
                                                              filename = "med_surface_field")

run!(simulation)

