using Printf
using Oceananigans
using Oceananigans.Units
using Oceananigans: architecture
using ClimaOcean
using ClimaOcean.ECCO2
using ClimaOcean.OceanSimulations
using ClimaOcean.OceanSeaIceModels

#####
##### Near - Global Ocean at 1/4th of a degree
#####

# 60 vertical levels
z_faces = exponential_z_faces(Nz=60, depth=6500)

Nx = 1440
Ny = 600
Nz = length(z_faces) - 1

arch = GPU() 

# A near-global grid from 75ᵒ S to 75ᵒ N
grid = LatitudeLongitudeGrid(arch; 
                             size = (Nx, Ny, Nz), 
                             halo = (7, 7, 7), 
                             z = z_faces, 
                             longitude = (0, 360),
                             latitude = (-75, 75))

# We retrieve the bathymetry from the ETOPO1 data by ensuring a 
# minimum depth of 10 meters (everyhting shallower is considered land)
# and removing all connected regions (inland lakes)
bottom_height = retrieve_bathymetry(grid, bathymetry_file; 
                                    minimum_depth = 10,
                                    dir = "./",
                                    interpolation_passes = 20,
                                    connected_regions_allowed = 0)
 
# An immersed boundary using a staircase representation of bathymetry
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map = true) 

#####
##### The Ocean component
#####                             

# We retain all the defaults for the ocean model
ocean = ocean_simulation(grid) 
model = ocean.model

# We interpolate the initial conditions from the ECCO2 dataset 
# (for the moment these are both 1st January 1992)
set!(model, 
     T = ECCO2Metadata(:temperature),
     S = ECCO2Metadata(:salinity))

#####
##### The atmosphere
#####

# The whole prescribed atmosphere is loaded in memory 
# 4 snapshots at the time.
backend    = JRA55NetCDFBackend(4) 
atmosphere = JRA55_prescribed_atmosphere(arch; backend)

# Tabulated ocean albedo from Payne (1982)
# ocean emissivity is the default 0.97
radiation  = Radiation(arch)

#####
##### The atmospheric-forced coupled ocean-seaice model
#####

# Simplistic sea ice that ensure "no-cooling-fluxes" where `T < T_minimum`
# to change with a thermodynamic sea-ice model
sea_ice = ClimaOcean.OceanSeaIceModels.MinimumTemperatureSeaIce()

# The complete coupled model
coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

wall_time = [time_ns()]

# We define a progress function that
# shows the maximum values of velocity and temperature
# to make sure everything proceedes as planned
function progress(sim) 
    u, v, w = sim.model.velocities  
    T = sim.model.tracers.T

    Tmax = maximum(interior(T))
    Tmin = minimum(interior(T))
    umax = maximum(interior(u)), maximum(interior(v)), maximum(interior(w))
    step_time = 1e-9 * (time_ns() - wall_time[1])

    @info @sprintf("Time: %s, Iteration %d, Δt %s, max(vel): (%.2e, %.2e, %.2e), max(trac): %.2f, %.2f, wtime: %s \n",
                   prettytime(sim.model.clock.time),
                   sim.model.clock.iteration,
                   prettytime(sim.Δt),
                   umax..., Tmax, Tmin, prettytime(step_time))

     wall_time[1] = time_ns()
end

ocean.callbacks[:progress] = Callback(progress, IterationInterval(10))

fluxes = (u = model.velocities.u.boundary_conditions.top.condition,
          v = model.velocities.v.boundary_conditions.top.condition,
          T = model.tracers.T.boundary_conditions.top.condition,
          S = model.tracers.S.boundary_conditions.top.condition)


output_kwargs = (; overwrite_existing = true, array_type = Array{Float32})

# We add a couple of outputs: the surface state every 12 hours, the surface fluxes
# every 12 hours, and the whole state every ten days
ocean.output_writers[:fluxes] = JLD2OutputWriter(model, fluxes;
                                                  schedule = TimeInterval(0.5days),
                                                  overwrite_existing = true,
                                                  filename = "surface_fluxes",
                                                  output_kwargs...)

ocean.output_writers[:surface] = JLD2OutputWriter(model, merge(model.tracers, model.velocities);
                                                  schedule = TimeInterval(0.5days),
                                                  filename = "surface",
                                                  indices = (:, :, grid.Nz),
                                                  output_kwargs...)

ocean.output_writers[:snapshots] = JLD2OutputWriter(model, merge(model.tracers, model.velocities);
                                                    schedule = TimeInterval(10days)
                                                    filename = "snapshots",
                                                    output_kwargs...)

# Checkpointer for restarting purposes
ocean.output_writers[:checkpoint] = Checkpointer(model;
                                                 schedule = TimeInterval(60days),
                                                 overwrite_existing = true,
                                                 prefix = "checkpoint")

#####
##### Run the simulation!
#####

# warm up the simulation to ensure that the model adapts 
# to the interpolated initial conditions without crashing
ocean.Δt = 10
ocean.stop_iteration = 1
wizard = TimeStepWizard(; cfl = 0.1, max_Δt = 90, max_change = 1.1)
ocean.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

# Finally, the coupled simulation!
coupled_simulation = Simulation(coupled_model; Δt=1, stop_time = 30days)

run!(coupled_simulation)

# Now that the initial conditions and the initialization shocks have been smoothen out by 
# the initial transient, it is time to start the actual simulation!
wizard = TimeStepWizard(; cfl = 0.35, max_Δt = 540, max_change = 1.1)
ocean.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# Let's reset the maximum number of iterations and specify a 20years stop time
coupled_model.ocean.stop_time = 7200days
coupled_simulation.stop_time = 7200days
coupled_model.ocean.stop_iteration = Inf
coupled_simulation.stop_iteration = Inf

# Let's run it!
run!(coupled_simulation)
