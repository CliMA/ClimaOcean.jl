using Printf
using Oceananigans
using Oceananigans.Units
using ClimaOcean
using OrthogonalSphericalShellGrids
using Oceananigans
using Oceananigans: architecture
using ClimaOcean
using ClimaOcean.ECCO2
using Oceananigans.Units
using ClimaOcean.OceanSimulations
using ClimaOcean.OceanSeaIceModels
using ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: Radiation
using ClimaOcean.VerticalGrids: exponential_z_faces
using ClimaOcean.JRA55
using ClimaOcean.JRA55: JRA55NetCDFBackend, JRA55_prescribed_atmosphere
using ClimaOcean.Bathymetry

#####
##### Near - Global Ocean at 1/4th of a degree
#####

bathymetry_file = nothing # "bathymetry_tmp.jld2"

# 60 vertical levels
z_faces = exponential_z_faces(Nz=60, depth=6500)

Nx = 1440
Ny = 600
Nz = length(z_faces) - 1

arch = GPU() 

grid = LatitudeLongitudeGrid(arch; 
                             size = (Nx, Ny, Nz), 
                             halo = (7, 7, 7), 
                             z = z_faces, 
                             longitude = (0, 360),
                             latitude = (-75, 75))

bottom_height = retrieve_bathymetry(grid, bathymetry_file; 
                                    minimum_depth = 10,
                                    dir = "./",
                                    interpolation_passes = 20,
                                    connected_regions_allowed = 0)
 
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map = true) 

#####
##### The Ocean component
#####                             

# We retain all the defaults for the ocean model
ocean = ocean_simulation(grid) 
model = ocean.model

set!(model, 
     T = ECCO2Metadata(:temperature),
     S = ECCO2Metadata(:salinity))

#####
##### The atmosphere
#####

backend    = JRA55NetCDFBackend(4) 
atmosphere = JRA55_prescribed_atmosphere(arch; backend)

# Tabulated albedo from Payne (1982)
radiation  = Radiation(arch)

# Simplistic sea ice that ensure "no-cooling-fluxes" where `T < T_minimum`
# to change with a thermodynamic sea-ice model
sea_ice = ClimaOcean.OceanSeaIceModels.MinimumTemperatureSeaIce()

coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

wall_time = [time_ns()]

function progress(sim) 
    u, v, w = sim.model.velocities  
    T, S = sim.model.tracers

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

ocean.output_writers[:fluxes] = JLD2OutputWriter(model, fluxes,
                                                  schedule = TimeInterval(0.5days),
                                                  overwrite_existing = true,
                                                  array_type = Array{Float32},
                                                  filename = "surface_fluxes")

ocean.output_writers[:surface] = JLD2OutputWriter(model, merge(model.tracers, model.velocities),
                                                  schedule = TimeInterval(0.5days),
                                                  overwrite_existing = true,
                                                  array_type = Array{Float32},
                                                  filename = "surface",
                                                  indices = (:, :, grid.Nz))

ocean.output_writers[:snapshots] = JLD2OutputWriter(model, merge(model.tracers, model.velocities),
                                                    schedule = TimeInterval(10days),
                                                    overwrite_existing = true,
                                                    array_type = Array{Float32},
                                                    filename = "snapshots")

ocean.output_writers[:checkpoint] = Checkpointer(model, 
                                                 schedule = TimeInterval(60days),
                                                 overwrite_existing = true,
                                                 prefix = "checkpoint")

# Simulation warm up!
ocean.Δt = 10
ocean.stop_iteration = 1
wizard = TimeStepWizard(; cfl = 0.1, max_Δt = 90, max_change = 1.1)
ocean.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

stop_time = 30days

coupled_simulation = Simulation(coupled_model; Δt=1, stop_time)

run!(coupled_simulation)

wizard = TimeStepWizard(; cfl = 0.4, max_Δt = 540, max_change = 1.1)
ocean.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# Let's reset the maximum number of iterations
coupled_model.ocean.stop_time = 7200days
coupled_simulation.stop_time = 7200days
coupled_model.ocean.stop_iteration = Inf
coupled_simulation.stop_iteration = Inf

run!(coupled_simulation)
