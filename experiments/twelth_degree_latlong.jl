using MPI
MPI.Init()

using Oceananigans
using Oceananigans: architecture
using ClimaOcean
using ClimaOcean.ECCO2
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: CATKEVerticalDiffusivity
using Oceananigans.Coriolis: ActiveCellEnstrophyConserving
using Oceananigans.Units
using ClimaOcean.OceanSimulations
using ClimaOcean.OceanSeaIceModels
using ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: Radiation
using ClimaOcean.VerticalGrids: exponential_z_faces
using ClimaOcean.JRA55
using ClimaOcean.JRA55: JRA55NetCDFBackend, JRA55_prescribed_atmosphere
using Printf
using JLD2

#####
##### Global Ocean at 1/12th of a degree
#####

# import Oceananigans.Fields: set!
# using Oceananigans.Grids: on_architecture
# using Oceananigans.DistributedComputations: DistributedField, global_size, partition_global_array
# using Oceananigans.Architectures: architecture

# # Automatically partition under the hood if sizes are compatible
# function set!(u::DistributedField, v::Array)
#      gsize = global_size(architecture(u), size(u))
 
#      @show gsize, size(v)
#      if size(v) == gsize
#          f = partition_global_array(architecture(u), v, size(u))
#          u .= f
#          return u
#      else
#          try
#              f = on_architecture(architecture(u), v)
#              u .= f
#              return u
     
#          catch
#              throw(ArgumentError("ERROR: DimensionMismatch: array could not be set to match destination field"))
#          end
#      end
#  end
 
bathymetry_file = "bathymetry12.jld2"

# 100 vertical levels
z_faces = exponential_z_faces(Nz=100, depth=6000)

Nx = 4320
Ny = 1800
Nz = length(z_faces) - 1

bottom = zeros(Nx, Ny, 1)

arch = Distributed(GPU(), partition = Partition(8))

grid = LatitudeLongitudeGrid(arch; 
			     size = (Nx, Ny, Nz),
                             z = z_faces,
                             latitude  = (-75, 75),
                             longitude = (0, 360),
                             halo = (7, 7, 7))

bottom .= jldopen(bathymetry_file)["bathymetry"]

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom); active_cells_map = true)

# grid = load_balanced_regional_grid(arch; 
#                                    size = (Nx, Ny, Nz), 
#                                    z = z_faces, 
#                                    latitude  = (-75, 75),
#                                    longitude = (0, 360),
#                                    halo = (7, 7, 7),
#                                    interpolation_passes = 20,
#                                    minimum_depth = 10,
#                                    connected_regions_allowed = 3, # We allow the oceans, the med, the bering sea
#                                    bathymetry_file)
 
@show grid                                   
                              
#####
##### The Ocean component
#####                             

free_surface = SplitExplicitFreeSurface(; grid, cfl=0.7, fixed_Δt=270)

ocean = ocean_simulation(grid; free_surface)
model = ocean.model

# Initializing the model
set!(model, 
     T = ECCO2Metadata(:temperature), 
     S = ECCO2Metadata(:salinity),
     e = 1e-6)

#####
##### The atmosphere
#####

backend    = JRA55NetCDFBackend(5) 
atmosphere = JRA55_prescribed_atmosphere(arch; backend)
radiation  = Radiation()

coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)

wall_time = [time_ns()]

function progress(sim) 
    u, v, w = sim.model.velocities  
    T, S = sim.model.tracers

    step_time = 1e-9 * (time_ns() - wall_time[1])

    @info @sprintf("Time: %s, Iteration %d, Δt %s, max(vel): (%.2e, %.2e, %.2e), max(trac): %.2f, %.2f, wtime: %s \n",
                   prettytime(sim.model.clock.time),
                   sim.model.clock.iteration,
                   prettytime(sim.Δt),
                   maximum(abs, u), maximum(abs, v), maximum(abs, w),
                   maximum(abs, T), maximum(abs, S), step_time)

     wall_time[1] = time_ns()
end

ocean.callbacks[:progress] = Callback(progress, IterationInterval(1))

ocean.output_writers[:surface] = JLD2OutputWriter(model, merge(model.tracers, model.velocities),
                                                  schedule = TimeInterval(30days),
                                                  overwrite_existing = true,
                                                  array_type = Array{Float32},
                                                  filename = "snapshots_$(arch.local_rank)")

ocean.output_writers[:checkpoint] = Checkpointer(model, 
                                                 schedule = TimeInterval(60days),
                                                 overwrite_existing = true,
                                                 prefix = "checkpoint_$(arch.local_rank)")

# Simulation warm up!
ocean.Δt = 10
ocean.stop_iteration = 1
wizard = TimeStepWizard(; cfl = 0.3, max_Δt = 90, max_change = 1.1)
ocean.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

stop_time = 5days

coupled_simulation = Simulation(coupled_model, Δt=1, stop_time=stop_time)

run!(coupled_simulation)

wizard = TimeStepWizard(; cfl = 0.35, max_Δt = 270, max_change = 1.1)
ocean.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# Let's reset the maximum number of iterations
coupled_model.ocean.stop_time = 7200days

run!(coupled_simulation)

