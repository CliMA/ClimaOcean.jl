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
using CUDA
using CairoMakie
using SixelTerm
CUDA.device!(1)

include("correct_oceananigans.jl")

#####
##### Global Ocean at 1/12th of a degree
#####

bathymetry_file = nothing # "bathymetry_tmp.jld2"

# 100 vertical levels
z_faces = exponential_z_faces(Nz=60, depth=6500)

Nx = 720
Ny = 300
Nz = length(z_faces) - 1

arch = GPU() #Distributed(GPU(), partition = Partition(2))

grid = load_balanced_regional_grid(arch; 
                                   size = (Nx, Ny, Nz), 
                                   z = z_faces, 
                                   latitude  = (-75, 75),
                                   longitude = (0, 360),
                                   halo = (7, 7, 7),
                                   interpolation_passes = 20,
                                   minimum_depth = 10,
                                   connected_regions_allowed = 3, # We allow the oceans, the med, the bering sea
                                   bathymetry_file)
 
@show grid                                   
                              
#####
##### The Ocean component
#####                             

free_surface = SplitExplicitFreeSurface(; grid, cfl=0.7, fixed_Δt = 20minutes)
vertical_diffusivity = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), κ = 5e-5, ν = 5e-3)

closure = (RiBasedVerticalDiffusivity(), vertical_diffusivity)

ocean = ocean_simulation(grid; free_surface, closure) 
model = ocean.model

# Initializing the model
set!(model, "checkpoint_iteration497897.jld2")

#####
##### The atmosphere
#####

backend    = JRA55NetCDFBackend(10) 
atmosphere = JRA55_prescribed_atmosphere(arch; backend)
radiation  = Radiation()

sea_ice = ClimaOcean.OceanSeaIceModels.MinimumTemperatureSeaIce()

coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

wall_time = [time_ns()]

function progress(sim) 
    u, v, w = sim.model.velocities  
    T, S = sim.model.tracers

    Tmax = maximum(abs, T)
    Tmin = minimum(T)
    umax = maximum(abs, u), maximum(abs, v), maximum(abs, w)
    step_time = 1e-9 * (time_ns() - wall_time[1])

    @info @sprintf("Time: %s, Iteration %d, Δt %s, max(vel): (%.2e, %.2e, %.2e), max(trac): %.2f, %.2f, wtime: %s \n",
                   prettytime(sim.model.clock.time),
                   sim.model.clock.iteration,
                   prettytime(sim.Δt),
                   umax..., Tmax, Tmin, prettytime(step_time))

     wall_time[1] = time_ns()
end

ocean.callbacks[:progress] = Callback(progress, IterationInterval(100))

ocean.output_writers[:surface] = JLD2OutputWriter(model, merge(model.tracers, model.velocities),
                                                  schedule = TimeInterval(30days),
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

stop_time = 10days

coupled_simulation = Simulation(coupled_model, Δt=1, stop_time=stop_time)

try 
   run!(coupled_simulation)
catch 

end

wizard = TimeStepWizard(; cfl = 0.35, max_Δt = 15minutes, max_change = 1.1)
ocean.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# Let's reset the maximum number of iterations
coupled_model.ocean.stop_time = 7200days
coupled_simulation.stop_time = 7200days
coupled_model.ocean.stop_iteration = Inf
coupled_simulation.stop_iteration = Inf

try 
    run!(coupled_simulation)
catch 

end

