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

#####
##### Global Ocean at 1/12th of a degree
#####

# 100 vertical levels
z_faces = exponential_z_faces(10, 6000)

Nx = 4320
Ny = 1800
Nz = length(z_faces) - 1

arch = Distributed(GPU(), partition = Partition(8))

@show grid = load_balanced_regional_grid(arch; 
                                         size = (Nx, Ny, Nz), 
                                         z = z_faces, 
                                         latitude  = (-75, 75),
                                         longitude = (0, 360),
                                         halo = (7, 7, 7),
                                         interpolation_passes = 25,
                                         bathymetry_file = "bathymetry1.jld2")
          
#####
##### The Ocean component
#####                             

ocean = ocean_simulation(grid)
model = ocean.model

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

ocean.callbacks[:progress] = Callback(progress, IterationInterval(10))

# Simulation warm up!
coupled_model.Δt = 10
coupled_model.stop_iteration = 1000
run!(coupled_model)

# Run the real simulation
#
# Now that the solution has adjusted to the bathymetry we can ramp up the time
# step size. We use a `TimeStepWizard` to automatically adapt to a cfl of 0.35
wizard = TimeStepWizard(; cfl = 0.35, max_Δt = 10minutes, max_change = 1.1)

coupled_model.ocean.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# Let's reset the maximum number of iterations
coupled_model.stop_iteration = Inf

run!(simulation)

