using MPI
MPI.Init()

using Printf
using Oceananigans
using Oceananigans.Units
using ClimaOcean
using OrthogonalSphericalShellGrids
using Oceananigans
using Oceananigans: architecture
using Oceananigans.DistributedComputations: child_architecture
using Oceananigans.Grids: on_architecture
using Oceananigans.Coriolis: ActiveCellEnstrophyConserving
using Oceananigans.Units
using ClimaOcean
using ClimaOcean.OceanSimulations
using ClimaOcean.OceanSeaIceModels
using ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: Radiation, SimilarityTheoryTurbulentFluxes
using ClimaOcean.VerticalGrids: exponential_z_faces
using ClimaOcean.JRA55
using ClimaOcean.ECCO
using ClimaOcean.JRA55: JRA55NetCDFBackend, JRA55_prescribed_atmosphere
using ClimaOcean.ECCO: ECCO_restoring_forcing, ECCO4Monthly, ECCO2Daily, ECCOMetadata
using ClimaOcean.Bathymetry
using ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: LatitudeDependentAlbedo

using CFTime
using Dates

import Oceananigans.OutputReaders: cpu_interpolating_time_indices

cpu_interpolating_time_indices(arch::Distributed, args...) = cpu_interpolating_time_indices(child_architecture(arch), args...)

#####
##### Global Ocean at 1/6th of a degree
#####

bathymetry_file = nothing # "bathymetry_tmp.jld2"

# 60 vertical levels
z_faces = exponential_z_faces(Nz=60, depth=6000)

Nx = 4320
Ny = 2160
Nz = length(z_faces) - 1

arch = Distributed(GPU(), partition = Partition(x = 2, y = Equal()))
rank = arch.local_rank

grid = TripolarGrid(arch; 
                    size = (Nx, Ny, Nz), 
                    halo = (7, 7, 7), 
                    z = z_faces, 
                    north_poles_latitude = 55,
                    first_pole_longitude = 75)

bottom_height = retrieve_bathymetry(grid, bathymetry_file; 
                                    minimum_depth = 10,
				    interpolation_passes = 10)
 
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map = true) 

#####
##### The atmosphere component
#####

@info "Uploading an atmosphere on rank $(rank)..."
backend    = JRA55NetCDFBackend(4) 
atmosphere = JRA55_prescribed_atmosphere(arch; backend)
radiation  = Radiation(arch)

#####
##### The Ocean component
#####                             

free_surface = SplitExplicitFreeSurface(grid; substeps = 55)

#####
##### Add restoring to ECCO fields for temperature and salinity in the artic and antarctic
#####

const z_surface = Array(grid.zᵃᵃᶠ.parent)[end-8]
@info z_surface
@show z_surface

@inline function mask_T(λ, φ, z)
   if z < z_surface
      return zero(λ)
   else
      if φ < -70 || φ > 70
          return one(λ)
      elseif φ < -60
          return 1 / (-10) * (φ + 60)
      elseif φ > 60
          return 1 / 10 * (φ - 60)
      else
          return zero(λ)
      end
   end
end

@inline function mask_S(λ, φ, z)
   if z < z_surface
      return zero(λ)
   else
      if φ < -70 || φ > 70
         return one(λ)
      else
         return one(λ) * 0.1
      end
   end
end

maskT = CenterField(grid)
maskS = CenterField(grid)

set!(maskT, mask_T)
set!(maskS, mask_S)

dates = DateTimeProlepticGregorian(1993, 1, 1) : Month(1) : DateTimeProlepticGregorian(1993, 12, 1)

temperature = ECCOMetadata(:temperature, dates, ECCO4Monthly())
salinity    = ECCOMetadata(:salinity,    dates, ECCO4Monthly())

FT = ECCO_restoring_forcing(temperature; mask = maskT, grid, architecture = arch, timescale = 5days)
FS = ECCO_restoring_forcing(salinity;    mask = maskS, grid, architecture = arch, timescale = 5days)

forcing = (; T = FT, S = FS)

closure = RiBasedVerticalDiffusivity(; horizontal_Ri_filter = Oceananigans.TurbulenceClosures.FivePointHorizontalFilter())

@info "Building an ocean on rank $(rank)..."
ocean = ocean_simulation(grid; free_surface, forcing, closure) 
model = ocean.model

initial_date = dates[1]

#####
##### Coupling the different models...
#####

sea_ice = ClimaOcean.OceanSeaIceModels.MinimumTemperatureSeaIce()

@info "Coupling ocean and atmosphere on rank $(rank)..."
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

@info "Attaching Outputwriters on rank $(rank)..."

ocean.output_writers[:fluxes] = JLD2OutputWriter(model, fluxes,
                                                  schedule = TimeInterval(1days),
                                                  overwrite_existing = false,
                                                  array_type = Array{Float32},
                                                  filename = "surface_fluxes_$(arch.local_rank)")

ocean.output_writers[:surface] = JLD2OutputWriter(model, merge(model.tracers, model.velocities),
                                                  schedule = TimeInterval(1days),
                                                  overwrite_existing = false,
                                                  array_type = Array{Float32},
                                                  filename = "surface_$(arch.local_rank)",
                                                  indices = (:, :, grid.Nz))

ocean.output_writers[:snapshots] = JLD2OutputWriter(model, merge(model.tracers, model.velocities),
                                                    schedule = TimeInterval(10days),
                                                    overwrite_existing = false,
                                                    array_type = Array{Float32},
                                                    filename = "snapshots_$(arch.local_rank)")

ocean.output_writers[:checkpoint] = Checkpointer(model, 
                                                 schedule = TimeInterval(60days),
                                                 overwrite_existing = true,
                                                 prefix = "checkpoint_$(arch.local_rank)")


restart = nothing 

coupled_simulation = Simulation(coupled_model; Δt = 1, stop_time = 25days)
ocean.Δt = 10

if isnothing(restart)

    @info "Restarting from scratch on rank $(rank)..."
    # Set simulation from ECCO2 fields
    set!(ocean.model, 
         T = ECCOMetadata(:temperature, initial_date, ECCO2Daily()),
         S = ECCOMetadata(:salinity,    initial_date, ECCO2Daily()))
    
    # Simulation warm up!
    wizard = TimeStepWizard(; cfl = 0.1, max_Δt = 1.5minutes, max_change = 1.1)
    ocean.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

    run!(coupled_simulation)
else
    @info "Setting model from checkpoint $(restart) on rank $(rank)..."
    # Set the ocean from the restart file
    set!(ocean.model, restart)
end

# Remove the wizard
pop!(ocean.callbacks, :wizard)

ocean.Δt = 270
coupled_simulation.Δt = 270

# Let's reset the maximum number of iterations
coupled_model.ocean.stop_time = 6530days
coupled_simulation.stop_time = 6530days
coupled_model.ocean.stop_iteration = Inf
coupled_simulation.stop_iteration = Inf

run!(coupled_simulation)
