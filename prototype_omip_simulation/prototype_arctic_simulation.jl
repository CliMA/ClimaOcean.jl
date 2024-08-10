using Printf
using Oceananigans
using Oceananigans.Units
using ClimaOcean
using ClimaSeaIce
using OrthogonalSphericalShellGrids
using Oceananigans
using Oceananigans: architecture
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
using ClimaSeaIce.SeaIceDynamics: ExplicitMomentumSolver
using CairoMakie
using CFTime
using Dates

include("restoring_mask.jl")
include("simulation_outputs.jl")

#####
##### Global Ocean at 1/6th of a degree
#####

# 60 vertical levels
z_faces = exponential_z_faces(Nz=60, depth=6000)

Nx = 1440
Ny = 600
Nz = length(z_faces) - 1

arch = GPU() 

grid = TripolarGrid(arch; 
                    size = (Nx, Ny, Nz), 
                    halo = (7, 7, 7), 
                    z = z_faces, 
                    north_poles_latitude = 55,
                    first_pole_longitude = 75)

sea_ice_grid = TripolarGrid(arch; 
                            size = (Nx, Ny, 1), 
                            halo = (7, 7, 7), 
                            z = (-1, 0), 
                            north_poles_latitude = 55,
                            first_pole_longitude = 75)

bottom_height = retrieve_bathymetry(grid; 
                                    minimum_depth = 2,
                                    dir = "./",
                                    interpolation_passes = 20,
                                    connected_regions_allowed = 0)
 
grid         = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map = true) 
sea_ice_grid = ImmersedBoundaryGrid(sea_ice_grid, GridFittedBottom(bottom_height); active_cells_map = true) 

#####
##### The Ocean component
#####                             

tracer_advection = WENO(; order = 7)
free_surface = SplitExplicitFreeSurface(grid; cfl=0.7, fixed_Δt = 600)

ocean = ocean_simulation(grid; Δt = 10, tracer_advection, free_surface) 
model = ocean.model

#####
##### Sea ice simulation
#####

u_surface_velocity = interior(model.velocities.u, :, :, grid.Nz)
v_surface_velocity = interior(model.velocities.v, :, :, grid.Nz)

ocean_velocities = (; u = u_surface_velocity,
                      v = v_surface_velocity)

sea_ice = sea_ice_simulation(sea_ice_grid; ice_dynamics = nothing, advection = nothing) #ocean_velocities, ocean_ice_drag_coefficient = 5.5) 

#####
##### Setting the simulation with ECCO data
#####

date = DateTimeProlepticGregorian(1993, 1, 1)

temperature       = ECCOMetadata(:temperature, date, ECCO4Monthly())
salinity          = ECCOMetadata(:salinity,    date, ECCO4Monthly())
ice_thickness     = ECCOMetadata(:sea_ice_thickness,     date, ECCO4Monthly())
ice_concentration = ECCOMetadata(:sea_ice_area_fraction, date, ECCO4Monthly())

set!(ocean.model, T = temperature, S = salinity, e = 1e-6)

set!(sea_ice.model.ice_thickness,     ice_thickness;     maxiter = 2)
set!(sea_ice.model.ice_concentration, ice_concentration; maxiter = 2)

#####
##### A prescribed atmosphere
#####

backend    = JRA55NetCDFBackend(4) 
atmosphere = JRA55_prescribed_atmosphere(arch; backend)
radiation  = Radiation(arch; ocean_albedo = LatitudeDependentAlbedo())

#####
##### Building the complete simulation
#####

coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
coupled_simulation = Simulation(coupled_model; Δt = 10, stop_time = 30days)

wall_time = [time_ns()]

function progress(sim) 
    u, v, w = sim.model.ocean.model.velocities  
    T, S = sim.model.ocean.model.tracers
    h    = sim.model.sea_ice.model.ice_thickness

    Tmax = maximum(interior(T))
    hmax = maximum(interior(h))
    umax = maximum(interior(u)), maximum(interior(v)), maximum(interior(w))
    step_time = 1e-9 * (time_ns() - wall_time[1])

    @info @sprintf("Time: %s, Iteration %d, Δt %s, max(vel): (%.2e, %.2e, %.2e), max(trac): %.2f, %.2f, wtime: %s \n",
                   prettytime(sim.model.clock.time),
                   sim.model.clock.iteration,
                   prettytime(sim.Δt),
                   umax..., Tmax, hmax, prettytime(step_time))

     wall_time[1] = time_ns()
end

function update_dt!(simulation)
    simulation.Δt = simulation.model.ocean.Δt
    return nothing
end

coupled_simulation.callbacks[:progress]  = Callback(progress,   IterationInterval(10))
coupled_simulation.callbacks[:update_dt] = Callback(update_dt!, IterationInterval(1)) 

#### Saving outputs
set_outputs!(coupled_simulation)

# Run for 30 days!
conjure_time_step_wizard!(ocean; max_Δt = 50, max_change = 1.1, cfl = 0.1)
run!(coupled_simulation)

# # Finished the 30 days, run for other 720 days
coupled_simulation.stop_time = 1080days
ocean.stop_time = 1080days
sea_ice.stop_time = 1080days

conjure_time_step_wizard!(ocean; max_Δt = 600, max_change = 1.1, cfl = 0.25)
run!(coupled_simulation)
