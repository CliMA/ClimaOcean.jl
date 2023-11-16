using GLMakie
using Oceananigans
using Oceananigans.Utils: prettytime
using ClimaOcean.Bathymetry: regrid_bathymetry
using ClimaOcean.ECCO2: ecco2_field, ecco2_center_mask
using ClimaOcean.VerticalGrids: stretched_vertical_faces, PowerLawStretching
using ClimaOcean.InitialConditions: three_dimensional_regrid!
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: CATKEVerticalDiffusivity
using Oceananigans.Coriolis: ActiveCellEnstrophyConserving

#####
##### Regional Mediterranean grid 
#####

function regrid_ecco_tracers_to_grid(grid)
    Tecco = ecco2_field(:temperature)
    Secco = ecco2_field(:salinity)
    
    ecco_tracers = (; Tecco, Secco)
    
    # Make sure all values are extended properly before regridding
    adjust_tracers!(ecco_tracers; mask = ecco2_center_mask())
    
    T = CenterField(grid)
    S = CenterField(grid)
    
    three_dimensional_regrid!(T, Tecco)
    three_dimensional_regrid!(S, Secco)
    
    return T, S
end

# A stretched vertical grid with a Δz of 1.5 meters in the first 50 meters
z = stretched_vertical_faces(minimum_depth = 5000, 
                             surface_layer_Δz = 1.75, 
                             stretching = PowerLawStretching(1.070), 
                             surface_layer_height = 50)

Nx = 20 * 55 # 1 / 20th of a degree
Ny = 20 * 25
Nz = length(z) - 1

grid = LatitudeLongitudeGrid(GPU();
                             size = (Nx, Ny, Nz),
                             latitude = (25, 50),
                             longitude = (0, 45),
                             z,
                             halo = (7, 7, 7))

h = regrid_bathymetry(grid, height_above_water=1)

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(h))

T, S = regrid_ecco_tracers_to_grid(grid)

mask_immersed_field!(T)
mask_immersed_field!(S)

fig = Figure()
ax  = Axis(fig[1, 1])
heatmap!(ax, interior(T, :, :, Nz), colorrange = (10, 20), colormap = :thermal)
ax  = Axis(fig[1, 2])
heatmap!(ax, interior(S, :, :, Nz), colorrange = (35, 40), colormap = :haline)

# Correct oceananigans
import Oceananigans.Advection: nothing_to_default

nothing_to_default(user_value; default) = isnothing(user_value) ? default : user_value

# Construct the model and run it, will it run or we have to diffuse?
model = HydrostaticFreeSurfaceModel(; grid,
                                      momentum_advection = WENOVectorInvariant(),
                                      tracer_advection = WENO(grid; order = 7),
                                      free_surface = SplitExplicitFreeSurface(; cfl = 0.75, grid),
                                      closure = CATKEVerticalDiffusivity(),
                                      buoyancy = SeawaterBuoyancy(),
                                      tracers = (:T, :S, :e),
                                      coriolis = HydrostaticSphericalCoriolis(scheme = ActiveCellEnstrophyConserving()))

set!(model, T = T, S = S)

# Probably we'll need to diffuse? Probably not let's see now

simulation = Simulation(model, Δt = 20, stop_time = 2days)

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

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

# warm up!
run!(simulation)

simulation.stop_time = 10*365days

wizard = TimeStepWizard(; cfl = 0.2, max_Δt = 2minutes, max_change = 1.1)

simulation.callbacks = Callback(wizard, IterationInterval(10))

run!(simulation)



