using GLMakie
using Oceananigans
using Oceananigans: architecture
using ClimaOcean
using ClimaOcean.ECCO2
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: CATKEVerticalDiffusivity
using Oceananigans.Coriolis: ActiveCellEnstrophyConserving
using Oceananigans.Units
using Printf

#####
##### Regional Mediterranean grid 
#####

# Domain and Grid size
#
# We construct a grid that represents the Mediterranean sea, 
# with a resolution of 1/10th of a degree (roughly 10 km resolution)
λ₁, λ₂ = ( 0, 42) # domain in longitude
φ₁, φ₂ = (30, 45) # domain in latitude
# A stretched vertical grid with a Δz of 1.5 meters in the first 50 meters
z_face = stretched_vertical_faces(depth = 5000, 
                             surface_layer_Δz = 2.5, 
                             stretching = PowerLawStretching(1.070), 
                             surface_layer_height = 50)

Nx = 4 * 42 # 1 / 4th of a degree resolution
Ny = 4 * 15 # 1 / 4th of a degree resolution
Nz = length(z) - 1

grid = LatitudeLongitudeGrid(CPU();
                             size = (Nx, Ny, Nz),
                             latitude  = (φ₁, φ₂),
                             longitude = (λ₁, λ₂),
                             z = z_face,
                             halo = (7, 7, 7))

# Interpolating the bathymetry onto the grid
#
# We regrid the bathymetry onto the grid we constructed.
# we allow a minimum depth of 10 meters (all shallower regions are 
# considered land) and we use 5 intermediate grids (interpolation_passes = 5)
# Note that more interpolation passes will smooth the bathymetry
heigth = regrid_bathymetry(grid, 
                           height_above_water = 1,
                           minimum_depth = 10,
                           interpolation_passes = 5)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(height))

# Correct oceananigans
import Oceananigans.Advection: nothing_to_default

nothing_to_default(user_value; default) = isnothing(user_value) ? default : user_value

# Constructing the model
#
# We construct the model with three tracers, temperature (:T), salinity (:S) and TKE (:e).
# The latter is used as a prognostic variable for the micro-scale turbulence closure: `CATKEVerticalDiffusivity`
# We select a linear equation of state for buoyancy calculation, and WENO schemes for both tracer and momentum advection.
# The free surface utilizes a split-explicit formulation with a barotropic CFL of 0.75 based on wave speed.
model = HydrostaticFreeSurfaceModel(; grid,
                            momentum_advection = WENOVectorInvariant(),
                              tracer_advection = WENO(grid; order = 7),
                                  free_surface = SplitExplicitFreeSurface(; cfl = 0.75, grid),
                                      closure  = CATKEVerticalDiffusivity(),
                                      buoyancy = SeawaterBuoyancy(),
                                      tracers  = (:T, :S, :e),
                                      coriolis = HydrostaticSphericalCoriolis(scheme = ActiveCellEnstrophyConserving()))

# Initializing the model
#
# the model can be initialized with custom values or with ecco2 fields.
# In this case, our ECCO2 dataset has access to a temperature and a salinity
# field, so we initialize T and S from ECCO2. We initialize TKE (:e) with a 
# constant value of 1e-6 m²/s² throughout the domain
@info "initializing model from ECCO2 fields"
initialize!(model, T = :ecco2_temperature, S = :ecco2_salinity, e = 1e-6)

fig = Figure()
ax  = Axis(fig[1, 1])
heatmap!(ax, interior(model.tracers.T, :, :, Nz), colorrange = (10, 20), colormap = :thermal)
ax  = Axis(fig[1, 2])
heatmap!(ax, interior(model.tracers.S, :, :, Nz), colorrange = (35, 40), colormap = :haline)

simulation = Simulation(model, Δt = 20, stop_iteration = 100, stop_time = 10*365days)

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

simulation.stop_iteration = Inf

wizard = TimeStepWizard(; cfl = 0.2, max_Δt = 10minutes, max_change = 1.1)

simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

simulation.output_writers[:surface_fields] = JLD2OutputWriter(model, merge(model.velocities, model.tracers);
                                                              indices = (:, :, Nz),
                                                              schedule = TimeInterval(1day),
                                                              overwrite_existing = true,
                                                              filename = "med_surface_field")

run!(simulation)



