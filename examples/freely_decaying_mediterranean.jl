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
λ₁, λ₂  = ( 0, 42) # domain in longitude
φ₁, φ₂  = (30, 45) # domain in latitude
# A stretched vertical grid with a Δz of 1.5 meters in the first 50 meters
z_faces = stretched_vertical_faces(depth = 5000, 
                             surface_layer_Δz = 2.5, 
                             stretching = PowerLawStretching(1.070), 
                             surface_layer_height = 50)

Nx = 15 * 42 # 1 / 15th of a degree resolution
Ny = 15 * 15 # 1 / 15th of a degree resolution
Nz = length(z_faces) - 1

grid = LatitudeLongitudeGrid(CPU();
                             size = (Nx, Ny, Nz),
                             latitude  = (φ₁, φ₂),
                             longitude = (λ₁, λ₂),
                             z = z_faces,
                             halo = (7, 7, 7))

# Interpolating the bathymetry onto the grid
#
# We regrid the bathymetry onto the grid.
# we allow a minimum depth of 10 meters (all shallower regions are 
# considered land) and we use 25 intermediate grids (interpolation_passes = 25)
# Note that more interpolation passes will smooth the bathymetry
bottom_height = regrid_bathymetry(grid, 
                                  height_above_water = 1,
                                  minimum_depth = 10,
                                  interpolation_passes = 25)

# Let's use an active cell map to elide computation in inactive cells
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map = true)

# Correct oceananigans
import Oceananigans.Advection: nothing_to_default

nothing_to_default(user_value; default) = isnothing(user_value) ? default : user_value

# Constructing the model
#
# We construct a model that evolves two tracers, temperature (:T), salinity (:S)
# We do not have convection since the simulation is just slumping to equilibrium so we do not need a turbulence closure
# We select a linear equation of state for buoyancy calculation, and WENO schemes for both tracer and momentum advection.
# The free surface utilizes a split-explicit formulation with a barotropic CFL of 0.75 based on wave speed.
model = HydrostaticFreeSurfaceModel(; grid,
                            momentum_advection = WENOVectorInvariant(),
                              tracer_advection = WENO(grid; order = 7),
                                  free_surface = SplitExplicitFreeSurface(; cfl = 0.75, grid),
                                      buoyancy = SeawaterBuoyancy(),
                                      tracers  = (:T, :S, :c),
                                      coriolis = HydrostaticSphericalCoriolis(scheme = ActiveCellEnstrophyConserving()))

# Initializing the model
#
# the model can be initialized with custom values or with ecco2 fields.
# In this case, our ECCO2 dataset has access to a temperature and a salinity
# field, so we initialize T and S from ECCO2. 
# We initialize our passive tracer with a surface blob near to the coasts of Libia
@info "initializing model"
libia_blob(x, y, z) = z > -20 || (x - 15)^2 + (y - 34)^2 < 1.5 ? 1 : 0

set!(model, T = ECCO2Metadata(:temperature), S = ECCO2Metadata(:salinity), c = libia_blob)

fig = Figure()
ax  = Axis(fig[1, 1])
heatmap!(ax, interior(model.tracers.T, :, :, Nz), colorrange = (10, 20), colormap = :thermal)
ax  = Axis(fig[1, 2])
heatmap!(ax, interior(model.tracers.S, :, :, Nz), colorrange = (35, 40), colormap = :haline)

simulation = Simulation(model, Δt = 10minutes, stop_time = 10*365days)

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

# Record a video
#
# Let's read the data and record a video of the Mediterranean Sea's surface
# (1) Zonal velocity (u)
# (2) Meridional velocity (v)
# (3) Temperature (T)
# (4) Salinity (S)
u_series = FieldTimeSeries("med_surface_field.jld2", "u")
v_series = FieldTimeSeries("med_surface_field.jld2", "v")
T_series = FieldTimeSeries("med_surface_field.jld2", "T")
S_series = FieldTimeSeries("med_surface_field.jld2", "S")
c_series = FieldTimeSeries("med_surface_field.jld2", "c")
iter = Observable(1)

u = @lift begin
    f = interior(u_series[$iter], :, :, 1)
    f[f .== 0] .= NaN
    f
end
v = @lift begin
    f = interior(v_series[$iter], :, :, 1)
    f[f .== 0] .= NaN
    f
end
T = @lift begin
    f = interior(T_series[$iter], :, :, 1)
    f[f .== 0] .= NaN
    f
end
S = @lift begin
    f = interior(S_series[$iter], :, :, 1)
    f[f .== 0] .= NaN
    f
end
c = @lift begin
    f = interior(c_series[$iter], :, :, 1)
    f[f .== 0] .= NaN
    f
end

fig = Figure()
ax  = Axis(fig[1, 1], title = "surface zonal velocity ms⁻¹")
heatmap!(u)
ax  = Axis(fig[1, 2], title = "surface meridional velocity ms⁻¹")
heatmap!(v)
ax  = Axis(fig[2, 1], title = "surface temperature ᵒC")
heatmap!(T)
ax  = Axis(fig[2, 2], title = "surface salinity psu")
heatmap!(S)
ax  = Axis(fig[2, 3], title = "passive tracer -")
heatmap!(c)

GLMakie.record(fig, "mediterranean_video.mp4", 1:length(u_series.times); framerate = 5) do i
    @info "recording iteration $i"
    iter[] = i    
end