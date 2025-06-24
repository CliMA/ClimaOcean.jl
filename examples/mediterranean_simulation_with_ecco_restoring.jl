# # Mediterranean simulation with restoring to ECCO
#
# This example is a comprehensive example of setting up and running a high-resolution ocean
# simulation for the Mediterranean Sea using the Oceananigans and ClimaOcean packages, with
# a focus on restoring temperature and salinity fields from the ECCO (Estimating the Circulation
# and Climate of the Ocean) dataset.
#
# The example is divided into several sections, each handling a specific part of the simulation
# setup and execution process.

# ## Initial Setup with Package Imports
#
# We begin by importing necessary Julia packages for visualization (CairoMakie), ocean modeling
# (Oceananigans, ClimaOcean), and handling of dates and times (CFTime, Dates).
# These packages provide the foundational tools for creating the simulation environment,
# including grid setup, physical processes modeling, and data visualization.

using CairoMakie
using Oceananigans
using Oceananigans: architecture
using ClimaOcean
using ClimaOcean.ECCO
using Oceananigans.Units
using Printf
using Dates

# ## Grid Configuration for the Mediterranean Sea
#
# The script defines a high-resolution grid to represent the Mediterranean Sea, specifying the domain in terms of longitude (λ₁, λ₂),
# latitude (φ₁, φ₂), and a stretched vertical grid to capture the depth variation (`z_faces`).
# The grid resolution is set to approximately 1/15th of a degree, which translates to a spatial resolution of about 7 km.
# This section demonstrates the use of the LatitudeLongitudeGrid function to create a grid that matches the
# Mediterranean's geographical and bathymetric features.

λ₁, λ₂  = ( 0, 42) # domain in longitude
φ₁, φ₂  = (30, 45) # domain in latitude

z = StretchedInterfaces(; depth = 5000,
                   surface_layer_Δz = 2.5,
                   stretching = PowerLawStretching(1.070),
                   surface_layer_height = 50)

Nx = 15 * Int(λ₂ - λ₁) # 1/1G5 th of a degree resolution
Ny = 15 * Int(φ₂ - φ₁) # 1/15 th of a degree resolution
Nz = length(z)

grid = LatitudeLongitudeGrid(PU();
                             size = (Nx, Ny, Nz),
                             latitude  = (φ₁, φ₂),
                             longitude = (λ₁, λ₂),
                             z,
                             halo = (7, 7, 7))

# ### Bathymetry Interpolation
#
# The script interpolates bathymetric data onto the grid, ensuring that the model accurately represents
# the sea floor's topography. Parameters such as `minimum_depth` and `interpolation_passes`
# are adjusted to refine the bathymetry representation.

bottom_height = regrid_bathymetry(grid,
                                  height_above_water = 1,
                                  minimum_depth = 10,
                                  interpolation_passes = 25,
                                  connected_regions_allowed = 1)

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

# ## Downloading ECCO data
#
# The model is initialized with temperature and salinity fields from the ECCO dataset,
# using the function `ECCO_restoring_forcing` to apply restoring forcings for these tracers.
# This allows us to nudge the model towards realistic temperature and salinity profiles.

start_date = DateTime(1993, 1, 1)
end_date   = DateTime(1993, 12, 1)

FT = ECCO_restoring_forcing(:temperature; start_date, end_date, architecture = GPU(), timescale = 2days)
FS = ECCO_restoring_forcing(:salinity;    start_date, end_date, architecture = GPU(), timescale = 2days)

# Constructing the Simulation
#
# We construct an ocean simulation that evolves two tracers, temperature (:T), salinity (:S)
# and we pass the previously defined forcing that nudge these tracers

ocean = ocean_simulation(grid; forcing = (T = FT, S = FS))

# Initializing the model
#
# The model can be initialized with custom values or with ecco fields.
# In this case, our ECCO dataset has access to a temperature and a salinity
# field, so we initialize temperature T and salinity S from ECCO.

set!(ocean.model, T = ECCOMetadatum(:temperature; date=start_date),
                  S = ECCOMetadatum(:salinity;    date=start_date))

fig = Figure()
ax  = Axis(fig[1, 1])
heatmap!(ax, interior(model.tracers.T, :, :, Nz), colorrange = (10, 20), colormap = :thermal)
ax  = Axis(fig[1, 2])
heatmap!(ax, interior(model.tracers.S, :, :, Nz), colorrange = (35, 40), colormap = :haline)

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

# ## Simulation warm up!
#
# We have regridded from the coarse solution of the ECCO dataset (half of a degree) to a
# fine grid (1/15th of a degree). The bathymetry might also have little mismatches
# that might crash the simulation. We warm up the simulation with a little
# time step for few iterations to allow the solution to adjust to the new grid
# bathymetry.

ocean.Δt = 10
ocean.stop_iteration = 1000
run!(ocean)

# ## Run the real simulation
#
# Now that the solution has adjusted to the bathymetry we can ramp up the time
# step size. We use a `TimeStepWizard` to automatically adapt to a CFL of 0.2.

wizard = TimeStepWizard(; cfl = 0.2, max_Δt = 10minutes, max_change = 1.1)

ocean.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# Let's reset the maximum number of iterations
ocean.stop_iteration = Inf
ocean.stop_time = 200days

ocean.output_writers[:surface_fields] = JLD2Writer(model, merge(model.velocities, model.tracers);
                                                   indices = (:, :, Nz),
                                                   schedule = TimeInterval(1days),
                                                   overwrite_existing = true,
                                                   filename = "med_surface_field")

run!(ocean)

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

CairoMakie.record(fig, "mediterranean_video.mp4", 1:length(u_series.times); framerate = 5) do i
    @info "recording iteration $i"
    iter[] = i
end
