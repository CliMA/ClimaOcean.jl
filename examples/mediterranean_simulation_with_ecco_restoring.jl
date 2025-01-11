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
using ClimaOcean.ECCO: ECCO4Monthly
using Oceananigans.Units
using Printf

using CFTime
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
z_faces = stretched_vertical_faces(depth = 5000, 
                                   surface_layer_Δz = 2.5, 
                                   stretching = PowerLawStretching(1.070), 
                                   surface_layer_height = 50)

Nx = 15 * Int(λ₂ - λ₁) # 1/15th of a degree resolution
Ny = 15 * Int(φ₂ - φ₁) # 1/15th of a degree resolution
Nz = length(z_faces) - 1

grid = LatitudeLongitudeGrid(GPU();
                             size = (Nx, Ny, Nz),
                             latitude  = (φ₁, φ₂),
                             longitude = (λ₁, λ₂),
                             z = z_faces,
                             halo = (7, 7, 7))

# ### Bathymetry Interpolation
#
# The script interpolates bathymetric data onto the grid, ensuring that the model accurately represents 
# the sea floor's topography. Parameters such as `minimum_depth` and `interpolation_passes`
# are adjusted to refine the bathymetry representation.

bottom_height = regrid_bathymetry(grid, 
                                  minimum_depth = 10,
                                  interpolation_passes = 10,
                                  major_basins = 1)

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

# ## ECCO Restoring
#
# The model is restored to at the surface to the temperature and salinity fields from the ECCO dataset.
# We build the restoring using the `ECCORestoring` functionality. 
# This allows us to nudge the model towards realistic temperature and salinity profiles.
# `ECCORestoring` accepts a `mask` keyword argument to restrict the restoring region.

@inline surface_mask(x, y, z, t) = z > - 50 

start = DateTimeProlepticGregorian(1993, 1, 1)
stop  = DateTimeProlepticGregorian(1993, 12, 1)
dates = range(start, stop; step=Month(1))

FT = ECCORestoring(:temperature, grid; dates, mask=surface_mask, rate=1/5days)
FS = ECCORestoring(:salinity, grid;    dates, mask=surface_mask, rate=1/5days)

# Constructing the Simulation
#
# We construct an ocean simulation that evolves two tracers, temperature (:T), salinity (:S)
# and we pass the previously defined forcing that nudge these tracers 

ocean = ocean_simulation(grid; forcing = (T=FT, S=FS))

# Initializing the model
#
# The model can be initialized with custom values or with ecco fields.
# In this case, our ECCO dataset has access to a temperature and a salinity
# field, so we initialize temperature T and salinity S from ECCO.

set!(ocean.model, T=ECCOMetadata(:temperature; dates=dates[1]), 
                  S=ECCOMetadata(:salinity;    dates=dates[1]))

fig = Figure()
ax  = Axis(fig[1, 1])
heatmap!(ax, view(ocean.model.tracers.T, :, :, Nz), colorrange = (10, 20), colormap = :thermal)
ax  = Axis(fig[1, 2])
heatmap!(ax, view(ocean.model.tracers.S, :, :, Nz), colorrange = (35, 40), colormap = :haline)

save("initial_conditions.png", fig)
# ![](initial_conditions.png)

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

# Let's reset the maximum number of iterations and we can not increase the time step size

ocean.Δt = 3minutes
ocean.stop_iteration = Inf
ocean.stop_time = 100days
model = ocean.model

ocean.output_writers[:surface_fields] = JLD2OutputWriter(model, merge(model.velocities, model.tracers);
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

# u_series = FieldTimeSeries("med_surface_field.jld2", "u"; backend=OnDisk())
# v_series = FieldTimeSeries("med_surface_field.jld2", "v"; backend=OnDisk())
# T_series = FieldTimeSeries("med_surface_field.jld2", "T"; backend=OnDisk())
# S_series = FieldTimeSeries("med_surface_field.jld2", "S"; backend=OnDisk())
# iter = Observable(1)

# u = @lift(u_series[$iter])
# v = @lift(v_series[$iter])
# T = @lift(T_series[$iter])
# S = @lift(S_series[$iter])

# fig = Figure()
# ax  = Axis(fig[1, 1], title = "surface zonal velocity ms⁻¹")
# heatmap!(ax, u)
# ax  = Axis(fig[1, 2], title = "surface meridional velocity ms⁻¹")
# heatmap!(ax, v)
# ax  = Axis(fig[2, 1], title = "surface temperature ᵒC")
# heatmap!(ax, T)
# ax  = Axis(fig[2, 2], title = "surface salinity psu")
# heatmap!(ax, S)

# CairoMakie.record(fig, "mediterranean_video.mp4", 1:length(u_series.times); framerate = 5) do i
#     @info "recording iteration $i"
#     iter[] = i    
# end
# # ![](mediterranean_video.mp4)
