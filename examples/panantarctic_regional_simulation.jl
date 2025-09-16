# # Panantarctic regional ocean simulation
#
# This example sets up and runs a regional ocean simulation for a domain around Antarctica
# using the Oceananigans.jl and ClimaOcean.jl. The simulation covers latitudes from 80°S to 20°S,
# with a horizontal resolution of 1/4 degree and 60 vertical levels.
#
# The simulation's results are visualized with the CairoMakie.jl package.
#
# The example showcases how we can use DatasetRestoring functionality to restore to a given dataset
# and also how to add custom masks on forcings.
#
# ## Initial setup with package imports
#
# We begin by importing the necessary Julia packages for visualization (CairoMakie),
# ocean modeling (Oceananigans, ClimaOcean), handling dates and times (Dates),
# and CUDA for running on CUDA-enabled GPUs.
# These packages provide the foundational tools for setting up the simulation environment,
# including grid setup, physical processes modeling, and data visualization.

using ClimaOcean
using ClimaOcean.ECCO
using Oceananigans
using Oceananigans.Units
using CUDA
using CairoMakie
using Dates
using Printf

# ### Grid and Bathymetry
#
# We start by constructing a the grid around Antarctica.

arch = GPU()
Nx = 1440
Ny = 240
Nz = 40

depth = 6000meters
z = ExponentialCoordinate(Nz, -depth, 0)

underlying_grid = LatitudeLongitudeGrid(arch;
                                        size = (Nx, Ny, Nz),
                                        halo = (7, 7, 7),
                                        z = z,
                                        latitude  = (-80, -20),
                                        longitude = (0, 360))

# ### Bathymetry and immersed boundary
#
# We add the bottom height from ETOPO1 data on our grid via `regrid_bathymetry`:

bottom_height = regrid_bathymetry(underlying_grid;
                                  minimum_depth = 10meters,
                                  interpolation_passes = 5,
                                  major_basins = 1)

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_height); active_cells_map=true)

# and visualise it:

fig, ax, hm = heatmap(grid.immersed_boundary.bottom_height,
                      colormap=:deep, colorrange=(-depth, 0))
Colorbar(fig[0, 1], hm, label="Bottom height (m)", vertical=false)
save("panantarctic_bathymetry.png", fig)
nothing #hide

# ![](panantarctic_bathymetry.png)

# ### Restoring force
#
# We need to add restoring forces (both for tracers and for velocities) at the
# northern part of the domain to mimic the effect of the part of the ocean north
# of 20°S that we are not including in our domain.
#
# First we create some mask that have the following form:
#
# ```julia
#               φS                   φN             -20
# -------------- | ------------------ | ------------ |
# no restoring   0    linear mask     1   mask = 1   1
# ```

const φN₁ = -23
const φN₂ = -25
const φS₁ = -78
const φS₂ = -75

@inline northern_mask(φ)    = min(max((φ - φN₂) / (φN₁ - φN₂), zero(φ)), one(φ))
@inline southern_mask(φ, z) = ifelse(z > -20,
                                     min(max((φ - φS₂) / (φS₁ - φS₂), zero(φ)), one(φ)),
                                     zero(φ))

@inline function tracer_mask(λ, φ, z, t)
     n = northern_mask(φ)
     s = southern_mask(φ, z)
     return max(s, n)
end

# Now we are ready to construct the forcing. We relax temperature, salinity to
# data from the ECCO4 dataset at a timescale of 5 days. For the velocities at the
# northern part of the domain, we apply a sponge layer (i.e., we relax them to zero).

start_date = DateTime(1993, 1, 1)
end_date   = DateTime(1993, 4, 1)

dataset = ECCO4Monthly()
T_meta = Metadata(:temperature; start_date, end_date, dataset)
S_meta = Metadata(:salinity;    start_date, end_date, dataset)

@inline function u_sponge(i, j, k, grid, clock, fields, p)
     φ = Oceananigans.Grids.φnode(i, j, k, grid, Face(), Center(), Center())
     return - p.rate * fields.u[i, j, k] * northern_mask(φ)
end

@inline function v_sponge(i, j, k, grid, clock, fields, p)
     φ = Oceananigans.Grids.φnode(i, j, k, grid, Center(), Face(), Center())
     return - p.rate * fields.v[i, j, k] * northern_mask(φ)
end

rate = 1/5days
forcing = (T = DatasetRestoring(T_meta, grid; rate, mask=tracer_mask),
           S = DatasetRestoring(S_meta, grid; rate, mask=tracer_mask),
           u = Forcing(u_sponge; discrete_form=true, parameters=(; rate)),
           v = Forcing(v_sponge; discrete_form=true, parameters=(; rate)))

# ### Ocean model configuration
#
# We build our ocean model using `ocean_simulation`,

momentum_advection = WENOVectorInvariant()
tracer_advection   = WENO(order=7)

ocean = ocean_simulation(grid; forcing, momentum_advection, tracer_advection)

# We initialize the ocean model with eddy-resolving ECCO2 temperature, salinity,
# and horizontal velocities.

dataset = ECCO2Monthly()
T_meta = Metadatum(:temperature; date=start_date, dataset)
S_meta = Metadatum(:salinity;    date=start_date, dataset)
u_meta = Metadatum(:u_velocity;  date=start_date, dataset)
v_meta = Metadatum(:v_velocity;  date=start_date, dataset)

set!(ocean.model, T=T_meta, S=S_meta, u=u_meta, v=v_meta)

# ### Prescribed atmosphere and radiation
#
# Next we build a prescribed atmosphere state and radiation model,
# which will drive the ocean simulation.

backend    = JRA55NetCDFBackend(41)
atmosphere = JRA55PrescribedAtmosphere(arch; backend)
radiation  = Radiation(arch)

# ## The coupled simulation

# We put all the pieces together (ocean, atmosphere, and radiation)
# into a coupled model and a coupled simulation.
# We start with a small-ish time step of 2 minutes.
# We run the simulation for 10 days with this small-ish time step.

coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)
simulation    = Simulation(coupled_model; Δt=2minutes, stop_time = 10days)

# A callback function to monitor the simulation's progress is always useful.

wall_time = [time_ns()]

function progress(sim)
    ocean = sim.model.ocean
    u, v, w = ocean.model.velocities
    T = ocean.model.tracers.T

    Tmax, Tmin = maximum(T), minimum(T)
    umax = maximum(abs, u), maximum(abs, v), maximum(abs, w)
    step_time = 1e-9 * (time_ns() - wall_time[1])

    @info @sprintf("Time: %s, Iteration %d, Δt %s, max(vel): (%.2e, %.2e, %.2e), max(T): %.2f, min(T): %.2f, wtime: %s \n",
                   prettytime(ocean.model.clock.time),
                   ocean.model.clock.iteration,
                   prettytime(ocean.Δt),
                   umax..., Tmax, Tmin, prettytime(step_time))

     wall_time[1] = time_ns()
end

simulation.callbacks[:progress] = Callback(progress, TimeInterval(5days))

# ### Output
#
# We use output writers to save the simulation data at regular intervals.

ocean.output_writers[:surface] = JLD2Writer(ocean.model, merge(ocean.model.tracers, ocean.model.velocities);
                                            schedule = TimeInterval(1days),
                                            filename = "panantarctic_surface_fields",
                                            indices = (:, :, grid.Nz),
                                            overwrite_existing = true,
                                            array_type = Array{Float32})

# ### Spinning up the simulation
#
# We spin up the simulation with a small time step to ensure that the interpolated initial
# conditions adapt to the model numerics and parameterization without causing instability.
# A 10-day integration with a time step of 1 minute should be sufficient to dissipate spurious
# initialization shocks.

run!(simulation)
nothing #hide

# ### Running the simulation
#
# Now that the simulation has spun up, we can run increase the timestep and run for longer;
# here we choose 60 days.

simulation.stop_time = 60days
simulation.Δt = 10minutes
run!(simulation)
nothing #hide

# ## Visualizing the results
#
# The simulation has finished, let's visualize the results.
# In this section we pull up the saved data and create visualizations using the CairoMakie.jl package.
# In particular, we generate an animation of the evolution of surface fields:
# surface speed (s), surface temperature (T), and turbulent kinetic energy (e).

u = FieldTimeSeries("panantarctic_surface_fields.jld2", "u"; backend = OnDisk())
v = FieldTimeSeries("panantarctic_surface_fields.jld2", "v"; backend = OnDisk())
T = FieldTimeSeries("panantarctic_surface_fields.jld2", "T"; backend = OnDisk())
e = FieldTimeSeries("panantarctic_surface_fields.jld2", "e"; backend = OnDisk())

times = u.times
Nt = length(times)

n = Observable(Nt)

land = interior(T.grid.immersed_boundary.bottom_height) .>= 0

Tn = @lift begin
    Tn = interior(T[$n])
    Tn[land] .= NaN
    view(Tn, :, :, 1)
end

en = @lift begin
    en = interior(e[$n])
    en[land] .= NaN
    view(en, :, :, 1)
end

un = Field{Face, Center, Nothing}(u.grid)
vn = Field{Center, Face, Nothing}(v.grid)

s = @at (Center, Center, Nothing) sqrt(un^2 + vn^2)
s = Field(s)

sn = @lift begin
    parent(un) .= parent(u[$n])
    parent(vn) .= parent(v[$n])
    compute!(s)
    sn = interior(s)
    sn[land] .= NaN
    view(sn, :, :, 1)
end

title = @lift string("Panantarctic regional ocean simulation after ",
                     prettytime(times[$n] - times[1]))

λ, φ, _ = nodes(T)

fig = Figure(size = (1000, 1500))

axs = Axis(fig[1, 1], xlabel="Longitude (deg)", ylabel="Latitude (deg)")
axT = Axis(fig[2, 1], xlabel="Longitude (deg)", ylabel="Latitude (deg)")
axe = Axis(fig[3, 1], xlabel="Longitude (deg)", ylabel="Latitude (deg)")

hm = heatmap!(axs, λ, φ, sn, colorrange = (0, 0.5), colormap = :deep, nan_color=:lightgray)
Colorbar(fig[1, 2], hm, label = "Surface Speed (m s⁻¹)")

hm = heatmap!(axT, λ, φ, Tn, colorrange = (-1, 30), colormap = :magma, nan_color=:lightgray)
Colorbar(fig[2, 2], hm, label = "Surface Temperature (ᵒC)")

hm = heatmap!(axe, λ, φ, en, colorrange = (0, 1e-3), colormap = :solar, nan_color=:lightgray)
Colorbar(fig[3, 2], hm, label = "Turbulent Kinetic Energy (m² s⁻²)")

Label(fig[0, :], title)

save("acc_snapshot.png", fig)
nothing #hide

# ![](snapshot.png)

# And now we make a movie:

CairoMakie.record(fig, "panantarctic_regional_surface.mp4", 1:Nt, framerate = 8) do nn
    n[] = nn
end
nothing #hide

# ![](panantarctic_regional_surface.mp4)
