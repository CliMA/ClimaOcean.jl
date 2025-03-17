# # Near-global ocean simulation
#
# This example sets up and runs a near-global ocean simulation using the Oceananigans.jl and
# ClimaOcean.jl. The simulation covers latitudes from 75°S to 75°N with a horizontal
# resolution of 1/4 degree and 40 vertical levels.
#
# The simulation's results are visualized with the CairoMakie.jl package.
#
# ## Initial setup with package imports
#
# We begin by importing the necessary Julia packages for visualization (CairoMakie),
# ocean modeling (Oceananigans, ClimaOcean), and handling dates and times (CFTime, Dates).
# These packages provide the foundational tools for setting up the simulation environment,
# including grid setup, physical processes modeling, and data visualization.

using ClimaOcean
using ClimaOcean.ECCO
using Oceananigans
using Oceananigans.Units
using CairoMakie
using CFTime
using Dates
using Printf

# ### Grid configuration 
#
# We define a global grid with a horizontal resolution of 1/4 degree and 40 vertical levels.
# The grid is a `LatitudeLongitudeGrid` spanning latitudes from 75°S to 75°N.
# We use an exponential vertical spacing to better resolve the upper-ocean layers.
# The total depth of the domain is set to 6000 meters.
# Finally, we specify the architecture for the simulation, which in this case is a GPU.

arch = GPU()

Nx = 1440
Ny = 600
Nz = 40

depth = 6000meters
z_faces = exponential_z_faces(; Nz, depth)

grid = LatitudeLongitudeGrid(arch;
                             size = (Nx, Ny, Nz),
                             halo = (7, 7, 7),
                             z = z_faces,
                             latitude  = (-75, 75),
                             longitude = (0, 360))

# ### Bathymetry and immersed boundary
#
# We use `regrid_bathymetry` to derive the bottom height from ETOPO1 data.
# To smooth the interpolated data we use 5 interpolation passes. We also fill in
# (i) all the minor enclosed basins except the 3 largest `major_basins`, as well as
# (ii) regions that are shallower than `minimum_depth`.

bottom_height = regrid_bathymetry(grid;
                                  minimum_depth = 10meters,
                                  interpolation_passes = 5,
                                  major_basins = 3)

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)

# Let's see what the bathymetry looks like:

h = grid.immersed_boundary.bottom_height

fig, ax, hm = heatmap(h, colormap=:deep, colorrange=(-depth, 0))
Colorbar(fig[0, 1], hm, label="Bottom height (m)", vertical=false)
save("bathymetry.png", fig)
nothing #hide

# ![](bathymetry.png)

# ### Ocean model configuration
#
# We build our ocean model using `ocean_simulation`,

ocean = ocean_simulation(grid)

# which uses the default `ocean.model`,

ocean.model

# We initialize the ocean model with ECCO4 temperature and salinity for January 1, 1992.

set!(ocean.model, T=ECCOMetadatum(:temperature),
                  S=ECCOMetadatum(:salinity))

# ### Prescribed atmosphere and radiation
#
# Next we build a prescribed atmosphere state and radiation model,
# which will drive the ocean simulation. We use the default `Radiation` model,

# The radiation model specifies an ocean albedo emissivity to compute the net radiative
# fluxes. The default ocean albedo is based on Payne (1982) and depends on cloud cover
# (calculated from the ratio of maximum possible incident solar radiation to actual
# incident solar radiation) and latitude. The ocean emissivity is set to 0.97.

radiation = Radiation(arch)

# The atmospheric data is prescribed using the JRA55 dataset.
# The JRA55 dataset provides atmospheric data such as temperature, humidity, and winds
# to calculate turbulent fluxes using bulk formulae, see [`CrossRealmFluxes`](@ref).
# The number of snapshots that are loaded into memory is determined by
# the `backend`. Here, we load 41 snapshots at a time into memory.

atmosphere = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(41))

# ## The coupled simulation

# Next we assemble the ocean, atmosphere, and radiation
# into a coupled model,

coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)

# We then create a coupled simulation. We start with a small-ish time step of 90 seconds.
# We run the simulation for 10 days with this small-ish time step.

simulation = Simulation(coupled_model; Δt=90, stop_time=10days)

# We define a callback function to monitor the simulation's progress,

wall_time = Ref(time_ns())

function progress(sim)
    ocean = sim.model.ocean
    u, v, w = ocean.model.velocities
    T = ocean.model.tracers.T

    Tmax = maximum(interior(T))
    Tmin = minimum(interior(T))

    umax = (maximum(abs, interior(u)),
            maximum(abs, interior(v)),
            maximum(abs, interior(w)))

    step_time = 1e-9 * (time_ns() - wall_time[])

    msg = @sprintf("Iter: %d, time: %s, Δt: %s", iteration(sim), prettytime(sim), prettytime(sim.Δt))
    msg *= @sprintf(", max|u|: (%.2e, %.2e, %.2e) m s⁻¹, extrema(T): (%.2f, %.2f) ᵒC, wall time: %s",
                    umax..., Tmax, Tmin, prettytime(step_time))

    @info msg 

    wall_time[] = time_ns()
end

simulation.callbacks[:progress] = Callback(progress, TimeInterval(5days))

# ### Set up output writers
#
# We define output writers to save the simulation data at regular intervals.
# In this case, we save the surface fluxes and surface fields at a relatively high frequency (every day).
# The `indices` keyword argument allows us to save only a slice of the three dimensional variable.
# Below, we use `indices` to save only the values of the variables at the surface, which corresponds to `k = grid.Nz`

outputs = merge(ocean.model.tracers, ocean.model.velocities)
ocean.output_writers[:surface] = JLD2OutputWriter(ocean.model, outputs;
                                                  schedule = TimeInterval(1days),
                                                  filename = "near_global_surface_fields",
                                                  indices = (:, :, grid.Nz),
                                                  with_halos = true,
                                                  overwrite_existing = true,
                                                  array_type = Array{Float32})

# ### Spinning up the simulation
#
# We spin up the simulation with a small-ish time-step to resolve the "initialization shock"
# associated with starting from ECCO2 initial conditions that are both interpolated and also
# satisfy a different dynamical balance than our simulation.

run!(simulation)

# ### Running the simulation for real

# After the initial spin up of 10 days, we can increase the time-step and run for longer.

simulation.stop_time = 60days
simulation.Δt = 10minutes
run!(simulation)

# ## A pretty movie
#
# It's time to make a pretty movie of the simulation. First we load the output we've been saving on
# disk and plot the final snapshot:

u = FieldTimeSeries("near_global_surface_fields.jld2", "u"; backend = OnDisk())
v = FieldTimeSeries("near_global_surface_fields.jld2", "v"; backend = OnDisk())
T = FieldTimeSeries("near_global_surface_fields.jld2", "T"; backend = OnDisk())
e = FieldTimeSeries("near_global_surface_fields.jld2", "e"; backend = OnDisk())

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

s = @at (Center, Center, Nothing) sqrt(un^2 + vn^2) # compute √(u²+v²) and interpolate back to Center, Center
s = Field(s)

sn = @lift begin
    parent(un) .= parent(u[$n])
    parent(vn) .= parent(v[$n])
    compute!(s)
    sn = interior(s)
    sn[land] .= NaN
    view(sn, :, :, 1)
end

title = @lift string("Near-global 1/4 degree ocean simulation after ",
                     prettytime(times[$n] - times[1]))

λ, φ, _ = nodes(T) # T, e, and s all live on the same grid locations

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

save("snapshot.png", fig)
nothing #hide

# ![](snapshot.png)

# And now we make a movie:

record(fig, "near_global_ocean_surface.mp4", 1:Nt, framerate = 8) do nn
    n[] = nn
end
nothing #hide

# ![](near_global_ocean_surface.mp4)
