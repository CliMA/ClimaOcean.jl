# # Near-global ocean simulation
#
# This example sets up and runs a near-global ocean simulation using the Oceananigans.jl and
# ClimaOcean.jl packages. The simulation covers latitudes from 75°S to 75°N with a horizontal
# resolution of 1/4 degree and 40 vertical levels.
#
# The simulation's results are visualized using the CairoMakie.jl package.
#
# ## Initial setup with package imports
#
# We begin by importing the necessary Julia packages for visualization (CairoMakie),
# ocean modeling (Oceananigans, ClimaOcean), and handling dates and times (CFTime, Dates).
# These packages provide the foundational tools for setting up the simulation environment,
# including grid setup, physical processes modeling, and data visualization.

using Printf
using Oceananigans
using Oceananigans.Units
using ClimaOcean
using CairoMakie

using CFTime
using Dates

# ### Grid configuration 
#
# We define a global grid with a horizontal resolution of 1/4 degree and 30 vertical levels.
# The grid is a `LatitudeLongitudeGrid` spanning latitudes from 75°S to 75°N.
# We use an exponential vertical spacing to better resolve the upper ocean layers.
# The total depth of the domain is set to 6000 meters.
# Finally, we specify the architecture for the simulation, which in this case is a GPU.

arch = GPU() 

z_faces = exponential_z_faces(Nz=30, depth=4000)

Nx = 1440
Ny = 600
Nz = length(z_faces) - 1

grid = LatitudeLongitudeGrid(arch;
                             size = (Nx, Ny, Nz),
                             halo = (7, 7, 7),
                             z = z_faces,
                             latitude  = (-75, 75),
                             longitude = (0, 360))

# ### Bathymetry and immersed boundary
#
# We use `regrid_bathymetry` to derive the bottom height from ETOPO1 data.
# To smooth the interpolated data we use 5 interpolation passes. We also fill in all
# sminor enclosed basins but the 3 largest `major_basins` as well as reasons
# that are shallower than `minimum_depth = 10`.

bathymetry_path = ClimaOcean.Bathymetry.download_bathymetry_cache
rm(bathymetry_path, force=true)

# try twice
tries = 3
_try = 0
built_bathymetry = false

while (_try < tries) || !built_bathymetry
    try
        bottom_height = regrid_bathymetry(grid; 
                                          minimum_depth = 10,
                                          interpolation_passes = 5,
                                          major_basins = 3)

        built_bathymetry = true
    catch err
        @warn "Building the bathymetry failed, trying again..." exception=err
        rm(bathymetry_path, force=true)
        _try += 1
    end
end
        
# Let's see what the bathymetry looks like:

zb = grid.immersed_boundary.bottom_height

fig, ax, hm = heatmap(zb, colormap=:deep, colorrange=(-6000, 0))
cb = Colorbar(fig[0, 1], hm, label="Bottom height (m)", vertical=false)
hidedecorations!(ax)
save("bathymetry.png", fig)
nothing #hide

# ![](bathymetry.png)

# ### Ocean model configuration
#
# We build our ocean model using `ocean_simulation`,

ocean = ocean_simulation(grid)

# which uses the default `ocean.model`,

ocean.model

# We initialize the ocean model to ECCO2 temperature and salinity for January 1, 1993.

date = DateTimeProlepticGregorian(1993, 1, 1)
set!(ocean.model, T=ECCOMetadata(:temperature; dates=date), S=ECCOMetadata(:salinity; dates=date))

# ### Prescribed atmosphere and radiation
#
# Next we build a prescribed atmosphere state and radiation model,
# which will drive the development of the ocean simulation.
# We use the default `Radiation` model,

## The radiation model specifies an ocean albedo emissivity to compute the net radiative
## fluxes. The default ocean albedo is based on Payne (1982) and depends on cloud cover
## (calculated from the ratio of maximum possible incident solar radiation to actual
## incident solar radiation) and latitude. The ocean emissivity is set to 0.97.

radiation = Radiation(arch)

# The atmospheric data is prescribed using the JRA55 dataset.
# The number of snapshots that are loaded into memory is determined by
# the `backend`

# into memory in 41 snapshots at a time. The JRA55 dataset provides atmospheric
# data such as temperature, humidity, and wind fields to calculate turbulent fluxes
# using bulk formulae, see [`CrossRealmFluxes`](@ref).

atmosphere = JRA55_prescribed_atmosphere(arch; backend=JRA55NetCDFBackend(41))

# ### Sea ice model 
#
# This simulation includes a simplified representation of ice cover where the
# air-sea fluxes are shut down whenever the sea surface temperature is below
# the freezing point,

sea_ice = ClimaOcean.OceanSeaIceModels.MinimumTemperatureSeaIce()

# ## The coupled simulation

# Next we assemble the ocean, sea ice, atmosphere, and radiation
# into a coupled model,

coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

# We then create a coupled simulation, starting with a time step of 10 seconds
# and running the simulation for 10 days.

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
# The `indices` keyword argument allows us to save down a slice at the surface, which is located at `k = grid.Nz`

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
# We spin up the simulation with a very short time-step to resolve the "initialization shock"
# associated with starting from ECCO initial conditions that are both interpolated and also
# satisfy a different dynamical balance than our simulation.

run!(simulation)

# ### Running the simulation
#
# Now that the simulation has spun up, we can run it for the full 60 days.
# We increase the maximum time step size to 10 minutes and let the simulation run for 60 days.

# simulation.stop_time = 60days
# simulation.Δt = 10minutes
# run!(simulation)
# nothing #hide

# ## A pretty movie
#
# It's time to make a pretty movie of the simulation. First we plot a snapshot:

u = FieldTimeSeries("near_global_surface_fields.jld2", "u"; backend = OnDisk())
v = FieldTimeSeries("near_global_surface_fields.jld2", "v"; backend = OnDisk())
T = FieldTimeSeries("near_global_surface_fields.jld2", "T"; backend = OnDisk())
e = FieldTimeSeries("near_global_surface_fields.jld2", "e"; backend = OnDisk())

times = u.times
Nt = length(times)

n = Observable(Nt)

Tn = @lift interior(T[$n], :, :, 1)
en = @lift interior(e[$n], :, :, 1)

un = Field{Face, Center, Nothing}(u.grid)
vn = Field{Center, Face, Nothing}(v.grid)
s = Field(sqrt(un^2 + vn^2))

sn = @lift begin
    parent(un) .= parent(u[$n])
    parent(vn) .= parent(v[$n])
    compute!(s)
end

fig = Figure(size = (800, 1200))

axs = Axis(fig[1, 1], xlabel="Longitude (deg)", ylabel="Latitude (deg)")
axT = Axis(fig[2, 1], xlabel="Longitude (deg)", ylabel="Latitude (deg)")
axe = Axis(fig[3, 1], xlabel="Longitude (deg)", ylabel="Latitude (deg)")

hm = heatmap!(axs, sn, colorrange = (0, 0.5), colormap = :deep)
Colorbar(fig[1, 2], hm, label = "Surface speed (m s⁻¹)")

hm = heatmap!(axT, Tn, colorrange = (-1, 30), colormap = :magma)
Colorbar(fig[2, 2], hm, label = "Surface Temperature (ᵒC)")

hm = heatmap!(axe, en, colorrange = (0, 1e-3), colormap = :solar)
Colorbar(fig[3, 2], hm, label = "Turbulent Kinetic Energy (m² s⁻²)")
save("snapshot.png", fig)
nothing #hide

# ![](snapshot.png)

# And now a movie:

record(fig, "near_global_ocean_surface.mp4", 1:Nt, framerate = 8) do nn
    n[] = nn
end
nothing #hide

# ![](near_global_ocean_surface.mp4)
