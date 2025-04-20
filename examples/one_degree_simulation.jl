# # One-degree global ocean simulation
#
# This example configures a global ocean--sea ice simulation at 1ᵒ horizontal resolution with
# realistic bathymetry, few closures. The simulation is forced by JRA55 atmospheric reanalysis
# and initialized by temperature and salinity from ECCO2 state estimate.
#
# For this example, we need Oceananigans, ClimaOcean, and
# CairoMakie to visualize the simulation. Also we need CFTime and Dates for date handling.

using ClimaOcean
using ClimaOcean.ECCO
using Oceananigans
using Oceananigans.Units
using Oceananigans.OrthogonalSphericalShellGrids
using Dates
using Printf
using ClimaOcean.ECCO: download_dataset

# ### Download necessary files to run the code

# ### ECCO files

start_date = DateTime(1993, 1, 1) 
stop_date = DateTime(1993, 12, 1) 
dates = range(start_date, step=Month(1), stop=stop_date)
ecco_temperature = Metadata(:temperature; dates, dataset=ECCO4Monthly())
ecco_salinity = Metadata(:salinity; dates, dataset=ECCO4Monthly())

download_dataset(ecco_temperature)
download_dataset(ecco_salinity)

# ### Grid and Bathymetry

arch = GPU()
Nx = 180
Ny = 85
Nz = 40

z = exponential_z_faces(; Nz, depth=4000, h=34)
underlying_grid = TripolarGrid(arch; size = (Nx, Ny, Nz), halo = (5, 5, 4), z)
# underlying_grid = LatitudeLongitudeGrid(arch; size = (Nx, Ny, Nz), halo = (5, 5, 4), z, longitude=(0, 360), latitude=(-80, 80))

## 75 interpolation passes smooth the bathymetry near Florida so that the Gulf Stream is able to flow:
bottom_height = regrid_bathymetry(underlying_grid;
                                  minimum_depth = 20,
                                  interpolation_passes = 100,
                                  major_basins = 2)

# For this bathymetry at this horizontal resolution we need to manually open the Gibraltar strait.
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_height); active_cells_map=true)

# ### Restoring
#
# We include temperature and salinity surface restoring to ECCO data near the surface.
restoring_rate  = 1 / 10days
mask = LinearlyTaperedPolarMask(southern=(-80, -70), northern=(70, 90), z=(-100, 0))
FT = ECCORestoring(ecco_temperature, arch; mask, rate=restoring_rate)
FS = ECCORestoring(ecco_salinity, arch; mask, rate=restoring_rate)
forcing = (T=FT, S=FS)

# ### Closures
#
# We include a Gent-McWilliams isopycnal diffusivity as a parameterization for the mesoscale
# eddy fluxes. For vertical mixing at the upper-ocean boundary layer we include the CATKE
# parameterization. We also include some explicit horizontal diffusivity.

eddy_closure = Oceananigans.TurbulenceClosures.IsopycnalSkewSymmetricDiffusivity(κ_skew=1e3, κ_symmetric=1e3)
vertical_mixing = ClimaOcean.OceanSimulations.default_ocean_closure()
horizontal_viscosity = HorizontalScalarDiffusivity(ν=2000)
closure = (eddy_closure, horizontal_viscosity, vertical_mixing)

# ### Ocean simulation
# Now we bring everything together to construct the ocean simulation.
# We use a split-explicit timestepping with 30 substeps for the barotropic
# mode.

free_surface = SplitExplicitFreeSurface(grid; substeps=50)
momentum_advection = VectorInvariant()
tracer_advection = WENO(order=5)

ocean = ocean_simulation(grid; momentum_advection, tracer_advection, closure, forcing, free_surface)

# ### Initial condition

# We initialize the ocean from the ECCO state estimate.

set!(ocean.model, T=ecco_temperature[1], S=ecco_salinity[1])

# ### Atmospheric forcing

# We force the simulation with an JRA55-do atmospheric reanalysis.
radiation  = Radiation(arch)
atmosphere = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(20))

# ### Coupled simulation

# Now we are ready to build the coupled ocean--sea ice model and bring everything
# together into a `simulation`.

# We use a relatively short time step initially and only run for a few days to
# avoid numerical instabilities from the initial "shock" of the adjustment of the
# flow fields.

coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)
simulation = Simulation(coupled_model; Δt=5minutes, stop_time=10days)

# ### A progress messenger
#
# We write a function that prints out a helpful progress message while the simulation runs.

wall_time = Ref(time_ns())

function progress(sim)
    ocean = sim.model.ocean
    u, v, w = ocean.model.velocities
    T = ocean.model.tracers.T
    e = ocean.model.tracers.e
    Tmin, Tmax = minimum(T), maximum(T)
    emin, emax = minimum(e), maximum(e)
    umax = (maximum(abs, u), maximum(abs, v), maximum(abs, w))
            

    step_time = 1e-9 * (time_ns() - wall_time[])

    msg1 = @sprintf("time: %s, iteration: %d, Δt: %s, ", prettytime(sim), iteration(sim), prettytime(sim.Δt))
    msg2 = @sprintf("max|u|: (%.2e, %.2e, %.2e) m s⁻¹, ", umax...)
    msg3 = @sprintf("extrema(T): (%.2f, %.2f) ᵒC, ", Tmin, Tmax)
    msg4 = @sprintf("extrema(e): (%.2f, %.2f) m2 s-2, ", emin, emax)
    msg5 = @sprintf("wall time: %s \n", prettytime(step_time))

    @info msg1 * msg2 * msg3 * msg4

     wall_time[] = time_ns()

     return nothing
end

# And add it as a callback to the simulation.
add_callback!(simulation, progress, IterationInterval(100))

# ### Output
#
# We are almost there! We need to save some output. Below we choose to save _only surface_
# fields using the `indices` keyword argument. We save all velocity and tracer components.
# Note, that besides temperature and salinity, the CATKE vertical mixing parameterization
# also uses a prognostic turbulent kinetic energy, `e`, to diagnose the vertical mixing length.

outputs = merge(ocean.model.tracers, ocean.model.velocities)
ocean.output_writers[:surface] = JLD2Writer(ocean.model, outputs;
                                            schedule = TimeInterval(5days),
                                            filename = "one_degree_surface_fields",
                                            indices = (:, :, grid.Nz),
                                            overwrite_existing = true)

# ### Ready to run

# We are ready to press the big red button and run the simulation.

# After we run for a short time (here we set up the simulation with `stop_time = 10days`),
# we increase the timestep and run for longer.

run!(simulation)

simulation.Δt = 20minutes
simulation.stop_time = 360days
run!(simulation)

# ### A pretty movie
#
# We load the saved output and make a pretty movie of the simulation. First we plot a snapshot:
using CairoMakie

u = FieldTimeSeries("one_degree_surface_fields.jld2", "u"; backend = OnDisk())
v = FieldTimeSeries("one_degree_surface_fields.jld2", "v"; backend = OnDisk())
T = FieldTimeSeries("one_degree_surface_fields.jld2", "T"; backend = OnDisk())
e = FieldTimeSeries("one_degree_surface_fields.jld2", "e"; backend = OnDisk())

times = u.times
Nt = length(times)
n = Observable(Nt)

# We create a land mask and use it to fill land points with `NaN`s.
land = interior(T.grid.immersed_boundary.bottom_height) .≥ 0

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

# We compute the surface speed.
un = Field{Face, Center, Nothing}(u.grid)
vn = Field{Center, Face, Nothing}(v.grid)
s = Field(sqrt(un^2 + vn^2))

sn = @lift begin
    parent(un) .= parent(u[$n])
    parent(vn) .= parent(v[$n])
    compute!(s)
    sn = interior(s)
    sn[land] .= NaN
    view(sn, :, :, 1)
end

# Finally, we plot a snapshot of the surface speed, temperature, and the turbulent
# eddy kinetic energy from the CATKE vertical mixing parameterization.
fig = Figure(size = (800, 1200))

title = @lift string("Global 1ᵒ ocean simulation after ", prettytime(times[$n] - times[1]))

axs = Axis(fig[1, 1])
axT = Axis(fig[2, 1])
axe = Axis(fig[3, 1])

hm = heatmap!(axs, sn, colorrange = (0, 0.5), colormap = :deep, nan_color=:lightgray)
Colorbar(fig[1, 2], hm, label = "Surface speed (m s⁻¹)")

hm = heatmap!(axT, Tn, colorrange = (-1, 30), colormap = :magma, nan_color=:lightgray)
Colorbar(fig[2, 2], hm, label = "Surface Temperature (ᵒC)")

hm = heatmap!(axe, en, colorrange = (0, 1e-3), colormap = :solar, nan_color=:lightgray)
Colorbar(fig[3, 2], hm, label = "Turbulent Kinetic Energy (m² s⁻²)")

for ax in (axs, axT, axe)
    hidedecorations!(ax)
end

Label(fig[0, :], title)

save("global_snapshot.png", fig)
nothing #hide

# ![](global_snapshot.png)

# And now a movie:

record(fig, "one_degree_global_ocean_surface.mp4", 1:Nt, framerate = 8) do nn
    n[] = nn
end
nothing #hide

# ![](one_degree_global_ocean_surface.mp4)
