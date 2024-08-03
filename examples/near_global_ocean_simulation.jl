# # Near-global Ocean Simulation
#
# This Julia script sets up and runs a near-global ocean simulation using the Oceananigans.jl and ClimaOcean.jl packages. 
# The simulation covers latitudes from 75°S to 75°N with a horizontal resolution of 1/4 degree and 40 vertical levels. 
#
# The simulation runs for one year, and the results are visualized using the CairoMakie.jl package.
#
# ## Initial Setup with Package Imports
#
# The script begins by importing the necessary Julia packages for visualization (CairoMakie), 
# ocean modeling (Oceananigans, ClimaOcean), and handling dates and times (CFTime, Dates). 
# These packages provide the foundational tools for setting up the simulation environment, 
# including grid setup, physical processes modeling, and data visualization.

using Printf
using Oceananigans
using Oceananigans.Units
using Oceananigans: architecture, on_architecture
using ClimaOcean
using ClimaOcean.ECCO
using ClimaOcean.OceanSimulations
using ClimaOcean.OceanSeaIceModels
using CairoMakie
using OrthogonalSphericalShellGrids

using CFTime
using Dates

# ### Grid Configuration 
#
# We define a global grid with a horizontal resolution of 1/4 degree and 40 vertical levels. 
# The grid is a `LatitudeLongitudeGrid` capped at 75°S to 75°N.
# We use an exponential vertical spacing to better resolve the upper ocean layers. The total depth of the domain is set to 6000 meters.
# Finally, we specify the architecture for the simulation, which in this case is a GPU.

arch = GPU() 

z_faces = exponential_z_faces(Nz=40, depth=6000)

Nx = 1440
Ny = 600
Nz = length(z_faces) - 1

grid = LatitudeLongitudeGrid(arch; 
                             size = (Nx, Ny, Nz), 
                             halo = (7, 7, 7), 
                             z = z_faces, 
                             latitude  = (-75, 75),
                             longitude = (0, 360))

# ### Bathymetry and Immersed Boundary
#
# We retrieve the bathymetry from the ETOPO1 data, ensuring a minimum depth of 10 meters (depths shallower than this are considered land).
# The `interpolation_passes` parameter specifies the number of passes to interpolate the bathymetry data. A larger number 
# results in a smoother bathymetry. We also remove all connected regions (such as inland lakes) from the bathymetry data by specifying
# `connected_regions_allowed = 2` (Mediterrean sea an North sea in addition to the ocean) 

bottom_height = retrieve_bathymetry(grid; 
                                    minimum_depth = 10,
                                    interpolation_passes = 5,
                                    connected_regions_allowed = 2)
 
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height)) 

bathymetry = deepcopy(Array(interior(bottom_height, :, :, 1)))
bathymetry[bathymetry .>= 0] .= NaN

fig = Figure(size = (800, 400))
axis = Axis(fig[1, 1], title = "Bathymetry [m]")
hm = heatmap!(axis, bathymetry, colormap = :deep, colorrange = (-6000, 0))
cb = Colorbar(fig[1, 2], hm)
hidedecorations!(axis)

save("bathymetry.png", fig)
nothing #hide

# ![](bathymetry.png)

# ### Ocean Model Configuration
#
# To configure the ocean simulation, we use the `ocean_simulation` function from ClimaOcean.jl. This function allows us to build
# an ocean simulation with default parameters and numerics. The defaults include:
# - CATKE turbulence closure for vertical mixing, see [`CATKEVerticalDiffusivity`](@ref)
# - WENO-based advection scheme for momentum in the vector invariant form, see [`WENOVectorInvariant`](@ref)
# - WENO-based advection scheme for tracers, see [`WENO`](@ref)
# - `SplitExplicitFreeSurfaceSolver` with 75 substeps, see [`SplitExplicitFreeSurface`](@ref)
# - TEOS-10 equation of state, see [`TEOS10EquationOfState`](@ref)
# - Quadratic bottom drag with a drag coefficient of 0.003
#
# The ocean model is then initialized with the ECCO2 temperature and salinity fields for January 1, 1993.

free_surface = SplitExplicitFreeSurface(grid; substeps = 75)

ocean = ocean_simulation(grid; free_surface) 
model = ocean.model

date  = DateTimeProlepticGregorian(1993, 1, 1)

set!(model, 
     T = ECCOMetadata(:temperature; date),
     S = ECCOMetadata(:salinity;    date))

nothing #hide

# ### Prescribed Atmosphere and Radiation
#
# The atmospheric data is prescribed using the JRA55 dataset, which is loaded into memory in 4 snapshots at a time.
# The JRA55 dataset provides atmospheric data such as temperature, humidity, and wind fields to calculate turbulent fluxes
# using bulk formulae, see [`CrossRealmFluxes`](@ref).
#
# The radiation model specifies an ocean albedo emissivity to compute the net radiative fluxes. 
# The default ocean albedo is based on Payne (1982) and depends on cloud cover (calculated from
# the ratio of maximum possible incident solar radiation to actual incident solar radiation) and latitude.
# The ocean emissivity is set to 0.97.

backend    = JRA55NetCDFBackend(4) 
atmosphere = JRA55_prescribed_atmosphere(arch; backend)
radiation  = Radiation(arch)
nothing #hide

# ### Sea Ice Model 
#
# This simulation includes a simplified representation of ice cover where the air-sea fluxes are shut down whenever the 
# sea surface temperature is below the freezing point. Only heating fluxes are allowed. This is not a full sea ice model,
# but it prevents the temperature from dropping excessively low by including atmosphere-ocean fluxes.

sea_ice = ClimaOcean.OceanSeaIceModels.MinimumTemperatureSeaIce()
nothing #hide

# ## The Coupled Simulation
#
# Finally, we define the coupled model, which includes the ocean, atmosphere, and radiation parameters.
# The model is constructed using the `OceanSeaIceModel` constructor.
#
# We then create a coupled simulation, starting with a time step of 10 seconds and running the simulation for 10 days.
# We will eventually increase the time step size and end time as the simulation progresses and initialization shocks dissipate.
#
# We also define a callback function to monitor the simulation's progress. This function prints the current time, iteration, time step,
# as well as the maximum velocities and tracers in the domain. The wall time is also printed to monitor the time taken for each iteration.

coupled_model      = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
coupled_simulation = Simulation(coupled_model; Δt=10, stop_time = 10days)

wall_time = [time_ns()]

function progress(sim) 
    ocean = sim.model.ocean
    u, v, w = ocean.model.velocities  
    T = ocean.model.tracers.T

    Tmax = maximum(interior(T))
    Tmin = minimum(interior(T))
    umax = maximum(abs, interior(u)), maximum(abs, interior(v)), maximum(abs, interior(w))
    step_time = 1e-9 * (time_ns() - wall_time[1])

    @info @sprintf("Time: %s, Iteration %d, Δt %s, max(vel): (%.2e, %.2e, %.2e), max(T): %.2f, min(T): %.2f, wtime: %s \n",
                   prettytime(ocean.model.clock.time),
                   ocean.model.clock.iteration,
                   prettytime(ocean.Δt),
                   umax..., Tmax, Tmin, prettytime(step_time))

     wall_time[1] = time_ns()
end

coupled_simulation.callbacks[:progress] = Callback(progress, IterationInterval(500))
nothing #hide

# ### Set up Output Writers
#
# We define output writers to save the simulation data at regular intervals. 
# In this case, we save the surface fluxes and surface fields at a relatively high frequency (every half day).

fluxes = (u = model.velocities.u.boundary_conditions.top.condition,
          v = model.velocities.v.boundary_conditions.top.condition,
          T = model.tracers.T.boundary_conditions.top.condition,
          S = model.tracers.S.boundary_conditions.top.condition)

output_kwargs = (; overwrite_existing = true, array_type = Array{Float32})

ocean.output_writers[:fluxes] = JLD2OutputWriter(model, fluxes;
                                                 schedule = TimeInterval(1days),
                                                 overwrite_existing = true,
                                                 filename = "surface_fluxes",
                                                 output_kwargs...)

ocean.output_writers[:surface] = JLD2OutputWriter(model, merge(model.tracers, model.velocities);
                                                  schedule = TimeInterval(1days),
                                                  filename = "surface",
                                                  indices = (:, :, grid.Nz),
                                                  output_kwargs...)
nothing #hide

# ### Warming Up the Simulation
#
# As an initial condition, we have interpolated ECCO tracer fields onto our custom grid.
# The bathymetry of the original ECCO data may differ from our grid, so the initialization of the velocity
# field might cause shocks if a large time step is used.
#
# Therefore, we warm up with a small time step to ensure that the interpolated initial conditions adapt to the model numerics and
# parameterization without causing instability. A 30-day integration with a maximum time step of 1.5 minutes should be sufficient to dissipate
# spurious initialization shocks.

ocean.stop_time = 10days
wizard = TimeStepWizard(; cfl = 0.1, max_Δt = 90, max_change = 1.1)
ocean.callbacks[:wizard] = Callback(wizard, IterationInterval(10))
run!(coupled_simulation)
nothing #hide

# ### Running the simulation
#
# Now that the simulation has been warmed up, we can run it for the full two years.
# We increase the maximum time step size to 10 minutes and let the simulation run for 720 days.

ocean.stop_time = 60days
coupled_simulation.stop_time = 60days
wizard = TimeStepWizard(; cfl = 0.25, max_Δt = 15minutes, max_change = 1.1)
ocean.callbacks[:wizard] = Callback(wizard, IterationInterval(10))
run!(coupled_simulation)
nothing #hide

# ## Visualizing the Results
# 
# The simulation has finished, let's visualize the results.
# In this section we pull up the saved data and create visualizations using the CairoMakie.jl package.
# In particular, we generate a video of the evolution of surface fields:
# zonal velocity (u), meridional velocity (v), surface temperature (T), and turbulent kinetic energy (e)

u = FieldTimeSeries("surface.jld2", "u"; backend = OnDisk())
v = FieldTimeSeries("surface.jld2", "v"; backend = OnDisk())
T = FieldTimeSeries("surface.jld2", "T"; backend = OnDisk())
e = FieldTimeSeries("surface.jld2", "e"; backend = OnDisk())

Nt = length(u.times)

iter = Observable(Nt)

ui = @lift(interior(u[$iter], :, :, 1))
vi = @lift(interior(v[$iter], :, :, 1))
Ti = @lift(interior(T[$iter], :, :, 1))
ei = @lift(interior(e[$iter], :, :, 1))

fig = Figure(size = (1000, 750))

ax = Axis(fig[1, 1], title = "Zonal velocity [ms⁻¹]")
heatmap!(ax, ui, colorrange = (-0.5, 0.5), colormap = :bwr)
hidedecorations!(ax)

ax = Axis(fig[1, 2], title = "Meridional velocity [ms⁻¹]")
heatmap!(ax, vi, colorrange = (-0.5, 0.5), colormap = :bwr)
hidedecorations!(ax)

ax = Axis(fig[2, 1], title = "Surface Temperature [Cᵒ]")
heatmap!(ax, Ti, colorrange = (-1, 30), colormap = :magma)
hidedecorations!(ax)

ax = Axis(fig[2, 2], title = "Turbulent Kinetic Energy [m²s⁻²]")
heatmap!(ax, ei, colorrange = (0, 1e-3), colormap = :solar)
hidedecorations!(ax)

CairoMakie.record(fig, "near_global_ocean_surface.mp4", 1:length(u.times), framerate = 8) do i
     @info "Generating frame $i of $(length(u.times))"
     iter[] = i
end
nothing #hide

# ![](near_global_ocean_surface.mp4)

