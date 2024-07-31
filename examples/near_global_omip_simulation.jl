# # Near-global Ocean simulation
#
# This Julia script sets up and runs a near-global ocean simulation using the Oceananigans.jl and ClimaOcean.jl packages. 
# The simulation spans from 75°S to 75°N with a horizontal resolution of 1/4th of a degree and 40 vertical levels. 
#
# The simulation is then ran for one year and the results are visualized using the CairoMakie.jl package.
#
# ## Initial Setup with Package Imports
#
# The script begins by importing necessary Julia packages for visualization (CairoMakie), 
# ocean modeling (Oceananigans, ClimaOcean), and handling of dates and times (CFTime, Dates). 
# These packages provide the foundational tools for creating the simulation environment, 
# including grid setup, physical processes modeling, and data visualization.

using Printf
using Oceananigans
using Oceananigans.Units
using Oceananigans: architecture, on_architecture
using ClimaOcean
using ClimaOcean.ECCO
using ClimaOcean.OceanSimulations
using ClimaOcean.OceanSeaIceModels

using CFTime
using Dates

# ## Grid Configuration 
#
# We define a near-global grid from 75°S to 75°N with a horizontal resolution of 1/4th of a degree and 40 vertical levels. 
# The grid is created using Oceananigans' `LatitudeLongitudeGrid`. We use an exponential vertical spacing to better resolve the upper ocean layers.
# The total depth of the domain is set to 6000 meters.
# Finally, we specify the architecture to be used for the simulation, which in this case is a GPU.

arch = GPU() 

z_faces = exponential_z_faces(Nz=40, depth=6000)

Nx = 1440
Ny = 600
Nz = length(z_faces) - 1

grid = LatitudeLongitudeGrid(arch; 
                             size = (Nx, Ny, Nz), 
                             halo = (7, 7, 7), 
                             z = z_faces, 
                             longitude = (0, 360),
                             latitude = (-75, 75))

# ## Bathymetry and Immersed Boundary
#
# We retrieve the bathymetry from the ETOPO1 data by ensuring a minimum depth of 10 meters (everything shallower is considered land)
# The `interpolation_passes` parameter specifies the number of passes to interpolate the bathymetry data. The larger the number, 
# the smoother the bathymetry will be. We also remove all connected regions (inland lakes) from the bathymetry data by specifying
# `connected_regions_allowed = 0`.

bottom_height = retrieve_bathymetry(grid; 
                                    minimum_depth = 10,
                                    dir = "./",
                                    interpolation_passes = 20,
                                    connected_regions_allowed = 0)
 
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height)) 

# ## Ocean Model Configuration
#
# To set the ocean simulation we use the `ocean_simulation` function from ClimaOcean.jl. This allows us to build
# an ocean simulation with default parameters and numerics. In this case, the defaults are
# - CATKE turbulence closure for vertical mixing, see [`CATKEVerticalDiffusivity`](@ref)
# - WENO-based advection scheme for momentum in the vector invariant form, see [`WENOVectorInvariant`](@ref)
# - WENO-based advection scheme for tracers, see [`WENO`](@ref)
# - `SplitExplicitFreeSurfaceSolver` with a Courant number of 0.7, see [`SplitExplicitFreeSurface`](@ref)
# - TEOS-10 equation of state, see [`TEOS10EquationOfState`](@ref)
# - Quadratic bottom drag with a drag coefficient of 0.003
#
# The ocean model is then initialized with the ECCO2 temperature and salinity fields for January 1, 1993.

ocean = ocean_simulation(grid) 
model = ocean.model

date  = DateTimeProlepticGregorian(1993, 1, 1)

set!(model, 
     T = ECCOMetadata(:temperature; date),
     S = ECCOMetadata(:salinity;    date))

# ## Prescribed Atmosphere and Radiation
#
# The atmospheric data is prescribed using the JRA55 dataset. The dataset is loaded into memory in 4 snapshots at a time.
# The JRA55 dataset provides atmospheric data such as temperature, humidity, and wind fields to calculate turbulent fluxes
# through bulk formulae, see [`CrossRealmFluxes`](@ref)
#
# The radiation model specified an ocean albedo emissivity to compute the net radiative fluxes. 
# The default ocean albedo is based on the Payne (1982) and depends on cloud cover (computed from
# the maximum possible incident solar radiation divided by the actual incident solar radiation), and the latitude.
# The ocean emissivity is set to 0.97.

backend    = JRA55NetCDFBackend(4) 
atmosphere = JRA55_prescribed_atmosphere(arch; backend)
radiation  = Radiation(arch)

# ## Sea Ice Model 
#
# This simulation includes a simplified representation of an ice cover where the air-sea fluxes are shut down whenever the 
# sea surface temperature is below the freezing point. Only heating fluxes are allowed. This is by no means a sea ice model
# but it allows to include atmosphere-ocean fluxes without having the temperature plummeting to - ∞.  

sea_ice = ClimaOcean.OceanSeaIceModels.MinimumTemperatureSeaIce()

# ## The Coupled Simulation
#
# Finally, we have everything it takes to define the coupled coupled, which includes the ocean, the atmosphere, and the radiation parameters.
# The model is constructed using the `OceanSeaIceModel` constructor.
#
# We can then construct a coupled simulation. In this case we start with a time step of 1 second and run the simulation for 720 days.
# We will eventually increase the time step size as the simulation progresses and the initialization shocks dissipate.
#
# We also define a callback function to monitor the progress of the simulation. This function prints the current time, iteration, time step,
# as well as the maximum velocities and tracers in the domain. The wall time is also printed to monitor the time taken for each iteration.

coupled_model      = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
coupled_simulation = Simulation(coupled_model; Δt=10, stop_time = 720days)

wall_time = [time_ns()]

function progress(sim) 
    ocean = sim.model.ocean
    u, v, w = ocean.model.velocities  
    T = ocean.model.tracers.T

    Tmax = maximum(interior(T))
    Tmin = minimum(interior(T))
    umax = maximum(interior(u)), maximum(interior(v)), maximum(interior(w))
    step_time = 1e-9 * (time_ns() - wall_time[1])

    @info @sprintf("Time: %s, Iteration %d, Δt %s, max(vel): (%.2e, %.2e, %.2e), max(trac): %.2f, %.2f, wtime: %s \n",
                   prettytime(sim.model.clock.time),
                   sim.model.clock.iteration,
                   prettytime(sim.Δt),
                   umax..., Tmax, Tmin, prettytime(step_time))

     wall_time[1] = time_ns()
end

coupled_simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

# ## Set up Output Writers
#
# We define output writers to save the simulation data at regular intervals. 
# In this case we save the surface fluxes, the surface fields, at a fairly high frequency (half a day)
# and snapshots of the fields every 10 days.
#
# We also add a checkpoint writer to save the model state every 60 days, so that we can easily restart the simulation if needed.

fluxes = (u = model.velocities.u.boundary_conditions.top.condition,
          v = model.velocities.v.boundary_conditions.top.condition,
          T = model.tracers.T.boundary_conditions.top.condition,
          S = model.tracers.S.boundary_conditions.top.condition)

output_kwargs = (; overwrite_existing = true, array_type = Array{Float32})

coupled_simulation.output_writers[:fluxes] = JLD2OutputWriter(model, fluxes;
                                                              schedule = TimeInterval(0.5days),
                                                              overwrite_existing = true,
                                                              filename = "surface_fluxes",
                                                              output_kwargs...)

coupled_simulation.output_writers[:surface] = JLD2OutputWriter(model, merge(model.tracers, model.velocities);
                                                               schedule = TimeInterval(0.5days),
                                                               filename = "surface",
                                                               indices = (:, :, grid.Nz),
                                                               output_kwargs...)

coupled_simulation.output_writers[:snapshots] = JLD2OutputWriter(model, merge(model.tracers, model.velocities);
                                                                 schedule = TimeInterval(10days),
                                                                 filename = "snapshots",
                                                                 output_kwargs...)

coupled_simulation.output_writers[:checkpoint] = Checkpointer(model;
                                                              schedule = TimeInterval(60days),
                                                              overwrite_existing = true,
                                                              prefix = "checkpoint")

# ## Warming up the simulation!
#
# As an initial condition we have interpolated ECCO tracer fields on our custom grid.
# Most likely the bathymetry of the original ECCO data is different than our grid, so the initialization of the velocity
# field might lead to shocks if we use a large time step.
#
# Therefore, we warm up with a small time step to ensure that the interpolated initial conditions adapt to the model numerics and
# parameterization without crashing. 30 days of integration with a maximum time step of 1.5 minutes should be enough to dissipate
# spurious initialization shocks.

ocean.stop_time = 30days
wizard = TimeStepWizard(; cfl = 0.1, max_Δt = 90, max_change = 1.1)
ocean.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

run!(coupled_simulation)

# ## Running the simulation
#
# Now that the simulation has been warmed up, we can run it for the full year.
# We increase the maximum time step size to 10 minutes and let the simulation run for 720 days.

ocean.stop_time = 720days
wizard = TimeStepWizard(; cfl = 0.25, max_Δt = 10minutes, max_change = 1.1)
ocean.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

# ## Visualizing the Results
# 
# after the simulation has finished, we can visualize the results.

using CairoMakie

u, v, w = model.velocities
T, S, e = model.tracers

using Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField

ζ = VerticalVorticityField(model)
s = Field(sqrt(u^2 + v^2))

compute!(ζ)
compute!(s)

ζ = on_architecture(CPU(), ζ)
s = on_architecture(CPU(), s)
T = on_architecture(CPU(), T)
e = on_architecture(CPU(), e)

fig = Figure(size = (1000, 800))

ax = Axis(fig[1, 1], title = "Vertical vorticity [s⁻¹]")
heatmap!(ax, interior(ζ, :, :, grid.Nz), colorrange = (-4e-5, 4e-5), colormap = :bwr)

ax = Axis(fig[1, 2], title = "Surface speed [ms⁻¹]")
heatmap!(ax, interior(s, :, :, grid.Nz), colorrange = (0, 0.5), colormap = :deep)

ax = Axis(fig[2, 1], title = "Surface Temperature [Cᵒ]")
heatmap!(ax, interior(T, :, :, grid.Nz), colorrange = (-1, 30), colormap = :magma)

ax = Axis(fig[2, 2], title = "Turbulent Kinetic Energy [m²s⁻²]")
heatmap!(ax, interior(e, :, :, grid.Nz), colorrange = (0, 1e-3), colormap = :solar)

save("near_global_ocean_surface.png", fig)
nothing #hide

# ![](near_global_ocean_surface.png)

