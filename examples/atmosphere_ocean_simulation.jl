# # Global coupled atmosphere and ocean--sea ice simulation
#
# This example configures a global ocean--sea ice simulation at 1ᵒ horizontal resolution with
# realistic bathymetry and a few closures including the "Gent-McWilliams" `IsopycnalSkewSymmetricDiffusivity`.
# The atmosphere is represented by a SpeecdyWeather model at T63 resolution (approximately 1.875ᵒ).
# and initialized by temperature, salinity, sea ice concentration, and sea ice thickness
# from the ECCO state estimate.
#
# For this example, we need Oceananigans.HydrostaticFreeSurfaceModel (the ocean), ClimaSeaIce.SeaIceModel (the sea ice) and 
# SpeedyWeather.PrimitiveWetModel (the atmosphere), coupled and orchestrated by ClimaOcean.OceanSeaIceModel (the coupled system).
# The XESMF.jl package is used to regrid fields between the atmosphere and ocean--sea ice components.

using Oceananigans, SpeedyWeather, XESMF, ClimaOcean
using NCDatasets, CairoMakie
using Oceananigans.Units
using Printf, Statistics, Dates

# ## Ocean and sea-ice model configuration
# The ocean and sea-ice are a simplified version of the "one_degree_simulation" example
# https://clima.github.io/ClimaOceanDocumentation/dev/literated/one_degree_simulation/
# The first step is to create the grid with realistic bathymetry.

Nx = 360 
Ny = 180 
Nz = 10  

r_faces = ExponentialDiscretization(Nz, -2000, 0)
grid    = TripolarGrid(Oceananigans.CPU(); size=(Nx, Ny, Nz), z=r_faces, halo=(6, 6, 5))

# Regridding the bathymetry...

bottom_height = regrid_bathymetry(grid; major_basins=1, interpolation_passes=15)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)

# Now we can specify the numerical details and closures for the ocean simulation.

momentum_advection = WENOVectorInvariant(order=5)
tracer_advection   = WENO(order=5)
free_surface = SplitExplicitFreeSurface(grid; substeps=40)

catke_closure   = ClimaOcean.OceanSimulations.default_ocean_closure()
viscous_closure = Oceananigans.TurbulenceClosures.HorizontalScalarBiharmonicDiffusivity(ν=1e10)
eddy_closure    = Oceananigans.TurbulenceClosures.IsopycnalSkewSymmetricDiffusivity(κ_skew=1e3, κ_symmetric=1e3)
closures        = (catke_closure, eddy_closure, viscous_closure, VerticalScalarDiffusivity(ν=1e-4))

# The ocean simulation, complete with initial conditions for temperature and salinity from ECCO.

ocean = ocean_simulation(grid; 
                         momentum_advection,
                         tracer_advection,
                         free_surface,
                         timestepper = :SplitRungeKutta3,
                         closure = closures)

Oceananigans.set!(ocean.model, T=Metadatum(:temperature, dataset=ECCO4Monthly()), 
                               S=Metadatum(:salinity,    dataset=ECCO4Monthly()))

# The sea-ice simulation, complete with initial conditions for sea-ice thickness and concentration from ECCO.

sea_ice = sea_ice_simulation(grid, ocean; advection=WENO(order=7))

Oceananigans.set!(sea_ice.model, h=Metadatum(:sea_ice_thickness, dataset=ECCO4Monthly()), 
                                 ℵ=Metadatum(:sea_ice_concentration, dataset=ECCO4Monthly()))

# ## Atmosphere model configuration
# The atmosphere is provided by SpeedyWeather.jl. Here we configure a T63L8 model with a 3 hour output interval.
# The `atmosphere_simulation` function takes care of building an atmosphere model with appropriate 
# hooks for ClimaOcean to compute intercomponent fluxes. We also set the output interval to 3 hours.

spectral_grid = SpectralGrid(trunc=63, nlayers=8, Grid=FullClenshawGrid)
atmosphere = atmosphere_simulation(spectral_grid; output=true)
atmosphere.model.output.output_dt = Hour(3)

# ## The coupled model
# Now we can build the coupled model. We need to specify the time step for the coupled model.
# We decide to step the global model every 2 atmosphere time steps. (i.e. the ocean and the 
# sea-ice will be stepped every two atmosphere time steps).

Δt = 2 * convert(eltype(grid), atmosphere.model.time_stepping.Δt_sec)

# We build the complete model. Since radiation is idealized in this example, we set the emissivities to zero.

radiation = Radiation(ocean_emissivity=0.0, sea_ice_emissivity=0.0)
earth_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
earth = Oceananigans.Simulation(earth_model; Δt, stop_time=60days)

# ## Running the simulation
# We can now run the simulation. We add a callback to monitor the progress of the simulation
# and write outputs to disk every 3 hours.

wall_time = Ref(time_ns())

function progress(sim)
    sea_ice = sim.model.sea_ice
    ocean   = sim.model.ocean
    hmax  = maximum(sea_ice.model.ice_thickness)
    ℵmax  = maximum(sea_ice.model.ice_concentration)
    uimax = maximum(abs, sea_ice.model.velocities.u)
    vimax = maximum(abs, sea_ice.model.velocities.v)
    uomax = maximum(abs, ocean.model.velocities.u)
    vomax = maximum(abs, ocean.model.velocities.v)

    step_time = 1e-9 * (time_ns() - wall_time[])

    msg1 = @sprintf("time: %s, iteration: %d, Δt: %s, ", prettytime(sim), iteration(sim), prettytime(sim.Δt))
    msg2 = @sprintf("max(h): %.2e m, max(ℵ): %.2e ", hmax, ℵmax)
    msg3 = @sprintf("max uᵢ: (%.2f, %.2f) m s⁻¹, ", uimax, vimax)
    msg4 = @sprintf("max uₒ: (%.2f, %.2f) m s⁻¹, ", uomax, vomax)
    msg5 = @sprintf("wall time: %s \n", prettytime(step_time))

    @info msg1 * msg2 * msg3 * msg4 * msg5

    wall_time[] = time_ns()

     return nothing
end

outputs = merge(ocean.model.velocities, ocean.model.tracers)
sea_ice_fields = merge(sea_ice.model.velocities, sea_ice.model.dynamics.auxiliaries.fields, 
                       (; h=sea_ice.model.ice_thickness, ℵ=sea_ice.model.ice_concentration))

ocean.output_writers[:free_surf] = JLD2Writer(ocean.model, (; η=ocean.model.free_surface.η);
                                              overwrite_existing=true,
                                              schedule=TimeInterval(3600 * 3),
                                              filename="ocean_free_surface.jld2")

ocean.output_writers[:surface] = JLD2Writer(ocean.model, outputs;
                                            overwrite_existing=true,
                                            schedule=TimeInterval(3600 * 3),
                                            filename="ocean_surface_fields.jld2",
                                            indices=(:, :, grid.Nz))

sea_ice.output_writers[:fields] = JLD2Writer(sea_ice.model, sea_ice_fields;
                                             overwrite_existing=true,
                                             schedule=TimeInterval(3600 * 3),
                                             filename="sea_ice_fields.jld2")

Qcao = earth.model.interfaces.atmosphere_ocean_interface.fluxes.sensible_heat
Qvao = earth.model.interfaces.atmosphere_ocean_interface.fluxes.latent_heat
τxao = earth.model.interfaces.atmosphere_ocean_interface.fluxes.x_momentum
τyao = earth.model.interfaces.atmosphere_ocean_interface.fluxes.y_momentum
Qcai = earth.model.interfaces.atmosphere_sea_ice_interface.fluxes.sensible_heat
Qvai = earth.model.interfaces.atmosphere_sea_ice_interface.fluxes.latent_heat
τxai = earth.model.interfaces.atmosphere_sea_ice_interface.fluxes.x_momentum
τyai = earth.model.interfaces.atmosphere_sea_ice_interface.fluxes.y_momentum
Qoi  = earth.model.interfaces.net_fluxes.sea_ice_bottom.heat
Soi  = earth.model.interfaces.sea_ice_ocean_interface.fluxes.salt
fluxes = (; Qcao, Qvao, τxao, τyao, Qcai, Qvai, τxai, τyai, Qoi, Soi)

earth.output_writers[:fluxes] = JLD2Writer(earth.model.ocean.model, fluxes;
                                           overwrite_existing=true,
                                           schedule=TimeInterval(3600 * 3),
                                           filename="intercomponent_fluxes.jld2")

add_callback!(earth, progress, IterationInterval(100))

Oceananigans.run!(earth)

# ## Visualizing the results
# We can visualize some of the results. Here we plot the surface speeds in the atmosphere, ocean, and sea-ice
# as well as the 2m temperature in the atmosphere, the sea surface temperature, and the sensible and latent heat
# fluxes at the atmosphere-ocean interface. SpeedyWeather outputs are stored in a NetCDF file located in the `run_0001` folder,
# while ocean and sea-ice outputs are stored in JLD2 files that can be read by Oceananigans.jl using the `FieldTimeSeries` type.

SWO = Dataset("run_0001/output.nc")

Ta = reverse(SWO["temp"][:, :, 8, :], dims=2)
ua = reverse(SWO["u"][:, :, 8, :],    dims=2)
va = reverse(SWO["v"][:, :, 8, :],    dims=2)
sp = sqrt.(ua.^2 + va.^2)

SST = FieldTimeSeries("ocean_surface_fields.jld2", "T")
SSU = FieldTimeSeries("ocean_surface_fields.jld2", "u")
SSV = FieldTimeSeries("ocean_surface_fields.jld2", "v")

SIU = FieldTimeSeries("sea_ice_fields.jld2", "u")
SIV = FieldTimeSeries("sea_ice_fields.jld2", "v")
SIA = FieldTimeSeries("sea_ice_fields.jld2", "ℵ")

Qcao = FieldTimeSeries("intercomponent_fluxes.jld2", "Qcao")
Qvao = FieldTimeSeries("intercomponent_fluxes.jld2", "Qvao")

uotmp = Oceananigans.Field{Face, Center, Nothing}(SST.grid)
votmp = Oceananigans.Field{Center, Face, Nothing}(SST.grid)

uitmp = Oceananigans.Field{Face,   Center, Nothing}(SST.grid)
vitmp = Oceananigans.Field{Center, Face,   Nothing}(SST.grid)
atmp  = Oceananigans.Field{Center, Center, Nothing}(SST.grid)

sotmp = Oceananigans.Field(sqrt(uotmp^2 + votmp^2))
sitmp = Oceananigans.Field(sqrt(uitmp^2 + vitmp^2) * atmp)

iter = Observable(1)
san = @lift sp[:, :, $iter]
son  = @lift begin
    Oceananigans.set!(uotmp, SSU[$iter * 2])
    Oceananigans.set!(votmp, SSV[$iter * 2])
    Oceananigans.compute!(sotmp)
    Oceananigans.interior(sotmp, :, :, 1)
end

ssn  = @lift begin
    Oceananigans.set!(uitmp, SIU[$iter * 2])
    Oceananigans.set!(vitmp, SIV[$iter * 2])
    Oceananigans.set!(atmp,  SIA[$iter * 2])
    Oceananigans.compute!(sitmp)
    Oceananigans.interior(sitmp, :, :, 1)
end

fig = Figure(size = (1200, 400))
ax2 = Axis(fig[1, 1], title = "Surface speed, atmosphere (m/s)")
hm2 = heatmap!(ax2, san; colormap = :deep)
ax1 = Axis(fig[1, 2], title = "Surface speed, ocean (m/s)")
hm = heatmap!(ax1, son; colormap = :deep)
ax3 = Axis(fig[1, 3], title = "Surface speed, sea-ice (m/s)")
hm = heatmap!(ax3, ssn; colormap = :deep)

record(fig, "surface_speeds.mp4", 1:length(sp[1, 1, :]), framerate = 15) do i
    iter[] = i
end
nothing #hide

# ![](surface_speeds.mp4)

Tan = @lift Ta[:, :, $iter]
Ton = @lift interior(SST[$iter * 2], :, :, 1)
Qcn = @lift interior(Qcao[$iter * 2], :, :, 1)
Qvn = @lift interior(Qvao[$iter * 2], :, :, 1)

fig = Figure(size = (1200, 800))
ax1 = Axis(fig[1, 1], title = "2m Temperature, atmosphere (K)")
hm = heatmap!(ax1, Tan; colormap = :plasma)
ax2 = Axis(fig[1, 2], title = "Sea Surface Temperature (C)")
hm2 = heatmap!(ax2, Ton; colormap = :plasma)
ax3 = Axis(fig[2, 1], title = "Sensible heat flux (W/m²)")
hm3 = heatmap!(ax3, Qcn; colormap = :balance, colorrange = (-200, 200))
ax4 = Axis(fig[2, 2], title = "Latent heat flux (W/m²)")
hm4 = heatmap!(ax4, Qvn; colormap = :balance, colorrange = (-200, 200))

record(fig, "surface_temperature_and_heat_flux.mp4", 1:length(sp[1, 1, :])) do i
    iter[] = i
end
nothing #hide

# ![](surface_temperature_and_heat_flux.mp4)
