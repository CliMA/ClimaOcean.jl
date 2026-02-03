# # Global climate simulation
#
# This example configures a global ocean--sea ice simulation at 1.5ᵒ horizontal resolution with
# realistic bathymetry and a few closures including the "Gent-McWilliams" `IsopycnalSkewSymmetricDiffusivity`.
# The atmosphere is represented by a 4-layer [SpeedyWeather](https://github.com/SpeedyWeather/SpeedyWeather.jl)
# simulation on the T63 spectral grid (this grid has approximately 1.875ᵒ resolution).
#
# The atmosphere is initialized with the [jablonowski2006baroclinic](@citet) initial conditions,
# which consist of a zonal wind centered at mid-latitudes and higher altitudes and a temperature
# profile that is baroclinically unstable. The surface pressure is adjusted by orography for
# approximately globally constant mean-sea level pressure. The initial specific humidity is
# calculated from temperature for a constant relative humidity everywhere.
# The ocean and sea ice are initialized by ocean temperature, salinity, sea ice concentration,
# and sea ice thickness from the ECCO state estimate.
#
# For this example, we need `Oceananigans.HydrostaticFreeSurfaceModel` (the ocean), `ClimaSeaIce.SeaIceModel` (the sea ice) and
# `SpeedyWeather.PrimitiveWetModel` (the atmosphere). All these are coupled and orchestrated by the `ClimaOcean.OceanSeaIceModel`
# (the coupled system).
#
# The XESMF.jl package is used to regrid fields between the atmosphere and ocean--sea ice components.

using Oceananigans, SpeedyWeather, XESMF, ClimaOcean
using NCDatasets, CairoMakie
using Oceananigans.Units
using Printf, Statistics, Dates

# ## Ocean and sea-ice model configuration
# The ocean and sea-ice are a simplified versions of the [one-degree ocean-sea ice example](@ref one-degree-ocean-seaice).
#
# The first step is to create the grid with realistic bathymetry.

Nx = 240
Ny = 120
Nz = 10
z = ExponentialDiscretization(Nz, -2000, 0)
grid = TripolarGrid(Oceananigans.CPU(); size=(Nx, Ny, Nz), z, halo=(6, 6, 5))
nothing #hide

# We regrid the bathymetry.

bottom_height = regrid_bathymetry(grid; major_basins=1, interpolation_passes=15)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)
nothing #hide

# Now we can specify the numerical details and the closures for the ocean simulation.

momentum_advection   = VectorInvariant()
tracer_advection     = WENO(order=5)
free_surface         = SplitExplicitFreeSurface(grid; substeps=40)
catke_closure        = ClimaOcean.Oceans.default_ocean_closure()
eddy_closure         = Oceananigans.TurbulenceClosures.IsopycnalSkewSymmetricDiffusivity(κ_skew=1e3, κ_symmetric=1e3)
viscous_closure      = Oceananigans.TurbulenceClosures.HorizontalScalarBiharmonicDiffusivity(ν=1e12)
closures             = (catke_closure, eddy_closure, viscous_closure)
nothing #hide

# The ocean simulation, complete with initial conditions for temperature and salinity from ECCO.

ocean = ocean_simulation(grid;
                         momentum_advection,
                         tracer_advection,
                         free_surface,
                         closure = closures)

Oceananigans.set!(ocean.model, T=Metadatum(:temperature, dataset=ECCO4Monthly()),
                               S=Metadatum(:salinity,    dataset=ECCO4Monthly()))

# The sea-ice simulation, complete with initial conditions for sea-ice thickness and sea-ice concentration from ECCO.

sea_ice = sea_ice_simulation(grid, ocean; advection=WENO(order=5))

Oceananigans.set!(sea_ice.model, h=Metadatum(:sea_ice_thickness, dataset=ECCO4Monthly()),
                                 ℵ=Metadatum(:sea_ice_concentration, dataset=ECCO4Monthly()))

# ## Atmosphere model configuration
# The atmosphere is provided by SpeedyWeather.jl. Here, we configure a T63L4 model with a 3-hour output interval.
# The `atmosphere_simulation` function takes care of building an atmosphere model with appropriate
# hooks so that ClimaOcean can compute inter-component fluxes.
nlayers = 4
spectral_grid = SpeedyWeather.SpectralGrid(; trunc=63, nlayers, Grid=FullClenshawGrid)
atmosphere = atmosphere_simulation(spectral_grid, output=true)

# The atmosphere model already includes some initial conditions as described above:

atmosphere.model.initial_conditions

# We use a three hour time-step:

atmosphere.model.output.output_dt = Hour(3)
nothing #hide

# ## The coupled model
# We are now ready to blend everything together.
# We need to specify the time step for the coupled model.
# We decide to step the global model every 2 atmosphere time steps (i.e., the ocean and the
# sea-ice are stepped every two atmosphere time steps).

Δt = 2 * convert(eltype(grid), atmosphere.model.time_stepping.Δt_sec)
nothing #hide

# We build the complete coupled `earth_model` and the coupled simulation.
# Since radiation is idealized in this example, we set the emissivities to zero.

radiation = Radiation(ocean_emissivity=0.0, sea_ice_emissivity=0.0)
earth_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

# ## Building and running the simulation
#
# We are ready to build and run the coupled simulation.
# But before we do, we add callbacks to write outputs to disk every 3 hours.

earth = Oceananigans.Simulation(earth_model; Δt, stop_time=30days)
outputs = merge(ocean.model.velocities, ocean.model.tracers)
sea_ice_fields = merge(sea_ice.model.velocities, sea_ice.model.dynamics.auxiliaries.fields,
                       (; h=sea_ice.model.ice_thickness, ℵ=sea_ice.model.ice_concentration))

ocean.output_writers[:free_surf] = JLD2Writer(ocean.model, (; η=ocean.model.free_surface.displacement);
                                              overwrite_existing=true,
                                              schedule=TimeInterval(3hours),
                                              including = [:grid],
                                              filename="ocean_free_surface.jld2")

ocean.output_writers[:surface] = JLD2Writer(ocean.model, outputs;
                                            overwrite_existing=true,
                                            schedule=TimeInterval(3hours),
                                            including = [:grid],
                                            filename="ocean_surface_fields.jld2",
                                            indices=(:, :, grid.Nz))

sea_ice.output_writers[:fields] = JLD2Writer(sea_ice.model, sea_ice_fields;
                                             overwrite_existing=true,
                                             schedule=TimeInterval(3hours),
                                             including = [:grid],
                                             filename="sea_ice_fields.jld2")

Qcao = earth.model.interfaces.atmosphere_ocean_interface.fluxes.sensible_heat
Qvao = earth.model.interfaces.atmosphere_ocean_interface.fluxes.latent_heat
τxao = earth.model.interfaces.atmosphere_ocean_interface.fluxes.x_momentum
τyao = earth.model.interfaces.atmosphere_ocean_interface.fluxes.y_momentum
Qcai = earth.model.interfaces.atmosphere_sea_ice_interface.fluxes.sensible_heat
Qvai = earth.model.interfaces.atmosphere_sea_ice_interface.fluxes.latent_heat
τxai = earth.model.interfaces.atmosphere_sea_ice_interface.fluxes.x_momentum
τyai = earth.model.interfaces.atmosphere_sea_ice_interface.fluxes.y_momentum
Qoi  = earth.model.interfaces.net_fluxes.sea_ice.bottom.heat
Soi  = earth.model.interfaces.sea_ice_ocean_interface.fluxes.salt
fluxes = (; Qcao, Qvao, τxao, τyao, Qcai, Qvai, τxai, τyai, Qoi, Soi)

ocean.output_writers[:fluxes] = JLD2Writer(earth.model.ocean.model, fluxes;
                                           overwrite_existing=true,
                                           schedule=TimeInterval(3hours),
                                           including = [:grid],
                                           filename="intercomponent_fluxes.jld2")

# We also add a callback function that prints out a helpful progress message while the simulation runs.

wall_time = Ref(time_ns())

function progress(sim)
    atmos = sim.model.atmosphere
    ocean = sim.model.ocean

    ua, va     = atmos.diagnostic_variables.dynamics.u_mean_grid, atmos.diagnostic_variables.dynamics.v_mean_grid
    uo, vo, wo = ocean.model.velocities

    uamax = (maximum(abs, ua), maximum(abs, va))
    uomax = (maximum(abs, uo), maximum(abs, vo), maximum(abs, wo))

    step_time = 1e-9 * (time_ns() - wall_time[])

    msg1 = @sprintf("time: %s, iter: %d", prettytime(sim), iteration(sim))
    msg2 = @sprintf(", max|ua|: (%.1e, %.1e) m s⁻¹", uamax...)
    msg3 = @sprintf(", max|uo|: (%.1e, %.1e, %.1e) m s⁻¹", uomax...)
    msg4 = @sprintf(", wall time: %s \n", prettytime(step_time))

    @info msg1 * msg2 * msg3 * msg4

    wall_time[] = time_ns()

     return nothing
end

add_callback!(earth, progress, TimeInterval(2days))

# Let's run the coupled model!

Oceananigans.run!(earth)

# ## Visualizing the results
# We can visualize some of the results. Here, we plot the surface speeds in the atmosphere, ocean, and sea-ice
# as well as the 2-meter temperature in the atmosphere, the sea surface temperature, and the sensible and latent heat
# fluxes at the atmosphere-ocean interface. SpeedyWeather outputs are stored in a NetCDF file located in the `run_0001` folder,
# while ocean and sea-ice outputs are stored in JLD2 files that can be read by Oceananigans.jl using the `FieldTimeSeries` type.

SWO = Dataset("run_0001/output.nc")

Ta = reverse(SWO["temp"][:, :, nlayers, :], dims=2)
ua = reverse(SWO["u"][:, :, nlayers, :],    dims=2)
va = reverse(SWO["v"][:, :, nlayers, :],    dims=2)
sp = sqrt.(ua.^2 + va.^2)

SST = FieldTimeSeries("ocean_surface_fields.jld2", "T")
SSU = FieldTimeSeries("ocean_surface_fields.jld2", "u")
SSV = FieldTimeSeries("ocean_surface_fields.jld2", "v")

SIU = FieldTimeSeries("sea_ice_fields.jld2", "u")
SIV = FieldTimeSeries("sea_ice_fields.jld2", "v")
SIA = FieldTimeSeries("sea_ice_fields.jld2", "ℵ")

Qcao = FieldTimeSeries("intercomponent_fluxes.jld2", "Qcao")
Qvao = FieldTimeSeries("intercomponent_fluxes.jld2", "Qvao")

Nt = min(length(sp[1, 1, :]), length(Qcao))
times = Qcao.times

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
    Oceananigans.set!(uotmp, SSU[$iter])
    Oceananigans.set!(votmp, SSV[$iter])
    Oceananigans.compute!(sotmp)
    Oceananigans.interior(sotmp, :, :, 1)
end
ssn  = @lift begin
    Oceananigans.set!(uitmp, SIU[$iter])
    Oceananigans.set!(vitmp, SIV[$iter])
    Oceananigans.set!(atmp,  SIA[$iter])
    Oceananigans.compute!(sitmp)
    Oceananigans.interior(sitmp, :, :, 1)
end

fig = Figure(size = (1000, 1500))

ax1 = Axis(fig[1, 1], title = "Surface speed, atmosphere")
ax2 = Axis(fig[2, 1], title = "Surface speed, ocean")
ax3 = Axis(fig[3, 1], title = "Surface speed, sea-ice")

hm1 = heatmap!(ax1, san; colormap = :deep,  nan_color=:lightgray, colorrange = (0, 35))
hm2 = heatmap!(ax2, son; colormap = :magma, nan_color=:lightgray, colorrange = (0, 0.6))
hm3 = heatmap!(ax3, ssn; colormap = :ice,   nan_color=:lightgray, colorrange = (0, 0.6))

Colorbar(fig[1, 2], hm1, label="(m s⁻¹)")
Colorbar(fig[2, 2], hm2, label="(m s⁻¹)")
Colorbar(fig[3, 2], hm3, label="(m s⁻¹)")

for ax in (ax1, ax2, ax3)
    hidedecorations!(ax)
end

title = @lift string(prettytime(times[$iter] - times[1]))
Label(fig[0, :], title, fontsize=18)

record(fig, "surface_speeds.mp4", 1:Nt, framerate=8) do i
    iter[] = i
end
nothing #hide

# ![](surface_speeds.mp4)

Tan = @lift Ta[:, :, $iter]
Ton = @lift interior(SST[$iter], :, :, 1)
Qcn = @lift interior(Qcao[$iter], :, :, 1)
Qvn = @lift interior(Qvao[$iter], :, :, 1)

fig = Figure(size = (1000, 2000))

ax1 = Axis(fig[1, 1], title = "2m Temperature, atmosphere")
ax2 = Axis(fig[2, 1], title = "Sea Surface Temperature")
ax3 = Axis(fig[3, 1], title = "Sensible heat flux")
ax4 = Axis(fig[4, 1], title = "Latent heat flux")

hm1 = heatmap!(ax1, Tan; colormap = :plasma, nan_color=:lightgray, colorrange = (-45, 30))
hm2 = heatmap!(ax2, Ton; colormap = :plasma, nan_color=:lightgray, colorrange = (-2, 32))
hm3 = heatmap!(ax3, Qcn; colormap = :balance, colorrange = (-200, 200),  nan_color=:lightgray)
hm4 = heatmap!(ax4, Qvn; colormap = :balance, colorrange = (-200, 200),  nan_color=:lightgray)

Colorbar(fig[1, 2], hm1, label="(ᵒC)")
Colorbar(fig[2, 2], hm2, label="(ᵒC)")
Colorbar(fig[3, 2], hm3, label="(W/m²)")
Colorbar(fig[4, 2], hm4, label="(W/m²)")

for ax in (ax1, ax2, ax3, ax4)
    hidedecorations!(ax)
end

title = @lift string(prettytime(times[$iter] - times[1]))
Label(fig[0, :], title, fontsize=18)

record(fig, "surface_temperature_and_heat_flux.mp4", 1:Nt, framerate=8) do i
    iter[] = i
end
nothing #hide

# ![](surface_temperature_and_heat_flux.mp4)
