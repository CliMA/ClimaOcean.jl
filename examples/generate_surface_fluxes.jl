# Surface fluxes from prescribed ocean and atmosphere
#
# ClimaOcean uses bulk formulae to estimate the surface exchange of momentum,
# heat, and water vapor between the atmosphere and the ocean.
#
# This example demonstrates an example of the turbulent surface flux calculations performed in ClimaOcean
# using ECCO2 data for the ocean and JRA55 data for the atmosphere.
#
# For this example, we need ClimaOcean with its DataWrangling modules: ECCO2 and JRA55.
# # We also need Oceananigans for the ImmersedBoundaryGrid and Field utilities, and CairoMakie to plot.

using ClimaOcean
using ClimaOcean.ECCO
using ClimaOcean.JRA55
using ClimaOcean.OceanSimulations
using Oceananigans
using CairoMakie

# # Computing fluxes on the ECCO2 grid
#
# We start by building the ECCO2 grid, using `ECCO_bottom_height` to identify the bottom height.

grid = ECCO_immersed_grid()

fig = Figure()
ax  = Axis(fig[1, 1])
heatmap!(ax, interior(grid.immersed_boundary.bottom_height, :, :, 1))
save("ECCO_continents.png", fig) #hide

# ![](ECCO_continents.png)

# Next, we construct our atmosphere and ocean.
#
# The atmosphere is prescribed, downloaded from the JRA55 dataset.
# It contains:
# - zonal wind `u`
# - meridional wind `v`
# - surface temperature `T`
# - surface relative humidity `q`
# - surface pressure `p`
# - downwelling shortwave radiation
# - downwelling longwave radiation
#
# We invoke the constructor with only the first two time indices, corresponding to 
# January 1st (at 00:00 AM and 03:00 AM).

atmosphere = JRA55_prescribed_atmosphere(1:2; backend = InMemory())
ocean1 = ocean_simulation(grid)
ocean2 = ocean_simulation(grid)

# Now that we have an atmosphere and ocean, we `set!` the ocean temperature and salinity
# to the ECCO2 data by first creating T, S metadata objects,

T_metadata = ECCOMetadata(:temperature)
S_metadata = ECCOMetadata(:salinity)

# Note that if a date is not provided to `ECCOMetadata`, then the default Jan 1st, 1992 is used.
# To copy the ECCO state into `ocean.model`, we use `set!`,

set!(ocean1.model; T=T_metadata, S=S_metadata)
set!(ocean2.model; T=T_metadata, S=S_metadata)

# Finally, we construct a coupled model, which will compute fluxes during construction.
# We omit `sea_ice` so the model is ocean-only, and use the default `Radiation()` that
# uses the two-band shortwave (visible and UV) + longwave (mid and far infrared)
# decomposition of the radiation spectrum.

using ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: SkinTemperature

similarity_theory = SimilarityTheoryTurbulentFluxes(grid, surface_temperature_type = SkinTemperature(κ = 0.01, δ = 1.0))
coupled_model1 = OceanSeaIceModel(ocean1; atmosphere, radiation=Radiation())
coupled_model2 = OceanSeaIceModel(ocean2; atmosphere, radiation=Radiation(), similarity_theory)

Qnet1 = coupled_model1.ocean.model.tracers.T.boundary_conditions.top.condition;
Qnet2 = coupled_model2.ocean.model.tracers.T.boundary_conditions.top.condition;

Ts1 = coupled_model1.fluxes.turbulent.fields.T_surface
Ts2 = coupled_model2.fluxes.turbulent.fields.T_surface

fig = Figure()
ax = Axis(fig[1, 1])
hm = heatmap!((interior(Qnet1, :, :, 1) .- interior(Qnet2, :, :, 1)) .* 1020 .* coupled_model1.fluxes.ocean_heat_capacity, colormap = :balance, colorrange = (-2, 2))
Colorbar(fig[1, 2], hm)
ax = Axis(fig[1, 3])
hm = heatmap!((interior(Ts1, :, :, 1) .- interior(Ts2, :, :, 1)), colormap = :balance, colorrange = (-0.2, 0.2))
Colorbar(fig[1, 4], hm)

fig2 = Figure()
ax  = Axis(fig2[1, 1], title = "iterations for Prescribed temperature") 
hm  = heatmap!(ax, interior(coupled_model1.fluxes.turbulent.fields.x_momentum, :, :, 1), colorrange = (-1, 1)) #, colormap = :balance)
ax  = Axis(fig2[1, 2], title = "iterations for Diagnostic temperature")
hm  = heatmap!(ax, interior(coupled_model2.fluxes.turbulent.fields.x_momentum, :, :, 1), colorrange = (-1, 1)) #, colormap = :balance)
Colorbar(fig2[1, 3], hm)
ax = Axis(fig2[1, 4], title = "Difference in iterations (Prescribed - Diagnostic)")
hm  = heatmap!(ax, interior(coupled_model1.fluxes.turbulent.fields.x_momentum, :, :, 1) .- interior(coupled_model2.fluxes.turbulent.fields.x_momentum, :, :, 1), colorrange = (-1e-2, 1e-2), colormap = :balance)
Colorbar(fig2[1, 5], hm)

# # Now that the surface fluxes are computed, we can extract and visualize them.
# # The turbulent fluxes are stored in `coupled_model.fluxes.turbulent`.

# fluxes  = coupled_model.fluxes.turbulent.fields
# λ, φ, z = nodes(fluxes.sensible_heat)

# fig = Figure(size = (800, 800), fontsize = 15)

# ax = Axis(fig[1, 1], title = "Sensible heat flux (W m⁻²)", ylabel = "Latitude")
# heatmap!(ax, λ, φ, interior(fluxes.sensible_heat, :, :, 1); colormap = :bwr)

# ax = Axis(fig[1, 2], title = "Latent heat flux (W m⁻²)")
# heatmap!(ax, λ, φ, interior(fluxes.latent_heat, :, :, 1); colormap = :bwr)

# ax = Axis(fig[2, 1], title = "Zonal wind stress (N m)", ylabel = "Latitude")
# heatmap!(ax, λ, φ, interior(fluxes.x_momentum, :, :, 1); colormap = :bwr)

# ax = Axis(fig[2, 2], title = "Meridional wind stress (N m)", xlabel = "Longitude")
# heatmap!(ax, λ, φ, interior(fluxes.y_momentum, :, :, 1); colormap = :bwr)

# ax = Axis(fig[3, 1], title = "Water vapor flux (kg m⁻² s⁻¹)", xlabel = "Longitude", ylabel = "Latitude")
# heatmap!(ax, λ, φ, interior(fluxes.water_vapor, :, :, 1); colormap = :bwr)

# save("fluxes.png", fig)
# # ![](fluxes.png)
