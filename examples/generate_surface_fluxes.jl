# # Surface fluxes from prescribed ocean and atmosphere
#
# ClimaOcean uses bulk formulae to estimate the surface exchange of momentum,
# heat, and water vapor between the atmosphere and the ocean.
#
# This example demonstrates an example of the turbulent surface flux calculations performed in ClimaOcean
# using ECCO2 data for the ocean and JRA55 data for the atmosphere.
#
# For this example, we need ClimaOcean with its DataWrangling modules: ECCO2 and JRA55.
# We also need Oceananigans for the ImmersedBoundaryGrid and Field utilities, and CairoMakie to plot.

using ClimaOcean
using ClimaOcean.ECCO
using ClimaOcean.JRA55
using ClimaOcean.OceanSimulations
using Oceananigans
using CairoMakie

# # Computing fluxes on the ECCO2 grid
#
# We start by building the ECCO2 grid, using `ECCO_mask`
# to identify missing cells.

mask = ECCO_mask()
grid = mask.grid
grid = ImmersedBoundaryGrid(grid, GridFittedBoundary(mask))

fig = Figure()
ax  = Axis(fig[1, 1])
heatmap!(ax, interior(grid.immersed_boundary.mask, :, :, grid.Nz))
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
ocean = ocean_simulation(grid)

# Now that we have an atmosphere and ocean, we `set!` the ocean temperature and salinity
# to the ECCO2 data by first creating T, S metadata objects,

T_metadata = ECCOMetadata(:temperature)
S_metadata = ECCOMetadata(:salinity)

# Note that if a date is not provided to `ECCOMetadata`, then the default Jan 1st, 1992 is used.
# To copy the ECCO state into `ocean.model`, we use `set!`,

set!(ocean.model; T=T_metadata, S=S_metadata)

# Finally, we construct a coupled model, which will compute fluxes during construction.
# We omit `sea_ice` so the model is ocean-only, and use the default `Radiation()` that
# uses the two-band shortwave (visible and UV) + longwave (mid and far infrared)
# decomposition of the radiation spectrum.

coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation=Radiation())

# Now that the surface fluxes are computed, we can extract and visualize them.
# The turbulent fluxes are stored in `coupled_model.fluxes.turbulent`.

fluxes = coupled_model.fluxes.turbulent.fields

fig = Figure(size = (800, 400))

ax = Axis(fig[1, 1], title = "Sensible heat flux (W m⁻²)")
heatmap!(ax, fluxes.sensible_heat; colormap = :bwr)

ax = Axis(fig[1, 2], title = "Latent heat flux (W m⁻²)")
heatmap!(ax, fluxes.latent_heat; colormap = :bwr)

ax = Axis(fig[2, 1], title = "Zonal wind stress (N m)")
heatmap!(ax, fluxes.x_momentum; colormap = :bwr)

ax = Axis(fig[2, 2], title = "Meridional wind stress (N m)")
heatmap!(ax, fluxes.y_momentum; colormap = :bwr)

ax = Axis(fig[3, 1], title = "Water vapor flux (kg m⁻² s⁻¹)")
heatmap!(ax, Mv; colormap = :bwr)

save("fluxes.png", fig)
# ![](fluxes.png)
