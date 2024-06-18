# # Calculating surface fluxes with an ocean and an atmosphere
#
# ClimaOcean uses bulk formulae to estimate surface exchange of momentum,
# heat, and water vapor between the atmosphere and the ocean.
#
# This script shows an example of the turbulent surface flux calculations performed in ClimaOcean
# using ECCO2 data for the ocean and JRA55 data for the atmosphere
#
# For this example we need ClimaOcean with it's DataWrangling modules: ECCO2 and JRA55.
# We also need Oceananigans for the ImmersedBoundaryGrid and Field utilies and CairoMakie to plot

using ClimaOcean
using ClimaOcean.ECCO2
using ClimaOcean.JRA55
using ClimaOcean.OceanSimulations
using Oceananigans
using CairoMakie

# We start by defining a grid. The ECCO2 grid is a good starting point.
# The ECCO2 grid is not "immersed" by default, but we can use the ECCO mask
# to define the "bathymetry" of the ECCO fields.
# `ecco2_center_mask` produces a field with `true` values where data is missing (aka, in immersed cells). 
# We can use this mask as an immersed boundary for our grid
# Let's create the grid and visualize the mask

mask = ecco2_center_mask()
grid = mask.grid
grid = ImmersedBoundaryGrid(grid, GridFittedBoundary(mask))

fig = Figure()
ax  = Axis(fig[1, 1])
heatmap!(ax, interior(grid.immersed_boundary.mask, :, :, grid.Nz))

save("ecco_continents.png", fig)
nothing #hide

# ![](ecco_continents.png)

# We now construct our atmosphere and our ocean. 
# The atmosphere is prescribed, downloaded from the JRA55 dataset.
# It contains 
# - zonal wind `u`
# - meridional wind `v`
# - surface temperature `T`
# - surface relative humidity `q`
# - surface pressure `p`
# - downwelling shortwave radiation
# - downwelling longwave radiation
#
# We invoke the constructor with only the first two time indices, as they correspond to 
# January 1st (at 00:00AM and 03:00AM). 
# By passing the ECCO grid we automatically interpolate the atmospheric data on the grid.
# Note that this is recommended only for small simulations (in terms of grid size).
# By omitting the grid, the interpolation will be done on the fly.
#
# We construct the ocean simulation without bothering for advection, closures or coriolis given
# that we will not time-step the ocean but only use it to construct the fluxes

atmosphere  = JRA55_prescribed_atmosphere(1:2; backend = InMemory(), grid = grid.underlying_grid)

ocean = ocean_simulation(grid; momentum_advection = nothing,
                                 tracer_advection = nothing,
                                          closure = nothing,
                                         coriolis = nothing)

# Now that we have an atmosphere, and a container for the ocean, we need to populate
# out ocean with initial conditions. To so this we can use the ECCO2 dataset by 
# `set!`ting the model with the `ECCO2Metadata`. If no date is specified,
# the fields corresponding to the 1st of January 1992 (the first available date in
# ECCO2) is used.
# This command will download the fields on the local machine.

set!(ocean.model;
      T = ECCO2Metadata(:temperature), 
      S = ECCO2Metadata(:salinity), 
      u = ECCO2Metadata(:u_velocity), 
      v = ECCO2Metadata(:v_velocity))

# The final step is to construct a coupled model.
# The coupled model requires an ocean, which we have just constructed and initialized,
# an atmosphere, which we downloaded from the JRA55 dataset, the sea_ice simulation 
# (in this case we neglect the sea ice by defining `sea_ice = nothing`), and a 
# radiation model.
# The default radiation model is a model that assumes only two spectral bands: a shortwave and
# a longwave band. 
# By constructing the coupled model, the `update_state!` function, which calculates the fluxes
# will be triggered.

radiation = Radiation()
sea_ice   = nothing
coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

# Now that the surface fluxes are computed we can extract them and visualize them.
# The turbulent fluxes are stored in `coupled_model.fluxes.turbulent`.
# 
# Qs = coupled_model.fluxes.turbulent.fields.sensible_heat : the sensible heat flux 
# Ql = coupled_model.fluxes.turbulent.fields.latent_heat   : the latent heat flux 
# τx = coupled_model.fluxes.turbulent.fields.x_momentum    : the zonal wind stress
# τy = coupled_model.fluxes.turbulent.fields.y_momentum    : the meridional wind stress
# Mv = coupled_model.fluxes.turbulent.fields.water_vapor   : evaporation
#
# They are 3D Fields with one point in the vertical. To extract the data we use the 
# `interior` functionality from Oceananigans

turbulent_fluxes = coupled_model.fluxes.turbulent.fields

Qs = interior(turbulent_fluxes.sensible_heat, :, :, 1)
Ql = interior(turbulent_fluxes.latent_heat,   :, :, 1)
τx = interior(turbulent_fluxes.x_momentum,    :, :, 1)
τy = interior(turbulent_fluxes.y_momentum,    :, :, 1)
Mv = interior(turbulent_fluxes.water_vapor,   :, :, 1)
nothing

fig = Figure()

ax = Axis(fig[1, 1], title = "Sensible heat flux")
heatmap!(ax, Qs; colormap = :bwr)

ax = Axis(fig[1, 2], title = "Latent heat flux")
heatmap!(ax, Ql; colormap = :bwr)

ax = Axis(fig[2, 1], title = "Zonal wind stress")
heatmap!(ax, τx; colormap = :bwr)

ax = Axis(fig[2, 2], title = "Meridional wind stress")
heatmap!(ax, τy; colormap = :bwr)

ax = Axis(fig[3, 1], title = "Evaporation")
heatmap!(ax, Mv; colormap = :bwr)

save("turbulent_fluxes.png", fig)
nothing #hide

# ![](turbulent_fluxes.png)
