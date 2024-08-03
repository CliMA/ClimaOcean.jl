# # Generate bathymetry data for the Mediterranean Sea
#
# This script shows how to configure an Immersed boundary grid with realistic bathymetry using ClimaOcean.jl
# by generating the bathymetry data for the Mediterranean Sea.
#
# For this example, we need Oceananigans for the LatitudeLongitudeGrid and Field utilities, 
# ClimaOcean to donwload and regrid the bathymetry, and CairoMakie to visualize the grid.

using ClimaOcean
using Oceananigans
using CairoMakie

# We start by defining a gridded domain for the Mediterranean Sea using the `LatitudeLongitudeGrid` from Oceananigans. 
# To have a reasonable resolution, we set a grid size to 1/25ᵒ of a degree in both latitude and longitude.
#
# The Mediterranean sea is positioned roughly between 25ᵒ and 50ᵒ latitude and -10ᵒ and 45ᵒ longitude.

Nλ = 25 * 55
Nφ = 25 * 25

grid = LatitudeLongitudeGrid(size = (Nλ, Nφ, 1),
                             latitude = (25, 50),
                             longitude = (-10, 45),
                             z = (0, 1))

# Next, we generate the bathymetry data for the Mediterranean Sea using the `regrid_bathymetry` function from ClimaOcean.
# The function downloads the bathymetry data from the ETOPO1 dataset, regrids it to the provided grid, and returns the bathymetry field.
# The three different regidding procedures here show the effect of different parameters on the generated bathymetry.
#
# - `h_rough`  shows the output of the function with default parameters, which means only one interpolation passes and no restrictions on connected regions.
# - `h_smooth` shows the output of the function with 40 interpolation passes, which results in a smoother bathymetry.
# - `h_nolakes` shows the output of the function with `connected_regions_allowed = 0`, which means that the function does not allow connected regions in the bathymetry 
#   (e.g., lakes) and fills them with land.

h_rough   = regrid_bathymetry(grid)
h_smooth  = regrid_bathymetry(grid; interpolation_passes = 40)
h_nolakes = regrid_bathymetry(grid; connected_regions_allowed = 0)

# Finally, we visualize the generated bathymetry data for the Mediterranean Sea using CairoMakie.

λ, φ, z = nodes(h_smooth)

land_smooth = interior(h_smooth) .>= 0
interior(h_smooth)[land_smooth] .= NaN
land_rough = interior(h_rough) .>= 0
interior(h_rough)[land_rough] .= NaN
land_nolakes = interior(h_nolakes) .>= 0
interior(h_nolakes)[land_nolakes] .= NaN

fig = Figure(resolution=(1200, 400))
ax = Axis(fig[1, 1])
heatmap!(ax, λ, φ, interior(h_smooth,  :, :, 1), nan_color=:white) 

ax = Axis(fig[1, 2])
heatmap!(ax, λ, φ, interior(h_rough,   :, :, 1), nan_color=:white) 

ax = Axis(fig[1, 3])
heatmap!(ax, λ, φ, interior(h_nolakes, :, :, 1), nan_color=:white) 

save("different_bottom_heights.png", fig)
nothing #hide

# ![](different_bottom_heights.png) 