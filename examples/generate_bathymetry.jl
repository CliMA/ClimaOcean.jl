using CairoMakie
using Oceananigans
using ClimaOcean.Bathymetry: regrid_bathymetry

#####
##### Quarter degree near-global grid
#####

grid = LatitudeLongitudeGrid(CPU();
                             size = (4 * 360, 4 * 160, 1),
                             latitude = (-80, 80),
                             longitude = (-180, 180),
                             z = (0, 1),
                             halo = (4, 4, 4))

h = regrid_bathymetry(grid, height_above_water=1, minimum_depth=5)

λ, φ, z = nodes(h)

land = interior(h) .>= 0
interior(h)[land] .= NaN

fig = Figure(size=(2400, 1200))
ax = Axis(fig[1, 1])
heatmap!(ax, λ, φ, interior(h, :, :, 1), nan_color=:white, colorrange=(-5000, 0))

λp = -112.45
φp = 42.86
text = "😻"
text!(ax, λp, φp; text, fontsize=30)

display(fig)

#####
##### Regional Mediterranean grid 
#####

# 1/25th degree resolution
Nλ = 25 * 55
Nφ = 25 * 25

grid = LatitudeLongitudeGrid(CPU();
                             size = (Nλ, Nφ, 1),
                             latitude = (25, 50),
                             longitude = (-10, 45),
                             z = (0, 1),
                             halo = (4, 4, 4))

h_smooth  = regrid_bathymetry(grid, height_above_water=1, minimum_depth=10, interpolation_passes = 40)
h_rough   = regrid_bathymetry(grid, height_above_water=1, minimum_depth=10, interpolation_passes = 1)
h_nolakes = regrid_bathymetry(grid, height_above_water=1, minimum_depth=10, connected_regions_allowed = 0)

λ, φ, z = nodes(h_smooth)

land_smooth = interior(h_smooth) .>= 0
interior(h_smooth)[land_smooth] .= NaN
land_rough = interior(h_rough) .>= 0
interior(h_rough)[land_rough] .= NaN
land_nolakes = interior(h_nolakes) .>= 0
interior(h_nolakes)[land_nolakes] .= NaN

fig = Figure(resolution=(2400, 800))
ax = Axis(fig[1, 1])
heatmap!(ax, λ, φ, interior(h_smooth,  :, :, 1), nan_color=:white) #, colorrange=(-5000, 0))
ax = Axis(fig[1, 2])
heatmap!(ax, λ, φ, interior(h_rough,   :, :, 1), nan_color=:white) #, colorrange=(-5000, 0))
ax = Axis(fig[1, 3])
heatmap!(ax, λ, φ, interior(h_nolakes, :, :, 1), nan_color=:white) #, colorrange=(-5000, 0))

display(fig)
