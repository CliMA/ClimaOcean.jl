using GLMakie
using Oceananigans
using ClimaOcean.Bathymetry: regrid_bathymetry
using OrthogonalSphericalShellGrids

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

h_smooth = regrid_bathymetry(grid, height_above_water=1, interpolation_passes = 40)
h_rough  = regrid_bathymetry(grid, height_above_water=1, interpolation_passes = 1)

λ, φ, z = nodes(h_smooth)

land_smooth = interior(h_smooth) .> 0
interior(h_smooth)[land_smooth] .= NaN
land_rough = interior(h_rough) .> 0
interior(h_rough)[land_rough] .= NaN

fig = Figure(resolution=(2400, 800))
ax = Axis(fig[1, 1])
heatmap!(ax, λ, φ, interior(h_smooth, :, :, 1), nan_color=:white) #, colorrange=(-5000, 0))
ax = Axis(fig[1, 2])
heatmap!(ax, λ, φ, interior(h_rough, :, :, 1), nan_color=:white) #, colorrange=(-5000, 0))

display(fig)

#####
##### Quarter degree Global ocean grid with the north-pole singularity in Russia at λ = 230
#####

grid = WarpedLatitudeLongitudeGrid(CPU();
                                   size = (4 * 360, 4 * 170, 1),
                                   southermost_latitude = -80,
                                   singularity_longitude = 90,
                                   z = (0, 1),
                                   halo = (4, 4, 4))

h = regrid_bathymetry(grid, height_above_water=1, minimum_depth=5)

λ, φ, z = nodes(h)

land = interior(h) .> 0
interior(h)[land] .= NaN

fig = Figure(resolution=(2400, 1200))
ax = Axis(fig[1, 1])
heatmap!(ax, interior(h, :, :, 1), nan_color=:white, colorrange=(-5000, 0))

display(fig)
