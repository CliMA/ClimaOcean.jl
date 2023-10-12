using GLMakie
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

land = interior(h) .> 0
interior(h)[land] .= NaN

fig = Figure(resolution=(2400, 1200))
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

grid = LatitudeLongitudeGrid(CPU();
                             size = (25 * 50, 55 * 50, 1),
                             latitude = (25, 50),
                             longitude = (-10, 45),
                             z = (0, 1),
                             halo = (4, 4, 4))

h = regrid_bathymetry(grid, height_above_water=1)

λ, φ, z = nodes(h)

land = interior(h) .> 0
interior(h)[land] .= NaN

fig = Figure(resolution=(2400, 1200))
ax = Axis(fig[1, 1])
heatmap!(ax, λ, φ, interior(h, :, :, 1), nan_color=:white) #, colorrange=(-5000, 0))

display(fig)

