using GLMakie
using Oceananigans
using ClimaOcean.Bathymetry: regrid_bathymetry

grid = LatitudeLongitudeGrid(CPU();
                             size = (4 * 360, 4 * 160, 1),
                             latitude = (-80, 80),
                             longitude = (-180, 180),
                             z = (0, 1),
                             halo = (4, 4, 4))

one_degree_h = regrid_bathymetry(grid, height_above_water=1)

λ, φ, z = nodes(one_degree_h)

land = interior(one_degree_h) .> 0
interior(one_degree_h)[land] .= NaN

fig = Figure(resolution=(2400, 1200))
ax = Axis(fig[1, 1])
heatmap!(ax, λ, φ, interior(one_degree_h, :, :, 1), nan_color=:white, colorrange=(-5000, 0))

λp = -112.45
φp = 42.86
text = "😻"
text!(ax, λp, φp; text)

display(fig)

