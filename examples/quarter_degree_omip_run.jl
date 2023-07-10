using GLMakie
using Oceananigans
using ClimaOcean.Bathymetry: regrid_bathymetry
using ClimaOcean.VerticalGrids: stretched_vertical_faces, PowerLawStretching

#####
##### Quarter degree near-global grid
#####

arch = CPU()

z_faces = stretched_vertical_faces(surface_layer_Δz = 8,
                                   surface_layer_height = 128,
                                   stretching = PowerLawStretching(1.02),
                                   minimum_depth = 6000)

grid = LatitudeLongitudeGrid(arch;
                             size = (4 * 360, 4 * 160, 1),
                             latitude = (-80, 80),
                             longitude = (-180, 180),
                             z = z_faces,
                             halo = (4, 4, 4))

h = regrid_bathymetry(grid, height_above_water=1, minimum_depth=5)

λ, φ, z = nodes(h)

land = interior(h) .> 0
interior(h)[land] .= NaN

fig = Figure(resolution=(2400, 1200))
ax = Axis(fig[1, 1], xlabel="Longitude", ylabel="Latitude")
heatmap!(ax, λ, φ, interior(h, :, :, 1), nan_color=:white, colorrange=(-5000, 0))

display(fig)

