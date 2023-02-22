using ClimaOcean
using Oceananigans
using GLMakie
using JLD2
using DataDeps

one_degree_grid = LatitudeLongitudeGrid(size = (360, 150, 1),
                                        longitude = (0, 360),
                                        latitude = (-75, 75),
                                        z = (0, 1),
                                        topology = (Periodic, Bounded, Bounded))

bathymetry_path = "near_global_bathymetry_360_150.jld2" #datadep"near_global_one_degree/bathymetry_lat_lon_360_150.jld2"
file = jldopen(bathymetry_path)
bathymetry_data = file["bathymetry"]
close(file)

bathymetry = Field{Center, Center, Nothing}(one_degree_grid)

land = bathymetry_data .> 0
bathymetry_data[land] .= 0
bathymetry .= bathymetry_data

one_two_degree_grid = LatitudeLongitudeGrid(size = (180, 150, 1),
                                            longitude = (0, 360),
                                            latitude = (-75, 75),
                                            z = (0, 1),
                                            topology = (Periodic, Bounded, Bounded))

bathymetry_12 = Field{Center, Center, Nothing}(one_two_degree_grid)
regrid!(bathymetry_12, bathymetry) 

two_degree_grid = LatitudeLongitudeGrid(size = (180, 75, 1),
                                        longitude = (0, 360),
                                        latitude = (-75, 75),
                                        z = (0, 1),
                                        topology = (Periodic, Bounded, Bounded))

bathymetry_two_degree = Field{Center, Center, Nothing}(two_degree_grid)
regrid!(bathymetry_two_degree, bathymetry_12) 

fig = Figure()
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 2])
heatmap!(ax1, interior(bathymetry, :, :, 1))
heatmap!(ax2, interior(bathymetry_two_degree, :, :, 1))
display(fig)
