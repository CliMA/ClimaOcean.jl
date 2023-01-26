using ClimaOcean
using Oceananigans
using DataDeps
using JLD2
using GLMakie

using Oceananigans.ImmersedBoundaries: mask_immersed_field!

ic_path = datadep"near_global_one_degree/initial_conditions_month_01_360_150_48.jld2"
bathymetry_path = datadep"near_global_one_degree/bathymetry_lat_lon_360_150.jld2"

architecture = CPU()
z_faces = ClimaOcean.VerticalGrids.z_49_levels_10_to_400_meter_spacing

one_degree_grid = LatitudeLongitudeGrid(architecture;
                                        size = (360, 150, 48),
                                        longitude = (-180, 180),
                                        latitude = (-75, 75),
                                        halo = (5, 5, 5),
                                        z = z_faces)

file = jldopen(bathymetry_path)
bathymetry = file["bathymetry"]
close(file)

immersed_grid = ImmersedBoundaryGrid(one_degree_grid, GridFittedBottom(bathymetry))

file = jldopen(ic_path)
T_data = file["T"]
S_data = file["S"]
close(file)

T = CenterField(immersed_grid) #one_degree_grid)
S = CenterField(immersed_grid) #one_degree_grid)

T .= T_data
S .= S_data

mask_immersed_field!(T, NaN)
mask_immersed_field!(S, NaN)

ClimaOcean.DataWrangling.continue_downards!(T)
ClimaOcean.DataWrangling.continue_downards!(S)

fig = Figure(resolution=(2400, 1800))

Nz = one_degree_grid.Nz

axT = Axis(fig[1, 1], title="Temperature at k=Nz")
heatmap!(axT, interior(T, :, :, Nz))

axS = Axis(fig[1, 2], title="Salinity at k=Nz")
heatmap!(axS, interior(S, :, :, Nz))

axT = Axis(fig[2, 1], title="Temperature at k=1")
heatmap!(axT, interior(T, :, :, 1))

axS = Axis(fig[2, 2], title="Salinityat k=1")
heatmap!(axS, interior(S, :, :, 1))


