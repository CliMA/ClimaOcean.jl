using ClimaOcean
using Oceananigans
using Oceananigans.Units
using DataDeps
using JLD2
using GLMakie

using Oceananigans.ImmersedBoundaries: mask_immersed_field!

architecture = CPU()
initial_conditions_path = datadep"near_global_one_degree/near_global_initial_conditions_360_150_48.jld2"
bathymetry_path = datadep"near_global_one_degree/near_global_bathymetry_360_150.jld2"

file = jldopen(initial_conditions_path)
grid = file["grid"]
T_data = file["T"]
S_data = file["S"]
close(file)

file = jldopen(bathymetry_path)
bathymetry = file["bathymetry"]
close(file)

immersed_grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bathymetry))

T = CenterField(immersed_grid)
S = CenterField(immersed_grid)

interior(T) .= T_data
interior(S) .= S_data

mask_immersed_field!(T, NaN)
mask_immersed_field!(S, NaN)

# Continue data below ocean bottom
ClimaOcean.DataWrangling.continue_downards!(T)
ClimaOcean.DataWrangling.continue_downards!(S)

# Inpaint data horizontally into continents
@time ClimaOcean.DataWrangling.inpaint_horizontally!(T)
@time ClimaOcean.DataWrangling.inpaint_horizontally!(S)

#####
##### Regrid
#####

fig = Figure(resolution=(2400, 1800))

i = 600
j = 300

axT = Axis(fig[1, 1], title="Temperature at k=Nz")
heatmap!(axT, interior(T, :, :, Nz))

axS = Axis(fig[1, 2], title="Salinity at k=Nz")
heatmap!(axS, interior(S, :, :, Nz))

axT = Axis(fig[2, 1], title="Temperature at k=1")
heatmap!(axT, interior(T, :, :, 1))

axS = Axis(fig[2, 2], title="Salinity at k=1")
heatmap!(axS, interior(S, :, :, 1))

axT = Axis(fig[3, 1], title="Temperature at i=$i")
heatmap!(axT, interior(T, i, :, :))

axS = Axis(fig[3, 2], title="Salinity at i=$i")
heatmap!(axS, interior(S, i, :, :))

axT = Axis(fig[4, 1], title="Temperature at j=$j")
heatmap!(axT, interior(T, :, j, :))

axS = Axis(fig[4, 2], title="Salinityat j=$j")
heatmap!(axS, interior(S, :, j, :))

display(fig)

