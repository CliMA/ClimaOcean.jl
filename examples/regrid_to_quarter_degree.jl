using ClimaOcean
using Oceananigans
using Oceananigans.Units
using Oceananigans.BuoyancyModels: ∂z_b
using DataDeps
using JLD2
using GLMakie

using ClimaOcean.DataWrangling: diffuse_tracers!
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using Oceananigans.ImmersedBoundaries: mask_immersed_field!

architecture = CPU()
initial_conditions_path = datadep"near_global_one_degree/near_global_initial_conditions_360_150_48.jld2"
bathymetry_path = datadep"near_global_one_degree/near_global_bathymetry_360_150.jld2"

file = jldopen(initial_conditions_path)
grid = file["grid"]
T_one_coarse_data = file["T"]
S_one_coarse_data = file["S"]
close(file)

T_one_coarse = CenterField(grid)
S_one_coarse= CenterField(grid)
interior(T_one_coarse) .= T_one_coarse_data
interior(S_one_coarse) .= S_one_coarse_data

#####
##### Regrid to high resolution in vertical
#####

# New z-grid generation, 5m spacing at surface
# Fixed spacing in the upper ocean
Δz₀ = 5.0         # surface layer grid spacing
h₀ = 100          # surface layer extent

# Generate surface layer grid
z = [-Δz₀ * (k-1) for k = 1:ceil(h₀ / Δz₀)]

# Generate stretched interior grid
γ = 1.02 # stretching parameter
Lz₀ = 5kilometers # minimum domain depth

while z[end] > - Lz₀
    Δz = (z[end-1] - z[end])^γ
    push!(z, round(z[end] - Δz, digits=1))
end

# Reverse grid to be right-side-up
z_refined = reverse(z)

# Infer domain parameters
Lz = z_refined[1]
Nz = length(z_refined) - 1

vertically_refined_one_degree_grid = LatitudeLongitudeGrid(architecture;
                                                           size = (360, 150, Nz),
                                                           longitude = (-180, 180),
                                                           latitude = (-75, 75),
                                                           halo = (5, 5, 5),
                                                           z = z_refined)

T_one = CenterField(vertically_refined_one_degree_grid)
S_one = CenterField(vertically_refined_one_degree_grid)

regrid!(T_one, T_one_coarse)
regrid!(S_one, S_one_coarse)

#####
##### Regrid one degree initial condition to quarter degree grid with high vertical resolution
#####

# Intermediate grid: quarter degree in x, one degree in y 
quarter_degree_x_grid = LatitudeLongitudeGrid(architecture;
                                              size = (1440, 150, Nz),
                                              longitude = (-180, 180),
                                              latitude = (-75, 75),
                                              halo = (5, 5, 5),
                                              z = z_refined)

T_quarter_x = CenterField(quarter_degree_x_grid)
S_quarter_x = CenterField(quarter_degree_x_grid)

regrid!(T_quarter_x, T_one)
regrid!(S_quarter_x, S_one)

# Quarter degree grid
quarter_degree_grid = LatitudeLongitudeGrid(architecture;
                                            size = (1440, 600, Nz),
                                            longitude = (-180, 180),
                                            latitude = (-75, 75),
                                            halo = (5, 5, 5),
                                            z = z_refined)

T_quarter = CenterField(quarter_degree_grid)
S_quarter = CenterField(quarter_degree_grid)

regrid!(T_quarter, T_quarter_x)
regrid!(S_quarter, S_quarter_x)

bathymetry_path = datadep"near_global_quarter_degree/near_global_bathymetry_1440_600.jld2"
file = jldopen(bathymetry_path)
quarter_degree_bathymetry = file["bathymetry"]
close(file)

immersed_quarter_degree_grid = ImmersedBoundaryGrid(quarter_degree_grid,
                                                    GridFittedBottom(quarter_degree_bathymetry))

T = CenterField(immersed_quarter_degree_grid, data=T_quarter.data)
S = CenterField(immersed_quarter_degree_grid, data=S_quarter.data)

step(x, d, c) = 1/2 * (tanh((x - c) / d) + 1)
vertical_scale(x, y, z, t) = 10 + 190 * step(abs(z), 200, 1000)

diffuse_tracers!(immersed_quarter_degree_grid;
                 tracers = (; T, S),
                 horizontal_scale = 0,
                 vertical_scale,
                 fractional_time_step = 0.1)

teos10 = TEOS10EquationOfState(reference_density=1020)
buoyancy = SeawaterBuoyancy(equation_of_state=teos10)
N²_op = KernelFunctionOperation{Center, Center, Face}(∂z_b, immersed_quarter_degree_grid;  
                                                      computed_dependencies=(buoyancy, (; T, S)))
N² = Field(N²_op)
compute!(N²)

#####
##### Visualize
#####

fig = Figure(resolution=(1800, 2400))

mask_immersed_field!(T, NaN)
mask_immersed_field!(S, NaN)
mask_immersed_field!(N², 0)

i = 600
j = 300
x, y, z = nodes(N²)
Nlims = (1e-6, 5e-4)

axT = Axis(fig[1, 1], title="Temperature at k=Nz")
heatmap!(axT, interior(T, :, :, Nz))

axS = Axis(fig[1, 2], title="Salinity at k=Nz")
heatmap!(axS, interior(S, :, :, Nz))

axN = Axis(fig[1, 3], title="Buoyancy gradient at k=Nz")
heatmap!(axN, x, y, interior(N², :, :, Nz), colorrange=Nlims)

axT = Axis(fig[2, 1], title="Temperature at k=1")
heatmap!(axT, interior(T, :, :, 1))

axS = Axis(fig[2, 2], title="Salinity at k=1")
heatmap!(axS, interior(S, :, :, 1))

axN = Axis(fig[2, 3], title="Buoyancy gradient at k=Nz-10")
heatmap!(axN, x, y, interior(N², :, :, Nz-10), colorrange=Nlims)

axT = Axis(fig[3, 1], title="Temperature at i=$i")
heatmap!(axT, interior(T, i, :, :))

axS = Axis(fig[3, 2], title="Salinity at i=$i")
heatmap!(axS, interior(S, i, :, :))

axN = Axis(fig[3, 3], title="Buoyancy gradient at i=$i")
heatmap!(axN, y, z, interior(N², i, :, :), colorrange=Nlims)

axT = Axis(fig[4, 1], title="Temperature at j=$j")
heatmap!(axT, interior(T, :, j, :))

axS = Axis(fig[4, 2], title="Salinity at j=$j")
heatmap!(axS, interior(S, :, j, :))

axN = Axis(fig[4, 3], title="Buoyancy gradient at j=$j")
heatmap!(axN, x, z, interior(N², :, j, :), colorrange=Nlims)

display(fig)

