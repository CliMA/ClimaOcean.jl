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

start_time = time_ns()
# T⁺ = T⁻ + κ * ∂z(∂z(T⁻))
for i = 1:10
    global start_time
    T_one .= T_one .+ 0.1 * ∂z(∂z(T_one))
    S_one .= S_one .+ 0.1 * ∂z(∂z(S_one))
    elapsed = 1e-9 * (time_ns() - start_time)
    elapsed_str = prettytime(elapsed)
    @info "One degree diffusion step $i, $elapsed_str"
    start_time = time_ns()
end

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

#####
##### Regrid
#####

fig = Figure(resolution=(1800, 2400))

mask_immersed_field!(T, NaN)
mask_immersed_field!(S, NaN)

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

