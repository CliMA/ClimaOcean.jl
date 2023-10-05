# Utility script for building an initial condition from BSOSE data on tartarus

using Oceananigans
using NCDatasets
using JLD2

# Initial condition
dir = "/storage2/alir/bsose_i122"
T_filename = "bsose_i122_2013to2017_1day_Theta.nc"
S_filename = "bsose_i122_2013to2017_1day_Salt.nc"
grid_filename = "grid.nc"

T_filepath = joinpath(dir, T_filename)
S_filepath = joinpath(dir, S_filename)
grid_filepath = joinpath(dir, grid_filename)

T_ds = Dataset(T_filepath)
S_ds = Dataset(S_filepath)
grid_ds = Dataset(grid_filepath)

z = reverse(grid_ds["RF"][:])
zb = - grid_ds["Depth"][:, :]

longitude = (0, 360)
φ = grid_ds["YG"][1, :]
Δφ = φ[end] - φ[end-1]
push!(φ, φ[end] + Δφ)

Nx = 6 * 360
Ny = length(φ) - 1
Nz = length(z) - 1

#=
grid = LatitudeLongitudeGrid(CPU(); z,
                             latitude = φ,
                             longitude = (0, 360),
                             size = (Nx, Ny, Nz),
                             halo = (7, 7, 7),
                             topology = (Periodic, Bounded, Bounded))

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(zb))

file["grid"] = grid
=#

time = T_ds["time"][:]
T = T_ds["THETA"][:, :, :, 1]
S = S_ds["SALT"][:, :, :, 1]

ic_filename = "sose_grid_initial_conditions.jld2"
file = jldopen(ic_filename, "a+")

file["T"] = T
file["S"] = S
file["time"] = time

file["longitude"] = longitude
file["latitude"] = φ
file["z"] = z
file["bottom_height"] = zb

close(file)
