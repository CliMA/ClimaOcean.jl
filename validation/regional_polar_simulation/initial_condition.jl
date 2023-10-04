using NCDatasets

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

# Tᵢ = T_ds["THETA"][:, :, :, 1]
# Sᵢ = S_ds["SALT"][:, :, :, 1]
