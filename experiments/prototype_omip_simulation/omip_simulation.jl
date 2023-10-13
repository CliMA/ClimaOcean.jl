using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using NCDatasets
using GLMakie
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState

using Downloads: download

temperature_filename = "THETA.1440x720x50.19920102.nc"
salinity_filename = "SALT.1440x720x50.19920102.nc"

# Downloaded from https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N

temperature_url = "https://www.dropbox.com/scl/fi/01h96yo2fhnnvt2zkmu0d/" *
                  "THETA.1440x720x50.19920102.nc?rlkey=ycso2v09gc6v2qb5j0lff0tjs&dl=0"

salinity_url = "https://www.dropbox.com/scl/fi/t068we10j5skphd461zg8/" *
               "SALT.1440x720x50.19920102.nc?rlkey=r5each0ytdtzh5icedvzpe7bw&dl=0"

isfile(temperature_filename) || download(temperature_url, temperature_filename)
isfile(salinity_filename)    || download(salinity_url,    salinity_filename)

temperature_ds = Dataset(temperature_filename)
salinity_ds = Dataset(salinity_filename)

# Construct vertical coordinate
depth = temperature_ds["DEPTH_T"][:]
zc = -reverse(depth)

# Interface depths from cell center depths
zf = (zc[1:end-1] .+ zc[2:end]) ./ 2
push!(zf, 0)

Δz = zc[2] - zc[1]
pushfirst!(zf, zf[1] - Δz)

T = temperature_ds["THETA"][:, :, :, 1]
S = salinity_ds["SALT"][:, :, :, 1]

T = convert(Array{Float32, 3}, T)
S = convert(Array{Float32, 3}, S)

T = reverse(T, dims=3)
S = reverse(S, dims=3)

missing_value = Float32(-9.9e22)

# Construct bottom_height depth by analyzing T
Nx, Ny′, Nz = size(T)
bottom_height = ones(Nx, Ny′) .* (zf[1] - Δz)

for i = 1:Nx, j = 1:Ny′
    @inbounds for k = Nz:-1:1
        if T[i, j, k] < -10
            bottom_height[i, j] = zf[k+1]
            break
        end
    end
end

# Grid construction
arch = CPU()
southern_limit = -80
northern_limit = 80
j₁ = 4 * (90 + southern_limit)
j₂ = 720 - 4 * (90 - northern_limit) + 1
Ny = j₂ - j₁ + 1

T = T[:, j₁:j₂, :]
S = S[:, j₁:j₂, :]
bottom_height = bottom_height[:, j₁:j₂]

grid = LatitudeLongitudeGrid(arch,
                             size = (Nx, Ny, Nz),
                             longitude = (0, 360),
                             latitude = (southern_limit, northern_limit),
                             z = zf,
                             topology = (Periodic, Bounded, Bounded))

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

using Oceananigans.Grids: φnodes, λnodes

λ = λnodes(grid, Center())
φ = φnodes(grid, Center())

fig = Figure()
ax = Axis(fig[1, 1])
heatmap!(ax, λ, φ, bottom_height)

# Model construction
teos10 = TEOS10EquationOfState()
buoyancy = SeawaterBuoyancy(equation_of_state=teos10)

model = HydrostaticFreeSurfaceModel(; grid, buoyancy,
                                    tracers = (:T, :S, :e),
                                    coriolis = HydrostaticSphericalCoriolis(),
                                    closure = CATKEVerticalDiffusivity())

