using Oceananigans
#using GLMakie
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using NCDatasets

# Load grid
dir = "/storage2/alir/bsose_i122"
grid_filename = "grid.nc"
grid_filepath = joinpath(dir, grid_filename)
grid_ds = Dataset(grid_filepath)

z = reverse(grid_ds["RF"][:])
zb = - grid_ds["Depth"][:, :]
φ = grid_ds["YG"][1, :]
Δφ = φ[end] - φ[end-1]
push!(φ, φ[end] + Δφ)

arch = GPU()
Nx = 6 * 360 # 1/6th degree
Ny = length(φ) - 1
Nz = length(z) - 1

grid = LatitudeLongitudeGrid(arch; z,
                             latitude = φ,
                             longitude = (0, 360),
                             size = (Nx, Ny, Nz),
                             halo = (5, 5, 5),
                             topology = (Periodic, Bounded, Bounded))

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(zb))

equation_of_state = TEOS10EquationOfState()

model = HydrostaticFreeSurfaceModel(; grid,
                                    tracers = (:T, :S, :e),
                                    buoyancy = SeawaterBuoyancy(; equation_of_state),
                                    coriolis = HydrostaticSphericalCoriolis(),
                                    free_surface = SplitExplicitFreeSurface(; grid, cfl=0.2),
                                    momentum_advection = VectorInvariant(),
                                    tracer_advection = WENO(),
                                    closure = CATKEVerticalDiffusivity())

# Initial condition
dir = "/storage2/alir/bsose_i122"

T_filename = "bsose_i122_2013to2017_1day_Theta.nc"
S_filename = "bsose_i122_2013to2017_1day_Salt.nc"
T_filepath = joinpath(dir, T_filename)
S_filepath = joinpath(dir, S_filename)
T_ds = Dataset(T_filepath)
S_ds = Dataset(S_filepath)
Tᵢ = T_ds["THETA"][:, :, :, 1]
Sᵢ = S_ds["SALT"][:, :, :, 1]

set!(model, T=Tᵢ, S=Sᵢ)

