using Oceananigans
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using MAT

vars = matread("grid.mat")

z = reverse(vars["RF"][:])
h = - vars["Depth"]
φ = vars["YC"][1, :]
Nx, Ny = size(h)
Nz = length(z) - 1

Δφ = 1/6
φ₁ = minimum(φ) - Δφ / 2
φ₂ = maximum(φ) + Δφ / 2

longitude = (0, 360)
latitude = (φ₁, φ₂)

grid = LatitudeLongitudeGrid(; longitude, latitude, z,
                             size = (Nx, Ny, Nz),
                             halo = (5, 5, 5),
                             topology = (Periodic, Bounded, Bounded))

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(h))

equation_of_state = TEOS10EquationOfState()

model = HydrostaticFreeSurfaceModel(; grid,
                                    tracers = (:T, :S, :e),
                                    buoyancy = SeawaterBuoyancy(; equation_of_state),
                                    coriolis = HydrostaticSphericalCoriolis(),
                                    free_surface = SplitExplicitFreeSurface(; grid, cfl=0.2),
                                    momentum_advection = VectorInvariant(),
                                    tracer_advection = WENO(),
                                    closure = CATKEVerticalDiffusivity())

