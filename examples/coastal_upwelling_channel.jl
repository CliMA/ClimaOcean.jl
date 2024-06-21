using Oceananigans
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Oceananigans.Units
using GLMakie

arch = CPU()
Nx = 64
Ny = 64
Nz = 4

N² = 1e-5

Lx = 100kilometers
Ly = 100kilometers
Lz = 50

grid = RectilinearGrid(arch, 
                       size = (Nx, Ny, Nz),
                       halo = (7, 7, 7),
                       x = (0, Lx),
                       y = (0, Ly),
                       z = (-Lz, 0),
                       topology = (Periodic, Bounded, Bounded))

#=
# We want the underwater slope to provide a depth of 500 m at y = 0 and the full 4 km by y =200. It follows
# a hyperbolic tangent curve with its center at y ~= 150 at a depth of ~ -4000 + (4000 - 500)/2 = -2250 m
y₀ = 150kilometers
Δy = 25kilometers

slope_depth = 500
basin_depth = 4000

""" Varies smoothly between 0 when y << y₀ to 1 when y >> y₀ + Δy, over a region of width Δy. """
step(y, y₀, Δy) = (1 + tanh((y - y₀) / Δy)) / 2

underwater_slope(x, y) = (-basin_depth - slope_depth + (slope_depth - basin_depth) * tanh((y - y₀) / Δy)) / 2
# TODO: add 50km troughs to the slope, dependent on x and y. Use appendix C to add detail here

ocean_grid = ImmersedBoundaryGrid(underlying_ocean_grid, GridFittedBottom(underwater_slope))
=#

# bottom_height(x, y) = -Lz * (x > 30kilometers)
# grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

model = HydrostaticFreeSurfaceModel(; grid,
                                    momentum_advection = WENOVectorInvariant(vorticity_order=5),
                                    tracer_advection = WENO(order=5),
                                    closure = CATKEVerticalDiffusivity(),
                                    tracers = (:b, :e),
                                    buoyancy = BuoyancyTracer())

bi(x, y, z) = N² * z
ui(x, y, z) = 1e-1 * (2rand() - 1)
set!(model, u=ui, v=ui, b=bi)

simulation = Simulation(model, Δt=1e-2, stop_iteration=10)

run!(simulation)

b = model.tracers.b
e = model.tracers.e
u, v, w = model.velocities

Nz = size(grid, 3)

bn = interior(b, :, :, Nz)
en = interior(e, :, :, Nz)
un = interior(u, :, :, Nz)

fig = Figure()

axb = Axis(fig[1, 1])
axe = Axis(fig[1, 2])
axu = Axis(fig[1, 3])

heatmap!(axb, bn)
heatmap!(axe, en)
heatmap!(axu, un)
heatmap!(axb, bn)
