using Oceananigans
using ClimaOcean.LimitedAreaSimulations: neverworld_simulation
using GLMakie

simulation = neverworld_simulation(CPU())

model = simulation.model
grid = model.grid
Nx, Ny, Nz = size(grid)
h = grid.immersed_boundary.bottom_height

fig = Figure()
ax = Axis(fig[1, 1])
heatmap!(ax, h[1:Nx, 1:Ny])
display(fig)

simulation.stop_iteration = 1
run!(simulation)

