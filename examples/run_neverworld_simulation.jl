using Oceananigans
using ClimaOcean.LimitedAreaSimulations: neverworld_simulation
using GLMakie

simulation = neverworld_simulation(CPU(), horizontal_size=(60, 70))

model = simulation.model
grid = model.grid

@show grid

Nx, Ny, Nz = size(grid)
j = round(Int, Ny/2)

fig = Figure()
ax = Axis(fig[1, 1])
h = grid.immersed_boundary.bottom_height
lines!(ax, h[1:Nx, j])
#heatmap!(ax, h[1:Nx, 1:Ny])
display(fig)

simulation.stop_iteration = 100
run!(simulation)

fig = Figure()
ax = Axis(fig[1, 1])

b = model.tracers.b
bⱼ = interior(b, :, j, :)
heatmap!(ax, bⱼ)

