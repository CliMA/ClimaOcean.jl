using Oceananigans
using ClimaOcean.LimitedAreaSimulations: neverworld_simulation
using GLMakie

#=
simulation = neverworld_simulation(GPU(), horizontal_size=(240, 280))

model = simulation.model
grid = model.grid

@show grid

Nx, Ny, Nz = size(grid)
j = round(Int, Ny/2)

#=
fig = Figure()
ax = Axis(fig[1, 1])
h = Field{Center, Center, Nothing}(grid)
parent(h) .= parent(grid.immersed_boundary.bottom_height)
#lines!(ax, Array(view(h, 1:Nx, j)))
heatmap!(ax, Array(interior(h, 1:Nx, 1:Ny, 1)))
display(fig)
=#

simulation.stop_iteration = 100
run!(simulation)
=#

fig = Figure()
ax = Axis(fig[1, 1])

b = model.tracers.b
bⱼ = Array(interior(b, :, j, :))
heatmap!(ax, bⱼ, colorrange=(0.05, 0.06))

