using ClimaOcean.NearGlobalSimulations: aqua_planet_simulation
using Oceananigans
using Oceananigans.Units
using GLMakie

# closure = RiBasedVerticalDiffusivity() 
simulation = aqua_planet_simulation(CPU())

@info "Running a simulation with Δt = $(prettytime(simulation.Δt))"

simulation.stop_iteration = 100
run!(simulation)

@info "Simulation took $(prettytime(simulation.run_wall_time))."

model = simulation.model
grid = model.grid
u, v, w = model.velocities
T, S, e = model.tracers

fig = Figure()

ax = Axis(fig[1, 1])
heatmap!(ax, interior(u, :, :, grid.Nz))

ax = Axis(fig[1, 2])
heatmap!(ax, interior(v, :, :, grid.Nz))

ax = Axis(fig[2, 1])
heatmap!(ax, interior(T, :, :, grid.Nz))

ax = Axis(fig[2, 2])
heatmap!(ax, interior(e, :, :, grid.Nz))

display(fig)
