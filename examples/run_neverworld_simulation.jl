using Oceananigans
using ClimaOcean.LimitedAreaSimulations: neverworld_simulation

simulation = neverworld_simulation(CPU())

simulation.stop_iteration = 1
run!(simulation)

