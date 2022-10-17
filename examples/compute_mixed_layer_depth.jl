using ClimaOcean.NearGlobalSimulations: one_degree_near_global_simulation
using ClimaOcean.Diagnostics
using Oceananigans
using Oceananigans

using GLMakie

# Build the simulation
simulation = one_degree_near_global_simulation(CPU())
model = simulation.model

h = MixedLayerDepthField(model)
compute!(h)
x, y, z = nodes(model.tracers.T) 

fig = Figure(resolution=(1800, 1200))
ax = Axis(fig[1, 1], xlabel="Longitude", ylabel="Latitude")
hm = heatmap!(ax, x, y, interior(h, :, :, 1), colorrange=(-1, 100))

Colorbar(fig[1, 2], hm;
         vertical = true,
         tellheight = false,
         flipaxis = true,
         label = "Mixed layer depth (m)")

display(fig)

