using Oceananigans
using Oceananigans.Units
using ClimaOcean.LimitedAreaSimulations: neverworld_simulation
using GLMakie
using Printf

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

function progress(sim) 
    b = sim.model.tracers.b
    e = sim.model.tracers.e
    u, v, w = sim.model.velocities

    msg = @sprintf("Iter: %d, time: %s, extrema(b): (%6.2e, %6.2e)",
                   iteration(sim), prettytime(sim), minimum(b), maximum(b))

    msg *= @sprintf(", extrema(e): (%6.2e, %6.2e)", minimum(e), maximum(e))

    msg *= @sprintf(", max|u|: %6.2e, max|w|: %6.2e",
                    maximum(maximum(abs, q) for q in (u, v, w)), maximum(abs, w))

    @info msg

    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

simulation.stop_time = 7days
run!(simulation)

fig = Figure()
axbxy = Axis(fig[1, 1])
axzxy = Axis(fig[1, 2])
axb = Axis(fig[2, 1])
axe = Axis(fig[2, 2])

u, v, w = model.velocities
using Oceananigans.Operators: ζ₃ᶠᶠᶜ
ζop = KernelFunctionOperation{Face, Face, Center}(ζ₃ᶠᶠᶜ, grid, u, v)


ζ = compute!(Field(ζop))
b = model.tracers.b
e = model.tracers.e
bⱼ = Array(interior(b, :, j, :))
eⱼ = Array(interior(e, :, j, :))
bₖ = Array(interior(b, :, :, Nz))
ζₖ = Array(interior(ζ, :, :, Nz))

heatmap!(axbxy, bₖ, colorrange=(0.05, 0.06))
heatmap!(axzxy, ζₖ)

heatmap!(axb, bⱼ, colorrange=(0.05, 0.06))
heatmap!(axe, eⱼ)

