using Oceananigans
using GLMakie

# Setup
grid = RectilinearGrid(size=(64, 64, 32), x=(-4, 4), y=(-4, 4), z=(0, 3), topology=(Periodic, Periodic, Bounded))
h = Field{Center, Center, Nothing}(grid)
mountain(x, y) = 1.5 * exp(-(x^2 + y^2) / 2) + 0.5 * rand()
set!(h, mountain)

# Fake model

Nx, Ny, Nz = size(h)
TX, TY, TZ = Oceananigans.Grids.topology(h.grid)
auxgrid = RectilinearGrid(size=(Nx, Ny, 1), x=(-4, 4), y=(-4, 4), z=(0, 3), topology=(TX, TY, TZ))

# κ ∼ dx^2 / dt
Δx = minimum_xspacing(auxgrid)
Δy = minimum_yspacing(auxgrid)
Δ² = 1 / (1/Δx^2 + 1/Δy^2)
κ = Δ² # Δt = 1

closure = HorizontalScalarDiffusivity(; κ)
velocities = PrescribedVelocityFields()

model = HydrostaticFreeSurfaceModel(grid=auxgrid; closure, velocities, tracers=:c)
set!(model, c=h)

Δt = 0.1
Oceananigans.TimeSteppers.first_time_step!(model, Δt)
for n = 2:100
    Oceananigans.TimeSteppers.time_step!(model, Δt)
end

fig = Figure()
axh = Axis(fig[1, 1])
axc = Axis(fig[2, 1])

c = model.tracers.c
∇h = sqrt(∂x(h)^2 + ∂y(h)^2)
∇c = sqrt(∂x(c)^2 + ∂y(c)^2)
@show maximum(∇h)
@show maximum(∇c)

heatmap!(axh, h)
heatmap!(axc, c)
display(fig)
