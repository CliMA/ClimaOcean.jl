using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using JLD2
using Printf

filename = "antarctic_circumpolar_grid_initial_conditions.jld2"
file = jldopen(filename)

Tᵢ = Array{Float64, 3}(file["T"])
Sᵢ = Array{Float64, 3}(file["S"])

longitude = file["longitude"]
latitude = file["latitude"]
z = file["z"]
zb = file["bottom_height"]

close(file)

arch = GPU()
Nx = 6 * 360
Ny = length(latitude) - 1
Nz = length(z) - 1

grid = LatitudeLongitudeGrid(arch; longitude, latitude, z,
                             size = (Nx, Ny, Nz),
                             halo = (7, 7, 7),
                             topology = (Periodic, Bounded, Bounded))

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(zb))

equation_of_state = TEOS10EquationOfState()

model = HydrostaticFreeSurfaceModel(; grid,
                                    tracers = (:T, :S, :e),
                                    buoyancy = SeawaterBuoyancy(; equation_of_state),
                                    coriolis = HydrostaticSphericalCoriolis(),
                                    free_surface = SplitExplicitFreeSurface(; grid, cfl=0.2),
                                    momentum_advection = VectorInvariant(),
                                    tracer_advection = WENO(),
                                    closure = CATKEVerticalDiffusivity())

set!(model, T=Tᵢ, S=Sᵢ)

simulation = Simulation(model, Δt=1minute, stop_iteration=200) #stop_time=1day)

function progress(sim)
    msg1 = string("Iter: ", iteration(sim), ", time: ", prettytime(sim))

    u, v, w = sim.model.velocities
    msg2 = @sprintf(", max|u|: (%.2e, %.2e, %.2e) m s⁻¹",
                    maximum(abs, u),
                    maximum(abs, v),
                    maximum(abs, w))

    @info msg1 * msg2

    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

output_dir = "." #/nobackup1/glwagner/"
outputs = merge(model.velocities, model.tracers)
Nz = size(grid, 3)
simulation.output_writers[:jld2] = JLD2OutputWriter(model, outputs,
                                                    schedule = IterationInterval(10),
                                                    indices = (:, :, Nz),
                                                    dir = output_dir,
                                                    filename = "antarctic_circumpolar_surface.jld2",
                                                    overwrite_existing = true)

run!(simulation)

