using Oceananigans
using Oceananigans.Units
using Oceananigans.Coriolis: ActiveCellEnstrophyConservingScheme
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using JLD2
using Printf

filename = "antarctic_circumpolar_grid_initial_conditions.jld2"
file = jldopen(filename)
Tᵢ = Array{Float64, 3}(file["T"])
Sᵢ = Array{Float64, 3}(file["S"])
Tᵢ = reverse(Tᵢ, dims=3)
Sᵢ = reverse(Sᵢ, dims=3)
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

@show grid

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(zb))

equation_of_state = TEOS10EquationOfState()

model = HydrostaticFreeSurfaceModel(; grid,
                                    tracers = (:T, :S, :e),
                                    buoyancy = SeawaterBuoyancy(; equation_of_state),
                                    coriolis = HydrostaticSphericalCoriolis(scheme = ActiveCellEnstrophyConservingScheme()),
                                    free_surface = SplitExplicitFreeSurface(; grid, cfl=0.2),
                                    momentum_advection = VectorInvariant(vorticity_scheme = WENO(),
                                                                         divergence_scheme = WENO(),
                                                                         vertical_scheme = WENO()),
                                    tracer_advection = WENO(),
                                    closure = CATKEVerticalDiffusivity())

@show model

set!(model, T=Tᵢ, S=Sᵢ)

simulation = Simulation(model, Δt=5minutes, stop_time=30days)

start_time = Ref(time_ns())

function progress(sim)
    elapsed = 1e-9 * (time_ns() - start_time[])

    msg1 = string("Iter: ", iteration(sim),
                  ", time: ", prettytime(sim),
                  ", wall time: ", prettytime(elapsed))

    u, v, w = sim.model.velocities
    msg2 = @sprintf(", max|u|: (%.2e, %.2e, %.2e) m s⁻¹",
                    maximum(abs, u),
                    maximum(abs, v),
                    maximum(abs, w))
    @info msg1 * msg2

    start_time[] = time_ns()

    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

output_dir = "." #/nobackup1/glwagner/"
outputs = merge(model.velocities, model.tracers)
Nz = size(grid, 3)
simulation.output_writers[:jld2] = JLD2OutputWriter(model, outputs,
                                                    schedule = TimeInterval(1day),
                                                    indices = (:, :, Nz),
                                                    dir = output_dir,
                                                    filename = "antarctic_circumpolar_surface.jld2",
                                                    overwrite_existing = true)

run!(simulation)

