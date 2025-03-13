#=
using ClimaOcean
using Oceananigans
using Oceananigans.Units
using CFTime
using Dates
using Printf

arch = CPU()

Nx = 144
Ny = 60
Nz = 40

grid = LatitudeLongitudeGrid(arch;
                             size = (Nx, Ny, Nz),
                             halo = (7, 7, 7),
                             z = (-6000, 0),
                             latitude  = (-75, 75),
                             longitude = (0, 360))

ocean = ocean_simulation(grid)

# date = DateTimeProlepticGregorian(1993, 1, 1)
# set!(ocean.model, T=ECCOMetadata(:temperature; dates=date),
#                   S=ECCOMetadata(:salinity; dates=date))

radiation = Radiation(arch)
=#
atmosphere = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(41))

coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)

simulation = Simulation(coupled_model; Δt=10, stop_iteration=50)

wall_time = Ref(time_ns())

function progress(sim)
    ocean = sim.model.ocean
    u, v, w = ocean.model.velocities
    T = ocean.model.tracers.T

    Tmax = maximum(interior(T))
    Tmin = minimum(interior(T))

    umax = (maximum(abs, interior(u)),
            maximum(abs, interior(v)),
            maximum(abs, interior(w)))

    step_time = 1e-9 * (time_ns() - wall_time[])

    msg = @sprintf("Iter: %d, simulation time: %s, atmosphere time: %s, Δt: %s", iteration(sim), prettytime(sim), prettytime(atmosphere.clock.time), prettytime(sim.Δt))
    msg *= @sprintf(", max|u|: (%.2e, %.2e, %.2e) m s⁻¹, extrema(T): (%.2f, %.2f) ᵒC, wall time: %s",
                    umax..., Tmax, Tmin, prettytime(step_time))

    @info msg

    wall_time[] = time_ns()
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(1))

outputs = merge(ocean.model.tracers, ocean.model.velocities)

simulation.output_writers[:surface] = JLD2OutputWriter(ocean.model, outputs;
                                                       schedule = IterationInterval(2),
                                                       filename = "checkpointer_mwe_surface",
                                                       indices = (:, :, grid.Nz),
                                                       with_halos = true,
                                                       overwrite_existing = true,
                                                       array_type = Array{Float32})

output_dir = "."
prefix = "checkpointer_mwe"

ocean.output_writers[:checkpoint] = Checkpointer(ocean.model;
                                                      schedule = IterationInterval(2),
                                                      prefix = prefix,
                                                      dir = output_dir,
                                                      verbose = true,
                                                      overwrite_existing = true)

# @show simulation

using ClimaOcean.OceanSeaIceModels: OceanSeaIceModel, OSIMSIM
import Oceananigans.OutputWriters: write_output!
import Oceananigans.Simulations: pickup!

write_output!(c::Checkpointer, model::OceanSeaIceModel) = write_output!(c::Checkpointer, model.ocean.model)

function pickup!(sim::OSIMSIM, pickup)
    @info "i try to properly pick up"
    set!(sim.model.ocean, pickup)
    return nothing
end

# run!(simulation)

# @info "simulation run for 10 iterations; you should have a checkpointer at 8"

# simulation.stop_iteration = 20

run!(simulation, pickup=true)
