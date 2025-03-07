using ClimaOcean
using Oceananigans
using Oceananigans.Units
using CFTime
using Dates
using Printf
using Glob

arch = CPU()

Nx = 144
Ny = 60
Nz = 40

depth = 6000meters
z_faces = exponential_z_faces(; Nz, depth)

grid = LatitudeLongitudeGrid(arch;
                             size = (Nx, Ny, Nz),
                             halo = (7, 7, 7),
                             z = z_faces,
                             latitude  = (-75, 75),
                             longitude = (0, 360))

ocean = ocean_simulation(grid)

# date = DateTimeProlepticGregorian(1993, 1, 1)
# set!(ocean.model, T=ECCOMetadata(:temperature; dates=date),
#                   S=ECCOMetadata(:salinity; dates=date))

radiation = Radiation(arch)

atmosphere = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(41))

coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)

simulation = Simulation(coupled_model; Δt=10, stop_iteration=10)

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

    msg = @sprintf("Iter: %d, time: %s, Δt: %s", iteration(sim), prettytime(sim), prettytime(sim.Δt))
    msg *= @sprintf(", max|u|: (%.2e, %.2e, %.2e) m s⁻¹, extrema(T): (%.2f, %.2f) ᵒC, wall time: %s",
                    umax..., Tmax, Tmin, prettytime(step_time))

    @info msg 

    wall_time[] = time_ns()
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(1))

outputs = merge(ocean.model.tracers, ocean.model.velocities)

ocean.output_writers[:surface] = JLD2OutputWriter(ocean.model, outputs;
                                                  schedule = IterationInterval(2),
                                                  filename = "checkpointer_mwe_surface",
                                                  indices = (:, :, grid.Nz),
                                                  with_halos = true,
                                                  overwrite_existing = true,
                                                  array_type = Array{Float32})

output_dir = "."
prefix = "checkpointer_mwe"

ocean.output_writers[:checkpoint] = Checkpointer(ocean.model;
                                                 schedule = IterationInterval(4),
                                                 prefix = prefix,
                                                #  cleanup = true,
                                                 dir = output_dir,
                                                 verbose = true,
                                                 overwrite_existing = true)

coupled_checkpointer = Checkpointer(coupled_model;
                                                 schedule = IterationInterval(4),
                                                 prefix = prefix,
                                                #  cleanup = true,
                                                 dir = output_dir,
                                                 verbose = true,
                                                 overwrite_existing = true)

#=

@show simulation

run!(simulation)

@info "simulation run for 50 iterations; you should have a checkpointer at 40"

checkpoint_file = prefix * "_iteration40.jld2"

set!(simulation, checkpoint_file)
    
coupled_model = OceanSeaIceModel(simulation.model.ocean; atmosphere, radiation)

simulation = Simulation(coupled_model; Δt=10, stop_iteration=20)
simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

run!(simulation)
=#