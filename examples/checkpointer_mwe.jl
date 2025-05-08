using ClimaOcean
using Oceananigans
using CFTime
using Dates
using Printf

arch = CPU()

Nx = 120
Ny = 50
Nz = 20

grid = LatitudeLongitudeGrid(arch;
                             size = (Nx, Ny, Nz),
                             halo = (7, 7, 7),
                             z = (-6000, 0),
                             latitude  = (-75, 75),
                             longitude = (0, 360))


function testbed_coupled_simulation(grid; stop_iteration=8)
    ocean = ocean_simulation(grid)
        
    radiation = Radiation(arch)
    
    atmosphere = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(4))
    
    coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)
    
    return Simulation(coupled_model; Δt=10, stop_iteration)
end

simulation = testbed_coupled_simulation(grid; stop_iteration=8)

wall_time = Ref(time_ns())

function progress(sim)
    ocean = sim.model.ocean
    atmosphere = sim.model.atmosphere

    step_time = 1e-9 * (time_ns() - wall_time[])

    msg = @sprintf("iteration: %d, sim time: %s, atmos time: %s, ocean time: %s, wall time: %s",
                   iteration(sim), sim.model.clock.time, atmosphere.clock.time, ocean.model.clock.time, prettytime(step_time))

    @info msg

    wall_time[] = time_ns()
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(1))

simulation.output_writers[:checkpoint] = Checkpointer(ocean.model;
                                                      schedule = IterationInterval(3),
                                                      prefix = "checkpointer",
                                                      dir = ".",
                                                      verbose = true,
                                                      overwrite_existing = true)

run!(simulation)

new_simulation = testbed_coupled_simulation(grid; stop_iteration=14)

new_simulation.output_writers[:checkpoint] = Checkpointer(ocean.model;
                                                      schedule = IterationInterval(3),
                                                      prefix = "checkpointer",
                                                      dir = ".",
                                                      verbose = true,
                                                      overwrite_existing = true)


run!(new_simulation, pickup=true)
