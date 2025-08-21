using ClimaOcean
using PythonCall
using Oceananigans
using Printf

VerosModule = Base.get_extension(ClimaOcean, :ClimaOceanPythonCallExt)
VerosModule.remove_outputs(:global_4deg)

ocean = VerosModule.veros_ocean_simulation("global_4deg", :GlobalFourDegreeSetup)
VerosModule.veros_settings_set!(ocean, "dt_tracer", 1800.0)

atmos = JRA55PrescribedAtmosphere(; backend = JRA55NetCDFBackend(10))
radiation = Radiation()
coupled_model = OceanSeaIceModel(ocean, nothing; atmosphere=atmos, radiation)
simulation = Simulation(coupled_model; Δt = 1800, stop_iteration = 100)

wall_time = Ref(time_ns())

function progress(sim)
    ocean   = sim.model.ocean
    umax = maximum(PyArray(ocean.setup.state.variables.u))
    vmax = maximum(PyArray(ocean.setup.state.variables.v))
    wmax = maximum(PyArray(ocean.setup.state.variables.w))

    step_time = 1e-9 * (time_ns() - wall_time[])

    msg1 = @sprintf("time: %s, iteration: %d, Δt: %s, ", prettytime(sim), iteration(sim), prettytime(sim.Δt))
    msg5 = @sprintf("maximum(u): (%.2f, %.2f, %.2f) m/s, ", umax, vmax, wmax)
    msg6 = @sprintf("wall time: %s \n", prettytime(step_time))

    @info msg1 * msg5 * msg6

     wall_time[] = time_ns()

     return nothing
end

add_callback!(simulation, progress, IterationInterval(10))

run!(simulation)
