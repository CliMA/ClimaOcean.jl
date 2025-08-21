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
simulation = Simulation(coupled_model; Δt = 1800, stop_iteration = 100000)

wall_time = Ref(time_ns())

s  = []
tx = []
ty = []

us = coupled_model.interfaces.exchanger.exchange_ocean_state.u
vs = coupled_model.interfaces.exchanger.exchange_ocean_state.v

stmp = Field(sqrt(us^2 + vs^2))

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

    compute!(stmp)
    push!(s,  deepcopy(interior(stmp, :, :, 1)))
    push!(tx, deepcopy(interior(coupled_model.interfaces.net_fluxes.ocean_surface.u, :, :, 1) .* 1020))
    push!(ty, deepcopy(interior(coupled_model.interfaces.net_fluxes.ocean_surface.v, :, :, 1) .* 1020))

    return nothing
end

add_callback!(simulation, progress, IterationInterval(5))

run!(simulation)
