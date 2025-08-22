using ClimaOcean
using PythonCall
using Oceananigans
using Printf

#####
##### A Prognostic Python Ocean (Veros) Simulation
#####

# We import the Veros 4 degree ocean simulation setup, which consists of a near-global ocean
# with a uniform resolution of 4 degrees in both latitude and longitude and a latitude range spanning
# from 80S to 80N. The setup is defined in the `veros.setups.global_4deg` module.

# Before importing the setup, we need to ensure that the Veros module is loaded
# and that every output is removed to avoid conflicts.

VerosModule = Base.get_extension(ClimaOcean, :ClimaOceanPythonCallExt)
VerosModule.remove_outputs(:global_4deg)

# Actually loading and instantiating the Veros setup in the variable `ocean`.
# This setup uses by default a different time-step for tracers and momentum, 
# so we set it to the same value (1800 seconds) for both.

ocean = VerosModule.VerosOceanSimulation("global_4deg", :GlobalFourDegreeSetup)

VerosModule.veros_settings_set!(ocean, "dt_tracer", 1800.0)
VerosModule.veros_settings_set!(ocean, "dt_mom",    1800.0)

##### 
##### A Prescribed Atmosphere (JRA55)
##### 

spectral_grid = SpectralGrid(trunc=63, nlayers=8, Grid=FullClenshawGrid)

humidity_flux_ocean = PrescribedOceanHumidityFlux(spectral_grid)
humidity_flux_land = SurfaceLandHumidityFlux(spectral_grid)
surface_humidity_flux = SurfaceHumidityFlux(ocean=humidity_flux_ocean, land=humidity_flux_land)

ocean_heat_flux = PrescribedOceanHeatFlux(spectral_grid)
land_heat_flux = SurfaceLandHeatFlux(spectral_grid)
surface_heat_flux = SurfaceHeatFlux(ocean=ocean_heat_flux, land=land_heat_flux)

atmosphere_model = PrimitiveWetModel(spectral_grid;
                                     surface_heat_flux,
                                     surface_humidity_flux,
                                     ocean = nothing,
                                     sea_ice = nothing) # This is provided by ClimaSeaIce

atmosphere = initialize!(atmosphere_model)
initialize!(atmosphere)

function initialize_atmospheric_state!(simulation::SpeedyWeather.Simulation)
    progn, diagn, model  = SpeedyWeather.unpack(simulation)

    (; time) = progn.clock                           # current time

    # set the tendencies back to zero for accumulation
    fill!(diagn.tendencies, 0, typeof(model))

    if model.physics                   
        # calculate all parameterizations
        SpeedyWeather.parameterization_tendencies!(diagn, progn, time, model)
    end
    
    return nothing
end

initialize_atmospheric_state!(atmosphere)
# atmosphere.model.feedback.verbose = false

#####
##### An ice-free ocean forced by a prescribed atmosphere
#####

radiation = Radiation(ocean_emissivity=0, sea_ice_emissivity=0)
coupled_model = OceanSeaIceModel(ocean, nothing; atmosphere, radiation)
simulation = Simulation(coupled_model; Δt = 1800, stop_iteration = 100000)

#####
##### A simple progress callback
#####

# We set up a progress callback that will print the current time, iteration, and maximum velocities
# at every 5 iterations. It also collects the surface velocity fields and the net fluxes
# into the arrays `s`, `tx`, and `ty` for later visualization.

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

#####
##### Let's go!!
#####

run!(simulation)
