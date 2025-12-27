# # A Python Ocean Simulation at 4ᵒ Resolution Forced by JRA55 Reanalysis and Initialized from ECCO
#
# This example showcases the use of ClimaOcean's PythonCall extension to run a
# near-global ocean simulation at 4-degree resolution using the Veros ocean model.
# The ocean is forced by the JRA55 reanalysis data and initialized from the ECCO
# state estimate.
#
# For this example, we need Oceananigans, ClimaOcean, Dates, CUDA, and
# CairoMakie to visualize the simulation.

using ClimaOcean
using PythonCall
using Oceananigans, Oceananigans.Units
using CairoMakie
using Printf

# We import the Veros 4 degree ocean simulation setup, which consists of a near-global ocean
# with a uniform resolution of 4 degrees in both latitude and longitude and a latitude range spanning
# from 80S to 80N. The setup is defined in the `veros.setups.global_4deg` module.

# Before importing the setup, we need to ensure that the Veros module is installed and loaded
# and that every output is removed to avoid conflicts.

VerosModule = Base.get_extension(ClimaOcean, :ClimaOceanVerosExt)

VerosModule.install_veros()
VerosModule.remove_outputs(:global_4deg)

# Actually loading and instantiating the Veros setup in the variable `ocean`.
# This setup uses by default a different time-step for tracers and momentum, 
# so we set it to the same value (1800 seconds) for both.

ocean = VerosModule.VerosOceanSimulation("global_4deg", :GlobalFourDegreeSetup)

# The loaded Veros setup contains a `set_forcing` method which computes the fluxes as restoring from climatology.
# We need to disable the hardcoded forcing function from the loaded Veros so that our forcings are not overwritten

ocean.setup.set_forcing = @pyeval("lambda state: None")

# The loaded setup had different time-steps for tracer and momentum. ClimaOcean handles only oceans with the
# same time steps, so we need to align the tracer variables' timestep with the momentum.

set!(ocean, "dt_tracer", 1800.0; path=:settings)
set!(ocean, "dt_mom",    1800.0; path=:settings)

# We force the 4-degree setup with a prescribed atmosphere based on the JRA-55 reanalysis data.
# This includes 2-meter wind velocity, temperature, humidity, downwelling longwave and shortwave
# radiation, as well as freshwater fluxes.

atmos = JRA55PrescribedAtmosphere(; backend = JRA55NetCDFBackend(10))

# The coupled ocean--atmosphere model.
# We use the default radiation model and we do not couple an ice model for simplicity.

radiation = Radiation()
coupled_model = OceanSeaIceModel(ocean, nothing; atmosphere=atmos, radiation)
simulation = Simulation(coupled_model; Δt = 1800, stop_time = 60days)

# We set up a progress callback that will print the current time, iteration, and maximum velocities
# every 10days. We also set up another callback that collects the surface velocity fields and the net fluxes
# into arrays for later visualization.

wall_time = Ref(time_ns())

us = coupled_model.interfaces.exchanger.ocean.state.u
vs = coupled_model.interfaces.exchanger.ocean.state.v

sptemp = Field(sqrt(us^2 + vs^2))

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

τx = []
τy = []
JS = []
JT = []
sp = []

function save_variables(sim)
    push!(sp, compute!(sptemp))
    push!(τx, deepcopy(sim.coupled_model.interfaces.net_fluxes.ocean.u))
    push!(τy, deepcopy(sim.coupled_model.interfaces.net_fluxes.ocean.v))
    push!(JS, deepcopy(sim.coupled_model.interfaces.net_fluxes.ocean.S))
    push!(JT, deepcopy(sim.coupled_model.interfaces.net_fluxes.ocean.T))
end

add_callback!(simulation, progress, TimeInterval(10days))
add_callback!(simulation, save_variables, IterationInterval(5))

# Let's run the simulation!

run!(simulation)

iter = Observable(1)
spi = @lift(sp[$iter])
τxi = @lift(τx[$iter])
τyi = @lift(τy[$iter])
JSi = @lift(JS[$iter])
JTi = @lift(JT[$iter])
Nt  = length(sp)

fig = Figure(resolution = (1000, 1500))
ax1 = Axis(fig[1, 1]; title = "Surface speed (m/s)", xlabel = "", ylabel = "Latitude")
ax2 = Axis(fig[2, 1]; title = "Zonal wind stress (N/m²)", xlabel = "", ylabel = "Latitude")
ax3 = Axis(fig[3, 1]; title = "Meridional wind stress (N/m²)", xlabel = "", ylabel = "Latitude")
ax4 = Axis(fig[4, 1]; title = "Temperature flux (ᵒC m/s)", xlabel = "", ylabel = "Latitude")
ax5 = Axis(fig[5, 1]; title = "Surface flux (psu m/s)", xlabel = "Longitude", ylabel = "Latitude")

grid = coupled_model.interfaces.exchanger.grid

λ = λnodes(grid, Center())
φ = φnodes(grid, Center())

hm1 = heatmap!(ax1, λ, φ, spi, colormap = :ice, colorrange = (0, 0.15))
hm2 = heatmap!(ax2, λ, φ, τxi, colormap = :bwr, colorrange = (-0.2, 0.2))
hm3 = heatmap!(ax3, λ, φ, τyi, colormap = :bwr, colorrange = (-0.2, 0.2))
hm4 = heatmap!(ax4, λ, φ, JTi, colormap = :thermal)
hm5 = heatmap!(ax5, λ, φ, JSi, colormap = :haline)

Colorbar(fig[1, 2], hm1)
Colorbar(fig[2, 2], hm2)
Colorbar(fig[3, 2], hm3)
Colorbar(fig[4, 2], hm4)
Colorbar(fig[5, 2], hm5)

CairoMakie.record(fig, "veros_ocean_surface.mp4", 1:Nt, framerate = 8) do nn
    iter[] = nn
end
nothing #hide

# ![](veros_ocean_surface.mp4)
