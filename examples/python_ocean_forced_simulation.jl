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
# at every 5 iterations. It also collects the surface velocity fields and the net fluxes
# into the arrays `s`, `tx`, and `ty` for later visualization.

wall_time = Ref(time_ns())

sp = []
S  = []
T  = []
tx = []
ty = []
Js = []
Jt = []

us = coupled_model.interfaces.exchanger.ocean.state.u
vs = coupled_model.interfaces.exchanger.ocean.state.v

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
    push!(sp, deepcopy(interior(stmp, :, :, 1)))
    push!(S,  deepcopy(interior(coupled_model.interfaces.exchanger.ocean.state.S, :, :, 1)))
    push!(T,  deepcopy(interior(coupled_model.interfaces.exchanger.ocean.state.T, :, :, 1)))
    push!(tx, deepcopy(interior(coupled_model.interfaces.net_fluxes.ocean.u, :, :, 1)))
    push!(ty, deepcopy(interior(coupled_model.interfaces.net_fluxes.ocean.v, :, :, 1)))
    push!(Js, deepcopy(interior(coupled_model.interfaces.net_fluxes.ocean.S, :, :, 1)))
    push!(Jt, deepcopy(interior(coupled_model.interfaces.net_fluxes.ocean.T, :, :, 1)))

    return nothing
end

add_callback!(simulation, progress, IterationInterval(5))

# Let's run the simulation!

run!(simulation)

# We can now visualize the surface speed and wind stress at the ocean surface
# over the course of the simulation.

iter = Observable(1)
spi  = @lift(sp[$iter])
Jsi  = @lift(tx[$iter])
Jti  = @lift(ty[$iter])
Si   = @lift(S[$iter])
Ti   = @lift(T[$iter])
Nt   = length(sp)

fig = Figure(resolution = (1200, 600))
ax1 = Axis(fig[1, 1]; title = "Surface speed (m/s)", xlabel = "Longitude", ylabel = "Latitude")
ax2 = Axis(fig[1, 2]; title = "Zonal wind stress (N/m²)", xlabel = "Longitude")
ax3 = Axis(fig[1, 3]; title = "Meridional wind stress (N/m²)", xlabel = "Longitude")
ax4 = Axis(fig[2, 1]; title = "Surface temperature (ᵒC)", xlabel = "Longitude", ylabel = "Latitude")
ax5 = Axis(fig[2, 2]; title = "Surface salinity (psu)", xlabel = "Longitude")
ax6 = Axis(fig[2, 3]; title = "Surface salinity flux", xlabel = "Longitude")

grid = coupled_model.interfaces.exchanger.grid

λ = λnodes(grid, Center())
φ = φnodes(grid, Center())

heatmap!(ax1, λ, φ, spi, colormap = :ice, colorrange = (0, 0.15))
heatmap!(ax2, λ, φ, txi, colormap = :bwr, colorrange = (-0.2, 0.2))
heatmap!(ax3, λ, φ, tyi, colormap = :bwr, colorrange = (-0.2, 0.2))

heatmap!(ax4, λ, φ, Ti,  colormap = :thermal, colorrange = (-1, 30))
heatmap!(ax5, λ, φ, Si,  colormap = :haline,  colorrange = (32, 37))
heatmap!(ax6, λ, φ, Jsi, colormap = :bwr)

CairoMakie.record(fig, "veros_ocean_surface.mp4", 1:Nt, framerate = 8) do nn
    iter[] = nn
end
nothing #hide

# ![](veros_ocean_surface.mp4)
