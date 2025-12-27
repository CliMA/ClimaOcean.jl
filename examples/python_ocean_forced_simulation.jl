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
# We replace it with a custom function that only computes the TKE forcing (which depends on the wind stresses
# that we set in ClimaOcean). This way our u, v, T, S forcings are not overwritten. 
# The `set_forcing_tke_only` method defined below is modified from the `set_forcing` method defined in 
# https://github.com/team-ocean/veros/blob/main/veros/setups/global_4deg/global_4deg.py

pyexec("""
def set_forcing_tke_only(state):
    from veros.core.operators import numpy as npx, update, at
    from veros import KernelOutput

    vs = state.variables
    settings = state.settings
    
    if settings.enable_tke:
        vs.forc_tke_surface = update(
            vs.forc_tke_surface,
            at[1:-1, 1:-1],
            npx.sqrt(
                (0.5 * (vs.surface_taux[1:-1, 1:-1] + vs.surface_taux[:-2, 1:-1]) / settings.rho_0) ** 2
                + (0.5 * (vs.surface_tauy[1:-1, 1:-1] + vs.surface_tauy[1:-1, :-2]) / settings.rho_0) ** 2
            ) ** 1.5,
        )

    return KernelOutput(
        surface_taux=vs.surface_taux,
        surface_tauy=vs.surface_tauy,
        forc_tke_surface=vs.forc_tke_surface,
        forc_temp_surface=vs.forc_temp_surface,
        forc_salt_surface=vs.forc_salt_surface,
    )

ocean.set_forcing = set_forcing_tke_only
""", Main, (ocean=ocean.setup,))

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
# every 10days. We also set up another callback that collects the surface prognostic variables
# into arrays for later visualization.

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

u = []
v = []
S = []
T = []

function save_variables(sim)
    push!(u, deepcopy(sim.model.interfaces.exchanger.ocean.state.u))
    push!(v, deepcopy(sim.model.interfaces.exchanger.ocean.state.v))
    push!(S, deepcopy(sim.model.interfaces.exchanger.ocean.state.S))
    push!(T, deepcopy(sim.model.interfaces.exchanger.ocean.state.T))
end

add_callback!(simulation, progress, TimeInterval(10days))
add_callback!(simulation, save_variables, IterationInterval(10))

# Let's run the simulation!

run!(simulation)

iter = Observable(1)
ui = @lift(u[$iter])
vi = @lift(v[$iter])
Si = @lift(S[$iter])
Ti = @lift(T[$iter])
Nt  = length(sp)

fig = Figure(resolution = (1000, 1500))
ax1 = Axis(fig[1, 1]; title = "Surface zonal velocity (m/s)", xlabel = "", ylabel = "Latitude")
ax2 = Axis(fig[2, 1]; title = "Surface meridional velocity (m/s)", xlabel = "", ylabel = "Latitude")
ax3 = Axis(fig[3, 1]; title = "Surface temperature (N/m²)", xlabel = "", ylabel = "Latitude")
ax4 = Axis(fig[4, 1]; title = "Surface salinity (psu)", xlabel = "", ylabel = "Latitude")

grid = coupled_model.interfaces.exchanger.grid

λ = λnodes(grid, Center())
φ = φnodes(grid, Center())

hm2 = heatmap!(ax1, λ, φ, ui, colormap = :bwr,     colorrange = (-0.2, 0.2))
hm3 = heatmap!(ax2, λ, φ, vi, colormap = :bwr,     colorrange = (-0.2, 0.2))
hm4 = heatmap!(ax3, λ, φ, Ti, colormap = :thermal, colorrange = (-1, 30))
hm5 = heatmap!(ax4, λ, φ, Si, colormap = :haline,  colorrange = (32, 37))

Colorbar(fig[1, 2], hm1)
Colorbar(fig[2, 2], hm2)
Colorbar(fig[3, 2], hm3)
Colorbar(fig[4, 2], hm4)

CairoMakie.record(fig, "veros_ocean_surface.mp4", 1:Nt, framerate = 8) do nn
    iter[] = nn
end
nothing #hide

# ![](veros_ocean_surface.mp4)
