# # Diurnal large eddy simulation
#
# This script assembles a large-eddy simulation that resolves the marine boundary layer
# under a clear diurnal cycle, transitioning from nighttime cooling to daytime heating.
# The setup is based on typical summer solstice conditions near San Diego (33°N latitude).
#
# The atmosphere is prescribed with time-varying temperature, humidity, wind, and
# radiation fields that follow a realistic diurnal pattern.

using ClimaOcean
using Oceananigans
using Oceananigans.Units
using Printf
using CairoMakie

# ## Parameters
#
# Physical and numerical parameters for the simulation.
# Solar radiation parameters:
# The solar constant (top-of-atmosphere irradiance) is ~1361 W/m².
# Clear-sky atmospheric transmittance is about 75%, so surface irradiance is ~1020 W/m².
# At 33°N on summer solstice, the solar elevation at noon is:
#
# ``θ = 90° - |ϕ - δ| = 90° - |33° - 23.5°| ≈ 80.5°``
#
# where ``ϕ`` is latitude and ``δ`` is the solar declination (23.5° at summer solstice).
# Thus peak surface irradiance is approximately ``1361 × sin(80.5°) × 0.75 ≈ 1006`` W/m².

day = 86400   # seconds in a day
latitude = 33   # degrees (San Diego)
solar_constant = 1361           # W/m², top-of-atmosphere
clear_sky_transmittance = 0.75  # atmospheric transmittance on clear day
surface_solar_irradiance = solar_constant * clear_sky_transmittance  # ≈ 1020 W/m²

# ## Diurnal cycle functions
#
# These functions describe the diurnal evolution of atmospheric state variables
# based on typical clear-sky summer solstice conditions for coastal Southern California.
#
# The simulation starts at midnight (t = 0). Key features:
# - Sunrise around 6 AM, sunset around 8 PM (14-hour day)
# - Maximum solar radiation ~1000 W/m² at solar noon (see solar parameters above)
# - Maximum temperature around 2-3 PM (thermal lag)
# - Minimum temperature around 6 AM (just before sunrise)

"""
    solar_elevation(t)

Compute the solar elevation angle (radians) at time `t` for a given `latitude`.
Simplified model assuming summer solstice conditions.
"""
function solar_elevation(t)
    hour_angle = 2π * (t / day - 0.5)  # 0 at noon
    δ = deg2rad(23.5)  # summer solstice declination
    ϕ = deg2rad(latitude)
    sin_elevation = sin(ϕ) * sin(δ) + cos(ϕ) * cos(δ) * cos(hour_angle)
    return asin(clamp(sin_elevation, -1, 1))
end

"""
    diurnal_shortwave(t)

Clear-sky downwelling shortwave radiation (W/m²) at time `t`.
"""
function diurnal_shortwave(x, y, t)
    elevation = solar_elevation(t)
    daytime_shortwave = surface_solar_irradiance * sin(elevation)
    return max(0, daytime_shortwave)
end

"""
    diurnal_longwave(t, day)

Downwelling longwave radiation (W/m²) at time `t`.
Varies with atmospheric temperature: ~300 W/m² at night, ~380 W/m² during day.
"""
function diurnal_longwave(x, y, t)
    phase = 2π * t / day - π/2 - 2 * 2π / 24
    return 340 + 40 * sin(phase)
end

"""
    diurnal_temperature(t, day)

Air temperature (Kelvin) at time `t`.
Ranges from ~17°C (290 K) before sunrise to ~24°C (297 K) in mid-afternoon.
"""
function diurnal_temperature(x, y, t)
    T_min, T_max = 290, 297
    T_mean = (T_min + T_max) / 2
    T_amp = (T_max - T_min) / 2
    phase = 2π * (t / day - 6 / 24) - π/2
    return T_mean + T_amp * sin(phase)
end

"""
    diurnal_humidity(t, day)

Specific humidity (kg/kg) at time `t`.
Higher at night (~0.012), lower during day (~0.008).
"""
function diurnal_humidity(x, y, t)
    q_mean, q_amp = 0.010, 0.002
    phase = 2π * (t / day - 6 / 24) - π/2
    return q_mean - q_amp * sin(phase)
end

"""
    diurnal_wind_u(t, day)

Eastward wind component (m/s) at time `t`.
Land-sea breeze pattern: offshore at night, onshore during day.
"""
function diurnal_wind_u(x, y, t)
    phase = 2π * t / day - π
    return 2 + 3 * sin(phase)
end

diurnal_wind_v(x, y, t) = -1  # slight southward component

# ## Ocean LES configuration
#
# 3D grid with 4m spacing: 512m × 4m × 256m (Nx × Ny × Nz = 128 × 1 × 64).
# The thin y-dimension makes this essentially a 2D slice for computational efficiency.

arch = GPU()
Lx, Ly, Lz = 512, 512, 256  # meters
Nx, Ny, Nz = 128, 128, 64   # 4m grid spacing

grid = RectilinearGrid(arch; size = (Nx, Ny, Nz), halo = (5, 5, 5),
                       topology = (Periodic, Periodic, Bounded),
                       x = (0, Lx), y = (0, Ly), z = (-Lz, 0))

coriolis = FPlane(; latitude)
ocean = nonhydrostatic_ocean_simulation(grid; coriolis)
conjure_time_step_wizard!(ocean, cfl=0.7)

# Initial conditions: warm mixed layer with small perturbations

Tᵢ(x, y, z) = 22 + 0.01 * randn()  # °C, typical summer SST
Sᵢ(x, y, z) = 33 + 0.001 * randn()  # psu, typical coastal salinity
set!(ocean.model, T=Tᵢ, S=Sᵢ)

# ## Prescribed diurnal atmosphere
#
# Set up the atmosphere with time-varying fields over a full diurnal cycle.
# We use hourly time snapshots for smooth interpolation.

atmosphere_grid = RectilinearGrid(arch; size = (Nx, Ny),
                                  topology = (Periodic, Periodic, Flat),
                                  x = (0, Lx), y = (0, Ly))

atmosphere_times = range(0, 24hours, length=25)
atmosphere = PrescribedAtmosphere(atmosphere_grid, collect(atmosphere_times))

# Use set! to fill the atmosphere with diurnal cycle functions.
# We use closures to capture the parameters.

set!(atmosphere;
     u = diurnal_wind_u,
     v = diurnal_wind_v,
     T = diurnal_temperature,
     q = diurnal_humidity,
     shortwave = diurnal_shortwave,
     longwave = diurnal_longwave)

# ## Coupled ocean–atmosphere model
#
# The ocean simulation has a `TimeStepWizard` callback (added by `nonhydrostatic_ocean_simulation`)
# that adaptively adjusts `ocean.Δt` based on CFL constraints. We use the
# `synchronize_coupled_time_step!` callback to synchronize the coupled simulation's time step
# with the ocean's adaptive time step before each iteration.

coupled_model = OceanSeaIceModel(ocean; atmosphere)
simulation = Simulation(coupled_model; Δt=10.0, stop_time=12hours)

# Add callback to synchronize the coupled time step with the ocean's adaptive time step.
add_callback!(simulation, ClimaOcean.OceanSeaIceModels.align_component_steps!) 

# ## Output: snapshots and time-averages
#
# We save:
# 1. Snapshots of velocity and temperature fields for visualization
# 2. Time-averaged profiles of u and T (horizontally averaged)
# 3. Time-series of surface fluxes (both pointwise and horizontally-averaged)

u, v, w = ocean.model.velocities
T = ocean.model.tracers.T
Q = coupled_model.interfaces.net_fluxes.ocean.T
τx = coupled_model.interfaces.net_fluxes.ocean.u
τy = coupled_model.interfaces.net_fluxes.ocean.v

# Snapshot output every 10 minutes

Nz = grid.Nz
uˢ = view(u, :, :, Nz)
vˢ = view(v, :, :, Nz)
Tˢ = view(T, :, :, Nz)
u_avg = Average(u, dims=(1, 2))
T_avg = Average(T, dims=(1, 2))
Q_avg = Average(Q, dims=(1, 2))
τx_avg = Average(τx, dims=(1, 2))

outputs = (; uˢ, vˢ, Tˢ, Q, τx, τy, u_avg, T_avg, Q_avg, τx_avg)

simulation.output_writers[:snapshots] = JLD2Writer(ocean.model, outputs;
                                                   filename = "diurnal_les",
                                                   schedule = TimeInterval(10minutes),
                                                   overwrite_existing = true)

# Progress callback

function progress(sim)
    t = time(sim)
    hour = t / 3600
    ocean_Δt = sim.model.ocean.Δt
    u, v, w = sim.model.ocean.model.velocities
    T = sim.model.ocean.model.tracers.T
    Tmax, Tmin = maximum(interior(T)), minimum(interior(T))
    umax, wmax = maximum(abs, u), maximum(abs, w)
    Q = sim.model.interfaces.net_fluxes.ocean_surface.Q
    Qmax, Qmin = maximum(Q), minimum(Q)
    @info @sprintf("Hour %5.2f | Δt: %s (ocean: %s) | SST: %.2f–%.2f°C | Q: %.0f–%.0f W/m² | max|u|: %.2e, max|w|: %.2e",
                   hour, prettytime(sim.Δt), prettytime(ocean_Δt), Tmin, Tmax, Qmin, Qmax, umax, wmax)
end

add_callback!(simulation, progress, IterationInterval(100))

# ## Run the simulation

run!(simulation)

#=
# ## Animation
#
# Create a multi-panel animation showing:
# - Main panels: x-z slices of u and T
# - Top panels: time-series of momentum and heat flux
# - Side panels: time-averaged profiles of u and T

using JLD2

# Load data

u_snapshots = FieldTimeSeries("diurnal_les.jld2", "uˢ")
T_snapshots = FieldTimeSeries("diurnal_les.jld2", "Tˢ")
Q_avg_ts = FieldTimeSeries("diurnal_les.jld2", "Q_avg")
τx_avg_ts = FieldTimeSeries("diurnal_les.jld2", "τx_avg")
u_profiles = FieldTimeSeries("diurnal_les.jld2", "u_avg")
T_profiles = FieldTimeSeries("diurnal_les.jld2", "T_avg")

times = u_snapshots.times
Nt = length(times)

# Compute time-averaged profiles

u_mean = zeros(Nz)
T_mean = zeros(Nz)
for n in 1:Nt
    u_mean .+= interior(u_profiles[n], 1, 1, :)
    T_mean .+= interior(T_profiles[n], 1, 1, :)
end
u_mean ./= Nt
T_mean ./= Nt

# Extract coordinates

xc = xnodes(u_snapshots)
zc = znodes(u_snapshots)

# Build time-series arrays for flux plots

t_hours = times ./ hours
Q_avg_series = [interior(Q_avg_ts[n], 1, 1, 1) for n in 1:Nt]
τx_avg_series = [interior(τx_avg_ts[n], 1, 1, 1) for n in 1:Nt]

# Create figure

fig = Figure(size = (1400, 900), fontsize = 14)

# Layout:
# Row 1: flux time-series (spans columns 2-3)
# Row 2: u profile | u heatmap | T heatmap | T profile

ax_flux = Axis(fig[1, 2:3], xlabel = "Hour", ylabel = "Flux",
               title = "Surface fluxes")

ax_u_profile = Axis(fig[2, 1], xlabel = "⟨u⟩ (m/s)", ylabel = "z (m)",
                    title = "Time-avg u")

ax_u = Axis(fig[2, 2], xlabel = "x (m)", ylabel = "z (m)",
            title = "Horizontal velocity u")

ax_T = Axis(fig[2, 3], xlabel = "x (m)", ylabel = "z (m)",
            title = "Temperature T")

ax_T_profile = Axis(fig[2, 4], xlabel = "⟨T⟩ (°C)", ylabel = "z (m)",
                    title = "Time-avg T")

# Time-averaged profile plots (static)

lines!(ax_u_profile, u_mean, zc, color = :blue, linewidth = 2)
lines!(ax_T_profile, T_mean, zc, color = :red, linewidth = 2)

# Flux time-series (will update with vertical line marker)

lines!(ax_flux, times ./ hours, Q_avg_series, label = "Heat flux Q (W/m²)", color = :red)
lines!(ax_flux, times ./ hours, τx_avg_series .* 1000, label = "Mom. flux τx (×10³ N/m²)", color = :blue)
axislegend(ax_flux, position = :rt)

time_marker = Observable(0)
vlines!(ax_flux, time_marker, color = :black, linestyle = :dash, linewidth = 2)

# Heatmap observables

n = Observable(1)

u_slice = @lift interior(u_snapshots[$n], :, 1, :)
T_slice = @lift interior(T_snapshots[$n], :, 1, :)

# Determine color limits from data

u_lim = maximum(abs, interior(u_snapshots[Nt]))
T_lim_min, T_lim_max = minimum(interior(T_snapshots[1])), maximum(interior(T_snapshots[Nt]))

hm_u = heatmap!(ax_u, xc, zc, u_slice, colormap = :balance, colorrange = (-u_lim, u_lim))
hm_T = heatmap!(ax_T, xc, zc, T_slice, colormap = :thermal, colorrange = (T_lim_min - 0.5, T_lim_max + 0.5))

Colorbar(fig[3, 2], hm_u, vertical = false, label = "u (m/s)")
Colorbar(fig[3, 3], hm_T, vertical = false, label = "T (°C)")

# Title with time

title = @lift "Diurnal LES: Hour " * @sprintf("%.1f", times[$n] / 3600)
Label(fig[0, :], title, fontsize = 20, font = :bold)

# Record animation

record(fig, "diurnal_les_animation.mp4", 1:Nt; framerate = 12) do i
    n[] = i
    time_marker[] = times[i] / 3600
end

# ![Diurnal LES Animation](diurnal_les_animation.mp4)
=#