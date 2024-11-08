# # Single column ocean simulation forced by JRA55 re-analysis
#
# In this example, we simulate the evolution of an ocean water column 
# forced by an atmosphere derived from the JRA55 re-analysis.
# The simulated column is located at ocean station
# Papa (144.9ᵒ W and 50.1ᵒ N)
#
# ## Install dependencies
#
# First let's make sure we have all required packages installed.

# ```julia
# using Pkg
# pkg"add Oceananigans, ClimaOcean, CairoMakie"
# ```

using ClimaOcean

using Oceananigans
using Oceananigans.Units
using Oceananigans.BuoyancyModels: buoyancy_frequency
using Oceananigans.Units: Time

using CairoMakie
using Printf

# # Construct the grid
#
# First, we construct a single column grid with 2 meter spacing
# located at ocean station Papa.

# Ocean station papa location
location_name = "ocean_station_papa"
λ★, φ★ = 35.1, 50.1

grid = RectilinearGrid(size = 200,
                       x = λ★,
                       y = φ★,
                       z = (-400, 0),
                       topology = (Flat, Flat, Bounded))

# # An "ocean simulation"
#
# Next, we use ClimaOcean's ocean_simulation constructor to build a realistic
# ocean simulation on the single column grid,

ocean = ocean_simulation(grid; Δt=10minutes, coriolis=FPlane(latitude = φ★))

# which wraps around the ocean model

ocean.model

# We set initial conditions from ECCO:

set!(ocean.model, T=ECCOMetadata(:temperature), S=ECCOMetadata(:salinity))

# # A prescribed atmosphere based on JRA55 re-analysis
#
# We build a PrescribedAtmosphere at the same location as the single colunm grid
# which is based on the JRA55 reanalysis.

simulation_days = 31
snapshots_per_day = 8 # corresponding to JRA55's 3-hour frequency
last_time = simulation_days * snapshots_per_day
atmosphere = JRA55_prescribed_atmosphere(1:last_time;
                                         longitude = λ★,
                                         latitude = φ★,
                                         #longitude = (λ★ - 1/4, λ★ + 1/4),
                                         #latitude  = (φ★ - 1/4, φ★ + 1/4),
                                         backend = InMemory(),
                                         include_rivers_and_icebergs = false)

# This builds a representation of the atmosphere on the small grid

atmosphere.grid

# Let's take a look at the atmospheric state

ua = interior(atmosphere.velocities.u, 1, 1, 1, :)
va = interior(atmosphere.velocities.v, 1, 1, 1, :)
Ta = interior(atmosphere.tracers.T, 1, 1, 1, :)
qa = interior(atmosphere.tracers.q, 1, 1, 1, :)
t_days = atmosphere.times / days

fig = Figure(size=(800, 600))
axu = Axis(fig[2, 1], xlabel="Days since Jan 1 1990", ylabel="Atmosphere \n velocity (m s⁻¹)")
axT = Axis(fig[3, 1], xlabel="Days since Jan 1 1990", ylabel="Atmosphere \n temperature (K)")
axq = Axis(fig[4, 1], xlabel="Days since Jan 1 1990", ylabel="Atmosphere \n specific humidity")
Label(fig[1, 1], "Atmospheric state over ocean station Papa", tellwidth=false)

lines!(axu, t_days, ua, label="Zonal velocity")
lines!(axu, t_days, va, label="Meridional velocity")
ylims!(axu, -6, 6)
axislegend(axu, framevisible=false, nbanks=2, position=:lb)

lines!(axT, t_days, Ta)
lines!(axq, t_days, qa)

display(fig)

radiation = Radiation()
coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)
simulation = Simulation(coupled_model, Δt=ocean.Δt, stop_time=30days)

wall_clock = Ref(time_ns())

function progress(sim)
    msg = "Ocean Station Papa"
    msg *= string(", iter: ", iteration(sim), ", time: ", prettytime(sim))

    elapsed = 1e-9 * (time_ns() - wall_clock[])
    msg *= string(", wall time: ", prettytime(elapsed))
    wall_clock[] = time_ns()

    u, v, w = sim.model.ocean.model.velocities
    msg *= @sprintf(", max|u|: (%.2e, %.2e)", maximum(abs, u), maximum(abs, v))

    T = sim.model.ocean.model.tracers.T
    S = sim.model.ocean.model.tracers.S
    e = sim.model.ocean.model.tracers.e

    τx = first(sim.model.fluxes.total.ocean.momentum.u)
    τy = first(sim.model.fluxes.total.ocean.momentum.v)
    Q = first(sim.model.fluxes.total.ocean.heat)

    u★ = sqrt(sqrt(τx^2 + τy^2))

    Nz = size(T, 3)
    msg *= @sprintf(", u★: %.2f m s⁻¹", u★)
    msg *= @sprintf(", Q: %.2f W m⁻²",  Q)
    msg *= @sprintf(", T₀: %.2f ᵒC", first(interior(T, 1, 1, Nz)))
    msg *= @sprintf(", extrema(T): (%.2f, %.2f) ᵒC", minimum(T), maximum(T))
    msg *= @sprintf(", S₀: %.2f g/kg", first(interior(S, 1, 1, Nz)))
    msg *= @sprintf(", e₀: %.2e m² s⁻²", first(interior(e, 1, 1, Nz)))

    @info msg
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

# Build flux outputs
τx = coupled_model.fluxes.total.ocean.momentum.u
τy = coupled_model.fluxes.total.ocean.momentum.v
JT = coupled_model.fluxes.total.ocean.tracers.T
Js = coupled_model.fluxes.total.ocean.tracers.S
E  = coupled_model.fluxes.turbulent.fields.water_vapor
Qc = coupled_model.fluxes.turbulent.fields.sensible_heat
Qv = coupled_model.fluxes.turbulent.fields.latent_heat
ρₒ = coupled_model.fluxes.ocean_reference_density
cₚ = coupled_model.fluxes.ocean_heat_capacity

Q = ρₒ * cₚ * JT
ρτx = ρₒ * τx
ρτy = ρₒ * τy
N² = buoyancy_frequency(ocean.model)
κc = ocean.model.diffusivity_fields.κc

fluxes = (; ρτx, ρτy, E, Js, Qv, Qc)
auxiliary_fields = (; N², κc)
fields = merge(ocean.model.velocities, ocean.model.tracers, auxiliary_fields)

# Slice fields at the surface
outputs = merge(fields, fluxes)

filename = "single_column_omip_$(location_name)"

simulation.output_writers[:jld2] = JLD2OutputWriter(ocean.model, outputs; filename,
                                                    schedule = TimeInterval(3hours),
                                                    overwrite_existing = true)

run!(simulation)

filename *= ".jld2"

u  = FieldTimeSeries(filename, "u")
v  = FieldTimeSeries(filename, "v")
T  = FieldTimeSeries(filename, "T")
S  = FieldTimeSeries(filename, "S")
e  = FieldTimeSeries(filename, "e")
N² = FieldTimeSeries(filename, "N²")
κ  = FieldTimeSeries(filename, "κc")

Qv = FieldTimeSeries(filename, "Qv")
Qc = FieldTimeSeries(filename, "Qc")
Js = FieldTimeSeries(filename, "Js")
Ev = FieldTimeSeries(filename, "E")
ρτx = FieldTimeSeries(filename, "ρτx")
ρτy = FieldTimeSeries(filename, "ρτy")

Nz = size(T, 3)
times = Qc.times

ua  = atmosphere.velocities.u
va  = atmosphere.velocities.v
Ta  = atmosphere.tracers.T
qa  = atmosphere.tracers.q
Qlw = atmosphere.downwelling_radiation.longwave
Qsw = atmosphere.downwelling_radiation.shortwave
Pr  = atmosphere.freshwater_flux.rain
Ps  = atmosphere.freshwater_flux.snow

Nt   = length(times)
uat  = zeros(Nt)
vat  = zeros(Nt)
Tat  = zeros(Nt)
qat  = zeros(Nt)
Qswt = zeros(Nt)
Qlwt = zeros(Nt)
Pt   = zeros(Nt)

for n = 1:Nt
    t = times[n]
    uat[n]  =  ua[1, 1, 1, Time(t)]
    vat[n]  =  va[1, 1, 1, Time(t)]
    Tat[n]  =  Ta[1, 1, 1, Time(t)]
    qat[n]  =  qa[1, 1, 1, Time(t)]
    Qswt[n] = Qsw[1, 1, 1, Time(t)]
    Qlwt[n] = Qlw[1, 1, 1, Time(t)]
    Pt[n]   =  Pr[1, 1, 1, Time(t)] + Ps[1, 1, 1, Time(t)]
end

set_theme!(Theme(linewidth=3))

fig = Figure(size=(2400, 1800))

axτ = Axis(fig[1, 1:2], xlabel="Days since Oct 1 1992", ylabel="Wind stress (N m⁻²)")
axu = Axis(fig[2, 1:2], xlabel="Days since Oct 1 1992", ylabel="Velocities (m s⁻¹)")
axQ = Axis(fig[1, 3:4], xlabel="Days since Oct 1 1992", ylabel="Heat flux (W m⁻²)")
axT = Axis(fig[2, 3:4], xlabel="Days since Oct 1 1992", ylabel="Surface temperature (ᵒC)")
axF = Axis(fig[1, 5:6], xlabel="Days since Oct 1 1992", ylabel="Freshwater volume flux (m s⁻¹)")
axS = Axis(fig[2, 5:6], xlabel="Days since Oct 1 1992", ylabel="Surface salinity (g kg⁻¹)")

axuz = Axis(fig[3, 1], xlabel="Velocities (m s⁻¹)",                ylabel="z (m)")
axTz = Axis(fig[3, 2], xlabel="Temperature (ᵒC)",                  ylabel="z (m)")
axSz = Axis(fig[3, 3], xlabel="Salinity (g kg⁻¹)",                 ylabel="z (m)")
axNz = Axis(fig[3, 4], xlabel="Buoyancy frequency (s⁻²)",          ylabel="z (m)")
axκz = Axis(fig[3, 5], xlabel="Eddy diffusivity (m² s⁻¹)",         ylabel="z (m)", xscale=log10)
axez = Axis(fig[3, 6], xlabel="Turbulent kinetic energy (m² s⁻²)", ylabel="z (m)", xscale=log10)

title = @sprintf("Single column simulation at %.2f, %.2f", φ★, λ★)
Label(fig[0, 1:6], title)

slider = Slider(fig[4, 1:6], range=1:Nt, startvalue=1)
n = slider.value

times = (times .- times[1]) ./days
Nt = length(times)
tn = @lift times[$n]

colors = Makie.wong_colors()

ρₒ = coupled_model.fluxes.ocean_reference_density
τx = interior(ρτx, 1, 1, 1, :) ./ ρₒ
τy = interior(ρτy, 1, 1, 1, :) ./ ρₒ
u★ = @. (τx^2 + τy^2)^(1/4)

lines!(axu, times, interior(u, 1, 1, Nz, :), color=colors[1], label="Zonal")
lines!(axu, times, interior(v, 1, 1, Nz, :), color=colors[2], label="Meridional")
lines!(axu, times, u★, color=colors[3], label="Ocean-side u★") 
vlines!(axu, tn, linewidth=4, color=(:black, 0.5))
axislegend(axu)

lines!(axτ, times, interior(ρτx, 1, 1, 1, :), label="Zonal")
lines!(axτ, times, interior(ρτy, 1, 1, 1, :), label="Meridional")
vlines!(axτ, tn, linewidth=4, color=(:black, 0.5))
axislegend(axτ)

lines!(axT, times, Tat[1:Nt] .- 273.15,      color=colors[1], linewidth=2, linestyle=:dash, label="Atmosphere temperature")
lines!(axT, times, interior(T, 1, 1, Nz, :), color=colors[2], linewidth=4, label="Ocean surface temperature")
vlines!(axT, tn, linewidth=4, color=(:black, 0.5))
axislegend(axT)

lines!(axQ, times, interior(Qv, 1, 1, 1, 1:Nt),    color=colors[2], label="Sensible",  linewidth=2)
lines!(axQ, times, interior(Qc, 1, 1, 1, 1:Nt),    color=colors[3], label="Latent",    linewidth=2)
lines!(axQ, times, - interior(Qsw, 1, 1, 1, 1:Nt), color=colors[4], label="Shortwave", linewidth=2)
lines!(axQ, times, - interior(Qlw, 1, 1, 1, 1:Nt), color=colors[5], label="Longwave",  linewidth=2)
vlines!(axQ, tn, linewidth=4, color=(:black, 0.5))
axislegend(axQ)

lines!(axF, times, Pt[1:Nt], label="Prescribed freshwater flux")
lines!(axF, times, - interior(Ev, 1, 1, 1, 1:Nt), label="Evaporation")
vlines!(axF, tn, linewidth=4, color=(:black, 0.5))
axislegend(axF)

lines!(axS, times, interior(S, 1, 1, Nz, :))
vlines!(axS, tn, linewidth=4, color=(:black, 0.5))

zc = znodes(T)
zf = znodes(κ)
un  = @lift interior(u[$n],  1, 1, :)
vn  = @lift interior(v[$n],  1, 1, :)
Tn  = @lift interior(T[$n],  1, 1, :)
Sn  = @lift interior(S[$n],  1, 1, :)
κn  = @lift interior(κ[$n],  1, 1, :)
en  = @lift interior(e[$n],  1, 1, :)
N²n = @lift interior(N²[$n], 1, 1, :)

scatterlines!(axuz, un,  zc, label="u") 
scatterlines!(axuz, vn,  zc, label="v") 
scatterlines!(axTz, Tn,  zc) 
scatterlines!(axSz, Sn,  zc) 
scatterlines!(axez, en,  zc) 
scatterlines!(axNz, N²n, zf) 
scatterlines!(axκz, κn,  zf) 

axislegend(axuz)

Tmax = maximum(interior(T))
Tmin = minimum(interior(T))
xlims!(axTz, Tmin - 0.1, Tmax + 0.1)

Nmax = maximum(interior(N²))
Nmin = minimum(interior(N²))
xlims!(axNz, Nmin / 2, Nmin * 1.1)

emax = maximum(interior(e))
xlims!(axez, 8e-7, emax * 1.1)
xlims!(axκz, 1e-7, 10)

Smax = maximum(interior(S))
Smin = minimum(interior(S))
xlims!(axSz, Smin - 0.2, Smax + 0.2)

record(fig, "single_column_profiles.mp4", 1:8:Nt, framerate=24) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end
nothing #hide

# ![](single_column_profiles.mp4)
