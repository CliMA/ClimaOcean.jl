# # Single column ocean simulation forced by JRA55 Reananlysis
#
# In this example, we simulate the evolution of an ocean water column 
# forced by an atmosphere prescribed by the JRA55 Reananlysis data
# Specifically, the column is positioned at the location of the Ocean station
# Papa measurements (144.9ᵒ W and 50.1ᵒ N)
#
# ## Install dependencies
#
# First let's make sure we have all required packages installed.

# ```julia
# using Pkg
# pkg"add Oceananigans, ClimaOcean, CairoMakie"
# ```

using Oceananigans
using Oceananigans.Units
using Oceananigans.BuoyancyModels: buoyancy_frequency
using Oceananigans.Units: Time

using ClimaOcean
using ClimaOcean.ECCO2: ECCO2Metadata
using ClimaOcean.OceanSimulations

using CairoMakie
using Printf

#####
##### Construct the grid
#####

# Since it is a single column, and therefore computationally
# inexpensive, we can run the simulation entirely on the CPU
arch = CPU()

Nz = 80
H  = 400

# Ocean station papa location
λ★, φ★ = 35.1, 50.1 
longitude = λ★ .+ (-0.25, 0.25)
latitude  = φ★ .+ (-0.25, 0.25)

# We use a SingleColumnGrid
grid = LatitudeLongitudeGrid(; size = (3, 3, Nz),
                               longitude,
                               latitude,
                               z = (-H, 0),
                               topology = (Periodic, Periodic, Bounded))

# Building the ocean simulation
momentum_advection = nothing
tracer_advection = nothing
coriolis = FPlane(latitude = φ★)

ocean = ocean_simulation(grid; 
                         coriolis,
                         tracer_advection,
                         momentum_advection,
                         bottom_drag_coefficient = 0)
model = ocean.model

start_time = time_ns()

# Initial conditions
set!(ocean.model, T = ECCO2Metadata(:temperature),
                  S = ECCO2Metadata(:salinity),
                  e = 1e-6)

elapsed = time_ns() - start_time
@info "Initial condition built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()
            
# Retrieving the atmosphere
last_time = floor(Int, 31 * 24 / 3) # 31 days in hours divided by JRA55's frequency in hours
backend = InMemory()
atmosphere = JRA55_prescribed_atmosphere(time_indices = 1:last_time; 
                                         longitude, latitude, backend,
                                         include_rivers_and_icebergs = false)

ocean.model.clock.time = 0
ocean.model.clock.iteration = 0
ocean.Δt = 10minutes

ua = atmosphere.velocities.u
va = atmosphere.velocities.v
Ta = atmosphere.tracers.T
qa = atmosphere.tracers.q
times = ua.times

fig = Figure(size=(1200, 1800))
axu = Axis(fig[1, 1])
axT = Axis(fig[2, 1])
axq = Axis(fig[3, 1])

lines!(axu, times ./ days, interior(ua, 1, 1, 1, :))
lines!(axu, times ./ days, interior(va, 1, 1, 1, :))
lines!(axT, times ./ days, interior(Ta, 1, 1, 1, :))
lines!(axq, times ./ days, interior(qa, 1, 1, 1, :))

display(fig)

radiation = Radiation()
coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)
coupled_simulation = Simulation(coupled_model, Δt=10minutes, stop_time=30days)

elapsed = time_ns() - start_time
@info "Coupled simulation built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

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

    τˣ = first(sim.model.fluxes.total.ocean.momentum.τˣ)
    τʸ = first(sim.model.fluxes.total.ocean.momentum.τʸ)
    Q = first(sim.model.fluxes.total.ocean.heat)

    u★ = sqrt(sqrt(τˣ^2 + τʸ^2))

    Nz = size(T, 3)
    msg *= @sprintf(", u★: %.2f m s⁻¹", u★)
    msg *= @sprintf(", Q: %.2f W m⁻²", Q)
    msg *= @sprintf(", T₀: %.2f ᵒC",     first(interior(T, 1, 1, Nz)))
    msg *= @sprintf(", extrema(T): (%.2f, %.2f) ᵒC", minimum(T), maximum(T))
    msg *= @sprintf(", S₀: %.2f g/kg",   first(interior(S, 1, 1, Nz)))
    msg *= @sprintf(", e₀: %.2e m² s⁻²", first(interior(e, 1, 1, Nz)))

    @info msg
end

coupled_simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

# Build flux outputs
Ju = coupled_model.fluxes.total.ocean.momentum.u
Jv = coupled_model.fluxes.total.ocean.momentum.v
JT = coupled_model.fluxes.total.ocean.tracers.T
Js = coupled_model.fluxes.total.ocean.tracers.S
E  = coupled_model.fluxes.turbulent.fields.water_vapor
Qc = coupled_model.fluxes.turbulent.fields.sensible_heat
Qv = coupled_model.fluxes.turbulent.fields.latent_heat
ρₒ = coupled_model.fluxes.ocean_reference_density
cₚ = coupled_model.fluxes.ocean_heat_capacity

Q  = ρₒ * cₚ * JT
τx = ρₒ * Ju
τy = ρₒ * Jv
N² = buoyancy_frequency(ocean.model)
κc = ocean.model.diffusivity_fields.κc

fluxes = (; τx, τy, E, Js, Qv, Qc)

auxiliary_fields = (; N², κc)
fields = merge(ocean.model.velocities, ocean.model.tracers, auxiliary_fields)

# Slice fields at the surface
outputs = merge(fields, fluxes)

filename = "single_column_omip_$location"

coupled_simulation.output_writers[:jld2] = JLD2OutputWriter(ocean.model, outputs; filename,
                                                            schedule = TimeInterval(3hours),
                                                            overwrite_existing = true)

run!(coupled_simulation)

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
τˣ = FieldTimeSeries(filename, "τx")
τʸ = FieldTimeSeries(filename, "τy")

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
Jᵘ = interior(τˣ, 1, 1, 1, :) ./ ρₒ
Jᵛ = interior(τʸ, 1, 1, 1, :) ./ ρₒ
u★ = @. (Jᵘ^2 + Jᵛ^2)^(1/4)

lines!(axu, times, interior(u, 1, 1, Nz, :), color=colors[1], label="Zonal")
lines!(axu, times, interior(v, 1, 1, Nz, :), color=colors[2], label="Meridional")
lines!(axu, times, u★, color=colors[3], label="Ocean-side u★") 
vlines!(axu, tn, linewidth=4, color=(:black, 0.5))
axislegend(axu)

lines!(axτ, times, interior(τˣ, 1, 1, 1, :), label="Zonal")
lines!(axτ, times, interior(τʸ, 1, 1, 1, :), label="Meridional")
vlines!(axτ, tn, linewidth=4, color=(:black, 0.5))
axislegend(axτ)

lines!(axT, times, Tat[1:Nt] .- 273.15,             color=colors[1], linewidth=2, linestyle=:dash, label="Atmosphere temperature")
lines!(axT, times, interior(T, 1, 1, Nz, :), color=colors[2], linewidth=4, label="Ocean surface temperature")
vlines!(axT, tn, linewidth=4, color=(:black, 0.5))
axislegend(axT)

lines!(axQ, times, interior(Qv, 1, 1, 1, 1:Nt),    color=colors[2], label="Sensible",  linewidth=2)
lines!(axQ, times, interior(Qc, 1, 1, 1, 1:Nt),    color=colors[3], label="Latent",    linewidth=2)
lines!(axQ, times, - interior(Qsw, 1, 1, 1, 1:Nt), color=colors[4], label="Shortwave", linewidth=2)
lines!(axQ, times, - interior(Qlw, 1, 1, 1, 1:Nt), color=colors[5], label="Longwave",  linewidth=2)
vlines!(axQ, tn, linewidth=4, color=(:black, 0.5))
axislegend(axQ)

#lines!(axF, times, interior(Jˢt, 1, 1, 1, :), label="Net freshwater flux")
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

display(fig)

record(fig, "$(location)_single_column_simulation.mp4", 1:Nt, framerate=24) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end
