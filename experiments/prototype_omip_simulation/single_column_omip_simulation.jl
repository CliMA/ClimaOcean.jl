using Oceananigans
using Oceananigans.Units
using Oceananigans.BuoyancyModels: buoyancy_frequency

using ClimaOcean
using ClimaOcean.OceanSeaIceModels: Radiation
using ClimaOcean.DataWrangling.JRA55: jra55_prescribed_atmosphere
using ClimaOcean.DataWrangling.ECCO2: ecco2_field

using GLMakie
using Printf
using Dates

start_time = time_ns()

include("single_column_omip_ocean_component.jl")

epoch = Date(1992, 1, 1)
date = Date(1992, 10, 01)
start_seconds = Second(date - epoch).value
uᵢ = ecco2_field(:u_velocity, date)
vᵢ = ecco2_field(:v_velocity, date)
Tᵢ = ecco2_field(:temperature, date)
Sᵢ = ecco2_field(:salinity, date)

land = interior(Tᵢ) .< -10
interior(Tᵢ)[land] .= NaN
interior(Sᵢ)[land] .= NaN

teos10 = TEOS10EquationOfState()
buoyancy = SeawaterBuoyancy(equation_of_state=teos10)
tracers = (T=Tᵢ, S=Sᵢ)
N²_op = buoyancy_frequency(buoyancy, Tᵢ.grid, tracers)
N² = Field(N²_op)
compute!(N²)

elapsed = time_ns() - start_time
@info "Initial condition built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

#####
##### Construct the grid
#####

zc = znodes(Tᵢ)
zf = znodes(N²)

arch = CPU()

Δ = 1/4 # resolution in degrees
φ₁ = -90 + Δ/2
φ₂ = +90 - Δ/2
λ₁ = 0   + Δ/2
λ₂ = 360 - Δ/2
φe = φ₁:Δ:φ₂
λe = λ₁:Δ:λ₂

φ★ = 50 # degrees latitude
λ★ = 180 + 35 # degrees longitude (?)

i★ = searchsortedfirst(λe, λ★)
j★ = searchsortedfirst(φe, φ★)

longitude = (λe[i★] - Δ/2, λe[i★] + Δ/2)
latitude  = (φe[j★] - Δ/2, φe[j★] + Δ/2)

# Column
uc = interior(uᵢ, i★:i★, j★:j★, :)
vc = interior(vᵢ, i★:i★, j★:j★, :)
Tc = interior(Tᵢ, i★:i★, j★:j★, :)
Sc = interior(Sᵢ, i★:i★, j★:j★, :)

# Find bottom
zm = -400
zf = znodes(Tᵢ.grid, Face())
kb = findlast(T -> T < -20, Tc[1, 1, :])
km = findlast(z -> z < zm, zf)
k★ = isnothing(kb) ? km : max(kb, km)

Nz = size(Tc, 3)
kf = k★:Nz+1
kc = k★:Nz
zf = zf[kf]
uc = uc[:, :, kc]
vc = vc[:, :, kc]
Tc = Tc[:, :, kc]
Sc = Sc[:, :, kc]
Nz′ = length(kc)

grid = LatitudeLongitudeGrid(arch; longitude, latitude,
                             size = (1, 1, Nz′),
                             z = zf,
                             topology = (Periodic, Periodic, Bounded))

elapsed = time_ns() - start_time
@info "Grid constructed. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

ocean = omip_ocean_component(grid)
elapsed = time_ns() - start_time
@info "Ocean component built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

Ndays = 365
Nt = 8 * Ndays
atmosphere = jra55_prescribed_atmosphere(grid, 1:Nt) #, 1:21)
elapsed = time_ns() - start_time
@info "Atmosphere built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

ocean.model.clock.time = start_seconds
ocean.model.clock.iteration = 0
set!(ocean.model, T=Tc, S=Sc, e=1e-6)

ua = atmosphere.velocities.u
va = atmosphere.velocities.v
Ta = atmosphere.tracers.T
qa = atmosphere.tracers.q
times = ua.times

#=
fig = Figure(resolution=(1200, 1800))
axu = Axis(fig[1, 1])
axT = Axis(fig[2, 1])
axq = Axis(fig[3, 1])

lines!(axu, times ./ days, interior(ua, 1, 1, 1, :))
lines!(axu, times ./ days, interior(va, 1, 1, 1, :))
lines!(axT, times ./ days, interior(Ta, 1, 1, 1, :))
lines!(axq, times ./ days, interior(qa, 1, 1, 1, :))

display(fig)
=#

sea_ice = nothing
radiation = Radiation()
coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

#=
coupled_simulation = Simulation(coupled_model, Δt=10minutes, stop_time=start_seconds + 90days)

elapsed = time_ns() - start_time
@info "Coupled simulation built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

wall_clock = Ref(time_ns())

function progress(sim)
    msg1 = string("Iter: ", iteration(sim), ", time: ", prettytime(sim))

    elapsed = 1e-9 * (time_ns() - wall_clock[])
    msg2 = string(", wall time: ", prettytime(elapsed))
    wall_clock[] = time_ns()

    u, v, w = sim.model.ocean.model.velocities
    msg3 = @sprintf(", max|u|: (%.2e, %.2e)", maximum(abs, u), maximum(abs, v))

    T = sim.model.ocean.model.tracers.T
    S = sim.model.ocean.model.tracers.S
    e = sim.model.ocean.model.tracers.e
    Nz = size(T, 3)
    msg4 = @sprintf(", T₀: %.2f ᵒC",     first(interior(T, 1, 1, Nz)))
    msg5 = @sprintf(", S₀: %.2f g/kg",   first(interior(S, 1, 1, Nz)))
    msg6 = @sprintf(", e₀: %.2e m² s⁻²", first(interior(e, 1, 1, Nz)))

    @info msg1 * msg2 * msg3 * msg4 * msg5 * msg6
end

coupled_simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

# Build flux outputs
Jᵘ = coupled_model.surfaces.ocean.momentum.u
Jᵛ = coupled_model.surfaces.ocean.momentum.v
Jᵀ = coupled_model.surfaces.ocean.tracers.T
F  = coupled_model.surfaces.ocean.tracers.S
ρₒ = coupled_model.ocean_reference_density
cₚ = coupled_model.ocean_heat_capacity

Q = ρₒ * cₚ * Jᵀ
τˣ = ρₒ * Jᵘ
τʸ = ρₒ * Jᵛ

fluxes = (; τˣ, τʸ, Q, F)
fields = merge(ocean.model.velocities, ocean.model.tracers)

# Slice fields at the surface
outputs = merge(fields, fluxes)

filename = "single_column_omip_surface_fields.jld2"

coupled_simulation.output_writers[:surface] = JLD2OutputWriter(ocean.model, outputs; filename,
                                                               schedule = TimeInterval(1hours),
                                                               overwrite_existing = true)

run!(coupled_simulation)

ut = FieldTimeSeries(filename, "u")
vt = FieldTimeSeries(filename, "v")
Tt = FieldTimeSeries(filename, "T")
St = FieldTimeSeries(filename, "S")
et = FieldTimeSeries(filename, "e")
Qt = FieldTimeSeries(filename, "Q")
Ft = FieldTimeSeries(filename, "F")
τˣt = FieldTimeSeries(filename, "τˣ")
τʸt = FieldTimeSeries(filename, "τʸ")

Nz = size(Tt, 3)
times = Qt.times

ua = atmosphere.velocities.u
va = atmosphere.velocities.v
Ta = atmosphere.tracers.T
qa = atmosphere.tracers.q
Qlw = atmosphere.downwelling_radiation.longwave
Qsw = atmosphere.downwelling_radiation.shortwave

using Oceananigans.Units: Time

Nt = length(times)
uat = zeros(Nt)
vat = zeros(Nt)
Tat = zeros(Nt)
qat = zeros(Nt)
Qswt = zeros(Nt)
Qlwt = zeros(Nt)

for n = 1:Nt
    t = times[n]
    uat[n]  = ua[1, 1, 1, Time(t)]
    vat[n]  = va[1, 1, 1, Time(t)]
    Tat[n]  = Ta[1, 1, 1, Time(t)]
    qat[n]  = qa[1, 1, 1, Time(t)]
    Qswt[n] = Qsw[1, 1, 1, Time(t)]
    Qlwt[n] = Qlw[1, 1, 1, Time(t)]
end

fig = Figure(resolution=(2400, 1800))

axu = Axis(fig[1, 1:4], xlabel="Time (days)", ylabel="Velocities (m s⁻¹)")
axτ = Axis(fig[2, 1:4], xlabel="Time (days)", ylabel="Wind stress (N m⁻²)")
axT = Axis(fig[3, 1:4], xlabel="Time (days)", ylabel="Temperature (K)")
axQ = Axis(fig[4, 1:4], xlabel="Time (days)", ylabel="Heat flux (W m⁻²)")
axF = Axis(fig[5, 1:4], xlabel="Time (days)", ylabel="Salt flux (...)")

axuz = Axis(fig[6, 1], xlabel="Velocities (m s⁻¹)", ylabel="z (m)")
axTz = Axis(fig[6, 2], xlabel="Temperature (K)", ylabel="z (m)")
axSz = Axis(fig[6, 3], xlabel="Salinity (g/kg)", ylabel="z (m)")
axez = Axis(fig[6, 4], xlabel="Turbulent kinetic energy (m² s⁻²)", ylabel="z (m)")

slider = Slider(fig[7, 1:4], range=1:Nt, startvalue=1)
n = slider.value

times ./= days
tn = @lift times[$n]

colors = Makie.wong_colors()

lines!(axu, times, uat, color=colors[1])
lines!(axu, times, interior(ut, 1, 1, Nz, :), color=colors[1], linestyle=:dash)
vlines!(axu, tn)

lines!(axu, times, vat, color=colors[2])
lines!(axu, times, interior(vt, 1, 1, Nz, :), color=colors[2], linestyle=:dash)

lines!(axτ, times, interior(τˣt, 1, 1, 1, :))
lines!(axτ, times, interior(τʸt, 1, 1, 1, :))
vlines!(axτ, tn)

lines!(axT, times, Tat, color=colors[1])
lines!(axT, times, interior(Tt, 1, 1, Nz, :) .+ 273.15, color=colors[1], linestyle=:dash)
vlines!(axT, tn)

lines!(axQ, times, interior(Qt, 1, 1, 1, :), color=colors[1], linestyle=:dash)
lines!(axQ, times, - Qswt, color=colors[2], linewidth=3)
lines!(axQ, times, - Qlwt, color=colors[3], linewidth=3)
vlines!(axQ, tn)

lines!(axF, times, interior(Ft, 1, 1, 1, :))
vlines!(axF, tn)

z = znodes(Tt)
un = @lift interior(ut[$n], 1, 1, :)
vn = @lift interior(vt[$n], 1, 1, :)
Tn = @lift interior(Tt[$n], 1, 1, :)
Sn = @lift interior(St[$n], 1, 1, :)
en = @lift interior(et[$n], 1, 1, :)
lines!(axuz, un, z) 
lines!(axuz, vn, z) 
lines!(axTz, Tn, z) 
lines!(axSz, Sn, z) 
lines!(axez, en, z) 

display(fig)
=#
