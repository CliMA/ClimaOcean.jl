using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: node
using Oceananigans.BuoyancyModels: buoyancy_frequency
using Oceananigans.Units: Time

using ClimaOcean
using ClimaOcean.OceanSeaIceModels: Radiation
using ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: interp_atmos_time_series
using ClimaOcean.DataWrangling.JRA55: JRA55_prescribed_atmosphere
using ClimaOcean.DataWrangling.ECCO2: ecco2_field

using GLMakie
using Printf
using Dates

include("omip_components.jl")

locations = (
    #eastern_mediterranean = (λ =  30, φ = 32), 
    ocean_station_papa = (λ = 215, φ = 50), 
    north_atlantic = (λ = 325, φ = 50), 
    drake_passage = (λ = 300, φ = -60), 
    weddell_sea = (λ = 325, φ = -70), 
    tasman_southern_ocean = (λ = 145, φ = -55), 
)

location = :ocean_station_papa

start_time = time_ns()

epoch = Date(1992, 1, 1)
date = Date(1992, 10, 1)
start_seconds = Second(date - epoch).value
uᵢ = ecco2_field(:u_velocity, date)
vᵢ = ecco2_field(:v_velocity, date)
Tᵢ = ecco2_field(:temperature, date)
Sᵢ = ecco2_field(:salinity, date)

land = interior(Tᵢ) .< -10
interior(Tᵢ)[land] .= NaN
interior(Sᵢ)[land] .= NaN

elapsed = time_ns() - start_time
@info "Initial condition built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

#####
##### Construct the grid
#####

arch = CPU()

Δ = 1/4 # resolution in degrees
φ₁ = -90 + Δ/2
φ₂ = +90 - Δ/2
λ₁ = 0   + Δ/2
λ₂ = 360 - Δ/2
φe = φ₁:Δ:φ₂
λe = λ₁:Δ:λ₂

λ★, φ★ = locations[location]

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
k★ = isnothing(kb) ? km : max(kb + 1, km)

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

Ndays = 7
Nt = 8 * Ndays
atmosphere = JRA55_prescribed_atmosphere(1:Nt, backend=InMemory(8)) #, 1:21)
elapsed = time_ns() - start_time
@info "Atmosphere built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

# ocean.model.clock.time = start_seconds
ocean.model.clock.iteration = 0
set!(ocean.model, T=Tc, S=Sc, e=1e-6)

ua = atmosphere.velocities.u
va = atmosphere.velocities.v
Ta = atmosphere.tracers.T
qa = atmosphere.tracers.q
times = ua.times

sea_ice = nothing
radiation = Radiation()
coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

coupled_model.clock.iteration = 0
coupled_model.clock.time = 0
set!(coupled_model.ocean.model, u=0, v=0, T=Tc, S=Sc, e=1e-6)
coupled_simulation = Simulation(coupled_model, Δt=10minutes, stop_time=14days)

elapsed = time_ns() - start_time
@info "Coupled simulation built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

wall_clock = Ref(time_ns())

atmos_grid = atmosphere.grid
const c = Center()

function progress(sim)
    msg = string("Iter: ", iteration(sim), ", time: ", prettytime(sim))

    #=
    elapsed = 1e-9 * (time_ns() - wall_clock[])
    msg *= string(", wall time: ", prettytime(elapsed))
    wall_clock[] = time_ns()
    =#

    t = time(sim)
    X = node(1, 1, 1, sim.model.ocean.model.grid, c, c, c)
    uai = interp_atmos_time_series(ua, X, Time(t), atmos_grid) 
    vai = interp_atmos_time_series(va, X, Time(t), atmos_grid) 
    Tai = interp_atmos_time_series(Ta, X, Time(t), atmos_grid) 
    qai = interp_atmos_time_series(qa, X, Time(t), atmos_grid) 

    msg *= @sprintf(", ua: %.2e, va: %.2e, Ta: %.2f, qa: %.2e", uai, vai, Tai, qai)

    u, v, w = sim.model.ocean.model.velocities
    msg *= @sprintf(", max|u|: (%.2e, %.2e)", maximum(abs, u), maximum(abs, v))

    T = sim.model.ocean.model.tracers.T
    S = sim.model.ocean.model.tracers.S
    e = sim.model.ocean.model.tracers.e

    τˣ = first(sim.model.fluxes.total.ocean.momentum.τˣ)
    τʸ = first(sim.model.fluxes.total.ocean.momentum.τʸ)
    u★ = (τˣ^2 + τʸ^2)^(1/4)

    Q = first(sim.model.fluxes.total.ocean.heat)

    t = time(sim)

    Nz = size(T, 3)
    msg *= @sprintf(", u★: %.2f m s⁻¹", u★)
    msg *= @sprintf(", Q: %.2f W m⁻²", Q)

    #=
    msg *= @sprintf(", T₀: %.2f ᵒC",     first(interior(T, 1, 1, Nz)))
    msg *= @sprintf(", extrema(T): (%.2f, %.2f) ᵒC", minimum(T), maximum(T))
    msg *= @sprintf(", S₀: %.2f g/kg",   first(interior(S, 1, 1, Nz)))
    msg *= @sprintf(", e₀: %.2e m² s⁻²", first(interior(e, 1, 1, Nz)))
    =#

    @info msg
end

coupled_simulation.callbacks[:progress] = Callback(progress, IterationInterval(1))

run!(coupled_simulation)

