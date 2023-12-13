using Oceananigans
using Oceananigans.Units

using ClimaOcean
using ClimaOcean.OceanSeaIceModels: SurfaceRadiation
using ClimaOcean.DataWrangling.JRA55: jra55_prescribed_atmosphere
using ClimaOcean.DataWrangling.ECCO2: ecco2_field

using NCDatasets
using GLMakie
using Printf

using Downloads: download

start_time = time_ns()

include("single_column_omip_ocean_component.jl")

Tᵢ = ecco2_field(:temperature)
Sᵢ = ecco2_field(:salinity)

elapsed = time_ns() - start_time
@info "Initial condition built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

#####
##### Construct the grid
#####

z = znodes(Tᵢ)

arch = CPU()

Δ = 1/4 # resolution in degrees
φ₁ = -90 + Δ/2
φ₂ = +90 - Δ/2
λ₁ = 0   + Δ/2
λ₂ = 360 - Δ/2
φe = φ₁:Δ:φ₂
λe = λ₁:Δ:λ₂

land = interior(Tᵢ) .< -10
interior(Tᵢ)[land] .= NaN
interior(Sᵢ)[land] .= NaN

Nz = size(Tᵢ, 3)
fig = Figure(resolution=(1200, 1200))
map = Axis(fig[1, 1:2], xlabel="λ (degrees)", ylabel="φ (degrees)")
hm = heatmap!(map, λe, φe, interior(Tᵢ, :, :, Nz), colorrange=(0, 30), nan_color=:gray)
Colorbar(fig[1, 3], hm, label="Surface temperature (ᵒC)")

axT = Axis(fig[2, 1], ylabel="z (m)", xlabel="Temperature (ᵒC)")
axS = Axis(fig[2, 2], ylabel="z (m)", xlabel="Salinity (psu)")

φs = [50,   55, 0,   -30]
λs = [215, 310, 210, 160]
Nc = length(φs)

for n = 1:Nc
    local φ★
    local λ★
    local i★
    local j★

    φ★ = φs[n]
    λ★ = λs[n]

    i★ = searchsortedfirst(λe, λ★)
    j★ = searchsortedfirst(φe, φ★)

    scatter!(map, λ★, φ★, strokewidth=4, strokecolor=:black,
             color=:pink, markersize=20)

    label = string("λ = ", λ★, ", φ = ", φ★)
    scatterlines!(axT, interior(Tᵢ, i★, j★, :), z; label)
    scatterlines!(axS, interior(Sᵢ, i★, j★, :), z; label)
end

xlims!(axT, 0, 30)
xlims!(axS, 32, 36)
ylims!(axT, -2000, 30)
ylims!(axS, -2000, 30)
axislegend(axT, position=:rb)

display(fig)

φ★ = 50 # degrees latitude
λ★ = 180 + 35 # degrees longitude (?)

i★ = searchsortedfirst(λe, λ★)
j★ = searchsortedfirst(φe, φ★)

longitude = (λe[i★] - Δ/2, λe[i★] + Δ/2)
latitude  = (φe[j★] - Δ/2, φe[j★] + Δ/2)

# Column
Tc = interior(Tᵢ, i★:i★, j★:j★, :)
Sc = interior(Sᵢ, i★:i★, j★:j★, :)

# Find bottom
zf = znodes(Tᵢ.grid, Face())
kb = findlast(T -> T < -20, Tc[1, 1, :])
km = findlast(z -> z < -2000, zf)
k★ = isnothing(kb) ? km : max(kb, km)

Nz = size(Tc, 3)
kf = k★:Nz+1
kc = k★:Nz
zf = zf[kf]
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

ocean.model.clock.time = 0
ocean.model.clock.iteration = 0
set!(ocean.model, T=Tc, S=Sc, e=1e-6)

days = 30
Nt = 8days
atmosphere = jra55_prescribed_atmosphere(grid, 1:Nt) #, 1:21)
elapsed = time_ns() - start_time
@info "Atmosphere built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

ua = atmosphere.velocities.u
va = atmosphere.velocities.v
Ta = atmosphere.tracers.T
qa = atmosphere.tracers.q
times = ua.times

fig = Figure()
axu = Axis(fig[1, 1])
axT = Axis(fig[2, 1])
axq = Axis(fig[3, 1])

lines!(axu, times ./ day, interior(ua, 1, 1, 1, :))
lines!(axu, times ./ day, interior(va, 1, 1, 1, :))
lines!(axT, times ./ day, interior(Ta, 1, 1, 1, :))
lines!(axq, times ./ day, interior(qa, 1, 1, 1, :))

display(fig)

sea_ice = nothing
surface_radiation = SurfaceRadiation()
coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, surface_radiation)
coupled_simulation = Simulation(coupled_model, Δt=5minutes, stop_time=30days)

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
    msg5 = @sprintf(", S₀: %.2f psu",    first(interior(S, 1, 1, Nz)))
    msg6 = @sprintf(", e₀: %.2e m² s⁻²", first(interior(e, 1, 1, Nz)))

    @info msg1 * msg2 * msg3 * msg4 * msg5 * msg6
end

coupled_simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

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
                                                               schedule = TimeInterval(1hour),
                                                               overwrite_existing = true)

run!(coupled_simulation)

Tt = FieldTimeSeries(filename, "T")
St = FieldTimeSeries(filename, "S")
Qt = FieldTimeSeries(filename, "Q")
Ft = FieldTimeSeries(filename, "F")
τˣt = FieldTimeSeries(filename, "τˣ")
τʸt = FieldTimeSeries(filename, "τʸ")

times = Qt.times

fig = Figure()
axτ = Axis(fig[1, 1])
axQ = Axis(fig[2, 1])
axF = Axis(fig[3, 1])

lines!(axτ, times, interior(τˣt, 1, 1, 1, :))
lines!(axτ, times, interior(τʸt, 1, 1, 1, :))
lines!(axQ, times, interior(Qt, 1, 1, 1, :))
lines!(axF, times, interior(Ft, 1, 1, 1, :))

display(fig)

