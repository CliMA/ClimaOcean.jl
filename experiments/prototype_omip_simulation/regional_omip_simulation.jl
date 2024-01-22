using Oceananigans
using Oceananigans.Units
using Oceananigans.BuoyancyModels: buoyancy_frequency
using Oceananigans.Units: Time

using ClimaOcean
using ClimaOcean.OceanSeaIceModels: Radiation
using ClimaOcean.DataWrangling.JRA55: jra55_prescribed_atmosphere
using ClimaOcean.DataWrangling.ECCO2: ecco2_field

using GLMakie
using Printf
using Dates

start_time = time_ns()

include("omip_ocean_component.jl")

epoch = Date(1992, 1, 1)
date = Date(1992, 10, 1)
start_seconds = Second(date - epoch).value
# uᵢ = ecco2_field(:u_velocity, date)
# vᵢ = ecco2_field(:v_velocity, date)
Te = ecco2_field(:temperature, date)
Se = ecco2_field(:salinity, date)

land = interior(Te) .< -10
interior(Te)[land] .= NaN
interior(Se)[land] .= NaN

elapsed = time_ns() - start_time
@info "Initial condition built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

#####
##### Construct the grid
#####

arch = CPU()

latitude = (-60, -50)
longitude = (300, 360)

i₁ = 4 * first(longitude) + 1
i₂ = 1440 - 4 * (360 - last(longitude))
Nx = i₂ - i₁ + 1

j₁ = 4 * (90 + first(latitude)) + 1
j₂ = 720 - 4 * (90 - last(latitude))
Ny = j₂ - j₁ + 1

zc = znodes(Te)
zf = znodes(Te.grid, Face())
Δz = first(zspacings(Te.grid, Center()))

Tᵢ = interior(Te, i₁:i₂, j₁:j₂, :)
Sᵢ = interior(Se, i₁:i₂, j₁:j₂, :)

# Construct bottom_height depth by analyzing T
Nx, Ny, Nz = size(Tᵢ)
bottom_height = ones(Nx, Ny) .* (zf[1] - Δz)

for i = 1:Nx, j = 1:Ny
    @inbounds for k = Nz:-1:1
        if isnan(Tᵢ[i, j, k])
            bottom_height[i, j] = zf[k+1]
            break
        end
    end
end

@show Nx Ny Nz zf

grid = LatitudeLongitudeGrid(arch; latitude, longitude,
                             size = (Nx, Ny, Nz),
                             halo = (7, 7, 7),
                             z = zf,
                             topology = (Periodic, Bounded, Bounded))

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

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
set!(ocean.model, T=Tᵢ, S=Sᵢ, e=1e-6)

ua = atmosphere.velocities.u
va = atmosphere.velocities.v
Ta = atmosphere.tracers.T
qa = atmosphere.tracers.q
times = ua.times

sea_ice = nothing
radiation = Radiation()
coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

coupled_simulation = Simulation(coupled_model, Δt=10minutes, stop_time=start_seconds + 60days)

elapsed = time_ns() - start_time
@info "Coupled simulation built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

wall_clock = Ref(time_ns())

function progress(sim)
    msg = string("(", location, ")")
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
    u★ = (τˣ^2 + τʸ^2)^(1/4)
    Q = first(sim.model.fluxes.total.ocean.heat)

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
Jᵘ  = coupled_model.fluxes.total.ocean.momentum.u
Jᵛ  = coupled_model.fluxes.total.ocean.momentum.v
Jᵀ  = coupled_model.fluxes.total.ocean.tracers.T
F   = coupled_model.fluxes.total.ocean.tracers.S
E   = coupled_model.fluxes.turbulent.fields.evaporation
Qse = coupled_model.fluxes.turbulent.fields.sensible_heat_flux
Qla = coupled_model.fluxes.turbulent.fields.latent_heat_flux
ρₒ  = coupled_model.fluxes.ocean_reference_density
cₚ  = coupled_model.fluxes.ocean_heat_capacity

Q = ρₒ * cₚ * Jᵀ
τˣ = ρₒ * Jᵘ
τʸ = ρₒ * Jᵛ
N² = buoyancy_frequency(ocean.model)
κᶜ = ocean.model.diffusivity_fields.κᶜ

fluxes = (; τˣ, τʸ, E, F, Q, Qse, Qla)

auxiliary_fields = (; N², κᶜ)
fields = merge(ocean.model.velocities, ocean.model.tracers, auxiliary_fields)

# Slice fields at the surface
outputs = merge(fields, fluxes)

filename = "single_column_omip_$location"

coupled_simulation.output_writers[:jld2] = JLD2OutputWriter(ocean.model, outputs; filename,
                                                            schedule = TimeInterval(3hours),
                                                            overwrite_existing = true)

#=
coupled_simulation.output_writers[:nc] = NetCDFOutputWriter(ocean.model, outputs; filename,
                                                            schedule = AveragedTimeInterval(1days),
                                                            overwrite_existing = true)
=#

# run!(coupled_simulation)

