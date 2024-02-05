using Oceananigans
using Oceananigans.Architectures: arch_array
using Oceananigans.Units
using Oceananigans.BuoyancyModels: buoyancy_frequency
using Oceananigans.Units: Time

using ClimaOcean
using ClimaOcean.OceanSeaIceModels: Radiation
using ClimaOcean.DataWrangling.JRA55: JRA55_prescribed_atmosphere
using ClimaOcean.DataWrangling.ECCO2: ecco2_field

# using GLMakie
using Printf
using Dates

start_time = time_ns()

include("omip_components.jl")

arch = GPU()
epoch = Date(1992, 1, 1)
date = Date(1992, 10, 1)
start_seconds = Second(date - epoch).value
# uᵢ = ecco2_field(:u_velocity, date)
# vᵢ = ecco2_field(:v_velocity, date)
Te = ecco2_field(:temperature, date)
Se = ecco2_field(:salinity, date)
ℋe = ecco2_field(:sea_ice_thickness, date)

land = interior(Te) .< -10
interior(Te)[land] .= NaN
interior(Se)[land] .= NaN

elapsed = time_ns() - start_time
@info "Initial condition built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

#####
##### Construct the grid
#####

latitude = (-75, -30)
longitude = (0, 360)

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
ℋᵢ = interior(ℋe, i₁:i₂, j₁:j₂, :)

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

Tᵢ = arch_array(arch, Tᵢ)
Sᵢ = arch_array(arch, Sᵢ)
ℋᵢ = arch_array(arch, ℋᵢ)

if longitude[2] - longitude[1] == 360
    TX = Periodic
else
    TX = Bounded
end

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

#=
Ndays = 30
Nt = 8 * Ndays
atmosphere = JRA55_prescribed_atmosphere(arch, 1:Nt; backend=InMemory(8))
elapsed = time_ns() - start_time
@info "Atmosphere built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()
=#

atmosphere = nothing

ocean.model.clock.time = start_seconds
ocean.model.clock.iteration = 0
set!(ocean.model, T=Tᵢ, S=Sᵢ, e=1e-6)

sea_ice = omip_sea_ice_component(ocean.model) #nothing
set!(sea_ice.model, h=ℋᵢ)

radiation = nothing # Radiation()
coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

stop_time = start_seconds + 90days
coupled_simulation = Simulation(coupled_model, Δt=10minutes, stop_time=stop_time)

elapsed = time_ns() - start_time
@info "Coupled simulation built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

wall_clock = Ref(time_ns())

function progress(sim)
    msg = string("Iter: ", iteration(sim), ", time: ", prettytime(sim))

    elapsed = 1e-9 * (time_ns() - wall_clock[])
    msg *= string(", wall time: ", prettytime(elapsed))
    wall_clock[] = time_ns()

    u, v, w = sim.model.ocean.model.velocities
    msg *= @sprintf(", max|u|: (%.2e, %.2e)", maximum(abs, u), maximum(abs, v))

    T = sim.model.ocean.model.tracers.T
    S = sim.model.ocean.model.tracers.S
    e = sim.model.ocean.model.tracers.e

    msg *= @sprintf(", extrema(T): (%.2f, %.2f) ᵒC", maximum(T), minimum(T))
    msg *= @sprintf(", extrema(S): (%.2f, %.2f) g kg⁻¹", minimum(S), maximum(S))
    msg *= @sprintf(", extrema(e): (%.2f, %.2f) m² s⁻²", minimum(e), maximum(e))

    @info msg
end

coupled_simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

# Build flux outputs
Jᵘ = coupled_model.fluxes.total.ocean.momentum.u
Jᵛ = coupled_model.fluxes.total.ocean.momentum.v
Jᵀ = coupled_model.fluxes.total.ocean.tracers.T
Jˢ = coupled_model.fluxes.total.ocean.tracers.S
Fv = coupled_model.fluxes.turbulent.fields.freshwater
Qc = coupled_model.fluxes.turbulent.fields.sensible_heat
Qv = coupled_model.fluxes.turbulent.fields.latent_heat
ρₒ = coupled_model.fluxes.ocean_reference_density
cₚ = coupled_model.fluxes.ocean_heat_capacity

ΣQ = ρₒ * cₚ * Jᵀ
τˣ = ρₒ * Jᵘ
τʸ = ρₒ * Jᵛ
N² = buoyancy_frequency(ocean.model)
κᶜ = ocean.model.diffusivity_fields.κᶜ

fluxes = (; τˣ, τʸ, Fv, Jˢ, ΣQ, Qc, Qv)

auxiliary_fields = (; N², κᶜ)
fields = merge(ocean.model.velocities, ocean.model.tracers, auxiliary_fields)

# Slice fields at the surface
outputs = merge(fields, fluxes)

filename = "regional_omip_simulation"

#=
coupled_simulation.output_writers[:fluxes] = JLD2OutputWriter(ocean.model, fluxes;
                                                              filename = filename * "_fluxes",
                                                              schedule = TimeInterval(1days),
                                                              overwrite_existing = true)
=#

coupled_simulation.output_writers[:fields] = JLD2OutputWriter(ocean.model, fields;
                                                              filename = filename * "_fields",
                                                              indices = (:, :, Nz),
                                                              schedule = TimeInterval(1days),
                                                              overwrite_existing = true)

coupled_simulation.output_writers[:seaice] = JLD2OutputWriter(sea_ice.model, (; h = sea_ice.model.ice_thickness);
                                                              filename = filename * "_sea_ice_thickness",
                                                              schedule = TimeInterval(1days),
                                                              overwrite_existing = true)

run!(coupled_simulation)
