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

#####
##### Construct initial conditions + grid
#####

epoch = Date(1992, 1, 1)
date = Date(1992, 1, 2)
start_seconds = Second(date - epoch).value
Te = ecco2_field(:temperature, date)
Se = ecco2_field(:salinity, date)

latitude = (-30, 30)
grid, (Tᵢ, Sᵢ) = regional_ecco2_grid(arch, Te, Se; latitude)

Nt = 8 * 30
atmosphere = JRA55_prescribed_atmosphere(arch, 1:Nt; backend=InMemory(8))
radiation = Radiation()
sea_ice = nothing

#closure = RiBasedVerticalDiffusivity(maximum_diffusivity=1e2, maximum_viscosity=1e2)
#closure = RiBasedVerticalDiffusivity()
closure = :default
ocean = omip_ocean_component(grid; closure)
set!(ocean.model, T=Tᵢ, S=Sᵢ)

if :e ∈ keys(ocean.model.tracers)
    set!(ocean.model, e=1e-6)
end

coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
coupled_model.clock.time = start_seconds
stop_time = start_seconds + 30days
coupled_simulation = Simulation(coupled_model, Δt=5minutes, stop_time=stop_time)

elapsed = 1e-9 * (time_ns() - start_time)
@info string("Coupled simulation built ", prettytime(elapsed), ".")
start_time = time_ns()

wall_clock = Ref(time_ns())

# Hm...

function clip_diffusivity(coupled_simulation)
    ocean_model = coupled_simulation.model.ocean.model
    κᶜ = parent(ocean_model.diffusivity_fields.κᶜ)
    κᵘ = parent(ocean_model.diffusivity_fields.κᵘ)
    @. κᶜ = min(κᶜ, 100)
    @. κᵘ = min(κᵘ, 100)
    return nothing
end

coupled_simulation.callbacks[:clip] = Callback(clip_diffusivity)

#####
##### Progress
#####

Jᵘ = coupled_model.fluxes.total.ocean.momentum.u
Jᵛ = coupled_model.fluxes.total.ocean.momentum.v
Jᵀ = coupled_model.fluxes.total.ocean.tracers.T
Jˢ = coupled_model.fluxes.total.ocean.tracers.S
Fv = coupled_model.fluxes.turbulent.fields.water_vapor
Qc = coupled_model.fluxes.turbulent.fields.sensible_heat
Qv = coupled_model.fluxes.turbulent.fields.latent_heat
ρₒ = coupled_model.fluxes.ocean_reference_density
cₚ = coupled_model.fluxes.ocean_heat_capacity

import Oceananigans.Fields: reduced_dimensions
reduced_dimensions(::Oceananigans.AbstractOperations.BinaryOperation) = tuple(3)

ΣQ = ρₒ * cₚ * Jᵀ
τˣ = ρₒ * Jᵘ
τʸ = ρₒ * Jᵛ
N² = buoyancy_frequency(ocean.model)
κᶜ = ocean.model.diffusivity_fields.κᶜ

function progress(sim)
    msg = string("Iter: ", iteration(sim), ", time: ", prettytime(sim))

    elapsed = 1e-9 * (time_ns() - wall_clock[])
    msg *= string(", wall time: ", prettytime(elapsed))
    wall_clock[] = time_ns()

    u, v, w = sim.model.ocean.model.velocities
    msg *= @sprintf(", max|u|: (%.2e, %.2e)", maximum(abs, u), maximum(abs, v))

    T = sim.model.ocean.model.tracers.T
    S = sim.model.ocean.model.tracers.S

    msg *= @sprintf(", extrema(T): (%.2f, %.2f) ᵒC",          minimum(T),  maximum(T))
    msg *= @sprintf(", extrema(S): (%.2f, %.2f) g kg⁻¹",      minimum(S),  maximum(S))
    msg *= @sprintf(", max|τ|: (%.2e, %.2e) N m⁻²",           maximum(τˣ), maximum(τʸ))
    msg *= @sprintf(", max|Qv|: %.2e W m⁻²",                  maximum(Qv))
    msg *= @sprintf(", max|Qc|: %.2e W m⁻²",                  maximum(Qc))
    msg *= @sprintf(", extrema(ΣQ): (%.2e, %.2e) W m⁻²",      minimum(ΣQ), maximum(ΣQ))
    msg *= @sprintf(", extrema(Fv): (%.2e, %.2e) kg s⁻¹ m⁻²", minimum(Fv), maximum(Fv))
    msg *= @sprintf(", extrema(κᶜ): (%.2e, %.2e) m² s⁻¹",     minimum(κᶜ), maximum(κᶜ))

    #e = sim.model.ocean.model.tracers.e
    #msg *= @sprintf(", extrema(e): (%.2f, %.2f) m² s⁻²", minimum(e), maximum(e))

    @info msg
end

coupled_simulation.callbacks[:progress] = Callback(progress, IterationInterval(1))

run!(coupled_simulation)

#=
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

coupled_simulation.output_writers[:fluxes] = JLD2OutputWriter(ocean.model, fluxes;
                                                              filename = filename * "_fluxes",
                                                              schedule = TimeInterval(1days),
                                                              overwrite_existing = true)

coupled_simulation.output_writers[:fields] = JLD2OutputWriter(ocean.model, fields;
                                                              filename = filename * "_fields",
                                                              indices = (:, :, Nz),
                                                              schedule = TimeInterval(1days),
                                                              overwrite_existing = true)

=#
