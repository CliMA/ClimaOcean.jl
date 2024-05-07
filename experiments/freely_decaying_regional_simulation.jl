using Oceananigans
using Oceananigans.Architectures: arch_array
using Oceananigans.Units
using Oceananigans.BuoyancyModels: buoyancy_frequency
using Oceananigans.Units: Time

using ClimaOcean
using ClimaOcean.OceanSeaIceModels: Radiation
using ClimaOcean.DataWrangling.JRA55: JRA55_prescribed_atmosphere
using ClimaOcean.DataWrangling.ECCO: ecco_field

# using GLMakie
using Printf
using Dates

start_time = time_ns()

include("omip_components.jl")

arch = CPU()
epoch = Date(1992, 1, 1)
date = Date(1992, 10, 1)
start_seconds = Second(date - epoch).value
Te = ecco_field(:temperature, date)
Se = ecco_field(:salinity, date)
# ℋe = ecco_field(:sea_ice_thickness, date)

land = interior(Te) .< -10
interior(Te)[land] .= NaN
interior(Se)[land] .= NaN

elapsed = time_ns() - start_time
@info "Initial condition built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

#####
##### Construct the grid
#####

latitude = (-80, -20)
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
# ℋᵢ = interior(ℋe, i₁:i₂, j₁:j₂, :)

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
# ℋᵢ = arch_array(arch, ℋᵢ)

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

#=
ocean = omip_ocean_component(grid)
sea_ice = omip_sea_ice_component(ocean.model)

coupled_model = OceanSeaIceModel(ocean, sea_ice)
stop_time = 4 * 360days
coupled_simulation = Simulation(coupled_model, Δt=10minutes, stop_time=stop_time)

set!(sea_ice.model, h=ℋᵢ)
#set!(ocean.model, T=Tᵢ, S=Sᵢ, e=1e-6)
set!(ocean.model, T=Tᵢ, S=Sᵢ)

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

    msg *= @sprintf(", extrema(T): (%.2f, %.2f) ᵒC", maximum(T), minimum(T))
    msg *= @sprintf(", extrema(S): (%.2f, %.2f) g kg⁻¹", minimum(S), maximum(S))

    # e = sim.model.ocean.model.tracers.e
    # msg *= @sprintf(", extrema(e): (%.2f, %.2f) m² s⁻²", minimum(e), maximum(e))

    @info msg
end

coupled_simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

N² = buoyancy_frequency(ocean.model)
κᶜ = ocean.model.diffusivity_fields.κᶜ
auxiliary_fields = (; N², κᶜ)
fields = merge(ocean.model.velocities, ocean.model.tracers, auxiliary_fields)

filename = "free_decay_heat_only_cold_ice_surface"

coupled_simulation.output_writers[:fields] = JLD2OutputWriter(ocean.model, fields;
                                                              filename = filename * "_fields",
                                                              indices = (:, :, Nz),
                                                              schedule = TimeInterval(10days),
                                                              overwrite_existing = true)

coupled_simulation.output_writers[:seaice] = JLD2OutputWriter(sea_ice.model, (; h = sea_ice.model.ice_thickness);
                                                              filename = filename * "_sea_ice_thickness",
                                                              schedule = TimeInterval(10days),
                                                              overwrite_existing = true)

run!(coupled_simulation)
=#
