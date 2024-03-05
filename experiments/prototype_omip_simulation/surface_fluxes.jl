using Oceananigans
using Oceananigans.Units
using Oceananigans.BuoyancyModels: buoyancy_frequency
using Oceananigans.Units: Time

using ClimaOcean
using ClimaOcean.OceanSeaIceModels: Radiation
using ClimaOcean.DataWrangling.JRA55: JRA55_prescribed_atmosphere
using ClimaOcean.DataWrangling.ECCO2: ecco2_field

using GLMakie
using Printf
using Dates

start_time = time_ns()

include("omip_components.jl")

arch = CPU()
epoch = Date(1992, 1, 1)
date = Date(1992, 10, 1)
start_seconds = Second(date - epoch).value
ue = ecco2_field(:u_velocity, date)
ve = ecco2_field(:v_velocity, date)
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

latitude = (-75, +65)
longitude = (0, 360)

i₁ = 4 * first(longitude) + 1
i₂ = 1440 - 4 * (360 - last(longitude))
Nx = i₂ - i₁ + 1

j₁ = 4 * (90 + first(latitude)) + 1
j₂ = 720 - 4 * (90 - last(latitude))
Ny = j₂ - j₁ + 1

Nz = size(Te, 3)

Tᵢ = Te[i₁:i₂, j₁:j₂,   Nz:Nz]
Sᵢ = Se[i₁:i₂, j₁:j₂,   Nz:Nz]
uᵢ = ue[i₁:i₂, j₁:j₂,   Nz:Nz]
vᵢ = ve[i₁:i₂, j₁:j₂+1, Nz:Nz]

# Construct bottom_height depth by analyzing T
Nx, Ny, Nz = size(Tᵢ)
bottom_height = zeros(Nx, Ny)

for i = 1:Nx, j = 1:Ny
    @inbounds bottom_height[i, j] = ifelse(isnan(Tᵢ[i, j, 1]), Inf, -Inf)
end

grid = LatitudeLongitudeGrid(arch; latitude, longitude,
                             size = (Nx, Ny, 1),
                             halo = (7, 7, 7),
                             z = (-10, 0),
                             topology = (Periodic, Bounded, Bounded))

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

elapsed = time_ns() - start_time
@info "Grid constructed. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

ocean = omip_ocean_component(grid)
elapsed = time_ns() - start_time
@info "Ocean component built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

atmosphere = JRA55_prescribed_atmosphere(1:2, backend=InMemory(), architecture=arch)
elapsed = time_ns() - start_time
@info "Atmosphere built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

#=
fig = Figure()
axu = Axis(fig[1, 1])
axv = Axis(fig[1, 2])
axT = Axis(fig[1, 3])
axq = Axis(fig[1, 4])

ua = atmosphere.velocities.u
va = atmosphere.velocities.v
Ta = atmosphere.tracers.T
qa = atmosphere.tracers.q

heatmap!(axu, interior(ua, :, :, 1, 1))
heatmap!(axv, interior(va, :, :, 1, 1))
heatmap!(axT, interior(Ta, :, :, 1, 1))
heatmap!(axq, interior(qa, :, :, 1, 1))

display(fig)

ocean.model.clock.time = start_seconds
ocean.model.clock.iteration = 0
set!(ocean.model, T=Tᵢ, S=Sᵢ, e=1e-6)

ua = atmosphere.velocities.u
va = atmosphere.velocities.v
Ta = atmosphere.tracers.T
qa = atmosphere.tracers.q
times = ua.times
=#

sea_ice = nothing
radiation = Radiation()
coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
elapsed = time_ns() - start_time
@info "Coupled model built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

# Build flux outputs
Jᵘ = coupled_model.fluxes.total.ocean.momentum.u
Jᵛ = coupled_model.fluxes.total.ocean.momentum.v
Jᵀ = coupled_model.fluxes.total.ocean.tracers.T
Jˢ = coupled_model.fluxes.total.ocean.tracers.S

E  = coupled_model.fluxes.turbulent.fields.water_vapor
Fʳ = atmosphere.freshwater_flux.rain
Fˢ = atmosphere.freshwater_flux.snow

Qᶜ = coupled_model.fluxes.turbulent.fields.sensible_heat
Qᵉ = coupled_model.fluxes.turbulent.fields.latent_heat
Qˡ = atmosphere.downwelling_radiation.longwave
Qˢ = atmosphere.downwelling_radiation.shortwave

ρₒ = coupled_model.fluxes.ocean_reference_density
cₚ = coupled_model.fluxes.ocean_heat_capacity

P = Field(- Fʳ[1] - Fˢ[1])
ΣQ = Field(ρₒ * cₚ * Jᵀ)
u★ = Field(sqrt(sqrt(Jᵘ^2 + Jᵛ^2)))

N² = buoyancy_frequency(ocean.model)
κᶜ = ocean.model.diffusivity_fields.κᶜ
To = ocean.model.tracers.T
qa = atmosphere.tracers.q

compute!(ΣQ)
compute!(P)
compute!(u★)

fig = Figure(size=(1200, 1800))

axτ = Axis(fig[1, 1], title="u★")
axT = Axis(fig[2, 1], title="T")
axq = Axis(fig[3, 1], title="q")
axE = Axis(fig[4, 1], title="E")
axP = Axis(fig[5, 1], title="P")

axQt = Axis(fig[1, 2], title="Net heat flux")
axQl = Axis(fig[2, 2], title="Incoming longwave heat flux")
axQs = Axis(fig[3, 2], title="Shortwave / solar heat flux")
axQe = Axis(fig[4, 2], title="Evaporative heat flux")
axQc = Axis(fig[5, 2], title="Conductive / sensible heat flux")

u★i = interior(u★, :, :, 1) 
Toi = interior(To, :, :, 1)
qai = interior(qa, :, :, 1, 1)
Ei  = interior(E,  :, :, 1)
Pi  = interior(P,  :, :, 1)

Flim = 1e-6
hmτ = heatmap!(axτ, u★i, colormap=:solar, colorrange=(0, 1e-1))
hmT = heatmap!(axT, Toi, colormap=:thermal, colorrange=(0, 10))
hmq = heatmap!(axq, qai, colormap=:grays)
hmE = heatmap!(axE, Ei, colormap=:balance, colorrange=(-Flim, Flim))
hmP = heatmap!(axP, Pi, colormap=:balance, colorrange=(-Flim, Flim))

Colorbar(fig[1, 0], hmτ, label="Friction velocity (m s⁻¹)")
Colorbar(fig[2, 0], hmT, label="Ocean surface temperature (ᵒC)")
Colorbar(fig[3, 0], hmq, label="Atmosphere specific humidity")
Colorbar(fig[4, 0], hmE, label="Evaporation freshwater flux (m s⁻¹)")
Colorbar(fig[5, 0], hmP, label="Precipitation freshwater flux (m s⁻¹)")

ΣQi = interior(ΣQ, :, :, 1)
Qˢi = - interior(Qˢ, :, :, 1, 1)
Qˡi = - interior(Qˡ, :, :, 1, 1)
Qᵉi = interior(Qᵉ, :, :, 1)
Qᶜi = interior(Qᶜ, :, :, 1)

Qlim = 1000
hmt = heatmap!(axQt, ΣQi, colormap=:balance, colorrange=(-Qlim, Qlim)) 
hms = heatmap!(axQs, Qˢi, colormap=:balance, colorrange=(-Qlim, Qlim)) 
hml = heatmap!(axQl, Qˡi, colormap=:balance, colorrange=(-Qlim, Qlim)) 
hme = heatmap!(axQe, Qᵉi, colormap=:balance, colorrange=(-Qlim, Qlim)) 
hmc = heatmap!(axQc, Qᶜi, colormap=:balance, colorrange=(-Qlim, Qlim)) 

Colorbar(fig[1, 3], hmt, label="Net heat flux")
Colorbar(fig[2, 3], hml, label="Incoming longwave heat flux")
Colorbar(fig[3, 3], hms, label="Shortwave / solar heat flux")
Colorbar(fig[4, 3], hme, label="Evaporative heat flux")
Colorbar(fig[5, 3], hmc, label="Conductive / sensible heat flux")

display(fig)

