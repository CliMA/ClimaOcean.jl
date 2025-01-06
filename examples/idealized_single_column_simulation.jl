using ClimaOcean
using ClimaSeaIce
using Oceananigans
using Oceananigans.Units
using SeawaterPolynomials

# Ocean state parameters
T₀ = 2    # Surface temperature, ᵒC
S₀ = 35   # Surface salinity
N² = 1e-4 # Buoyancy gradient due to temperature stratification
f = 0     # Coriolis parameter

# Atmospheric state parameters
Tₐ = 273.15 - 10 # Kelvin
u₁₀ = 10 # wind at 10 m, m/s
qₐ = 0.01 # specific humidity
Qℓ = 400 # shortwave radiation (W m⁻², positive means heating right now)

# Build the atmosphere
radiation = Radiation(ocean_albedo=0.1)
atmosphere_grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))
atmosphere_times = range(0, 1day, length=3)
atmosphere = PrescribedAtmosphere(atmosphere_grid, atmosphere_times)

parent(atmosphere.tracers.T) .= Tₐ     # K
parent(atmosphere.velocities.u) .= u₁₀ # m/s
parent(atmosphere.tracers.q) .= qₐ     # mass ratio
parent(atmosphere.downwelling_radiation.shortwave) .= Qℓ # W

# Build ocean model at rest with initial temperature stratification
grid = RectilinearGrid(size=100, z=(-400, 0), topology=(Flat, Flat, Bounded))
ocean = ocean_simulation(grid, coriolis=FPlane(; f))

eos = ocean.model.buoyancy.formulation.equation_of_state
g = ocean.model.buoyancy.formulation.gravitational_acceleration
α = SeawaterPolynomials.thermal_expansion(T₀, S₀, 0, eos)
dTdz = N² / (α * g)
Tᵢ(z) = T₀ + dTdz * z
set!(ocean.model, T=Tᵢ, S=S₀)

sea_ice_grid = RectilinearGrid(size=(),  topology=(Flat, Flat, Flat))
top_sea_ice_temperature = Field{Center, Center, Nothing}(sea_ice_grid)
top_heat_boundary_condition = PrescribedTemperature(top_sea_ice_temperature)
ice_thermodynamics = SlabSeaIceThermodynamics(sea_ice_grid; top_heat_boundary_condition)
                                              
top_sea_ice_heat_flux = Field{Center, Center, Nothing}(sea_ice_grid)
bottom_sea_ice_heat_flux = Field{Center, Center, Nothing}(sea_ice_grid)

sea_ice_model = SeaIceModel(sea_ice_grid;
                            top_heat_flux = top_sea_ice_heat_flux,
                            bottom_heat_flux = bottom_sea_ice_heat_flux,
                            ice_thermodynamics)

sea_ice = Simulation(sea_ice_model, Δt=10minutes)

model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
simulation = Simulation(model, Δt=10minutes, stop_time=10days)

Ql = model.fluxes.turbulent.fields.ocean.latent_heat
Qs = model.fluxes.turbulent.fields.ocean.sensible_heat
τx = model.fluxes.turbulent.fields.ocean.x_momentum
τy = model.fluxes.turbulent.fields.ocean.y_momentum
Fv = model.fluxes.turbulent.fields.ocean.water_vapor

# TODO: the total fluxes are defined on _interfaces_ between components:
# atmopshere_ocean, atmosphere_sea_ice, ocean_sea_ice. They aren't defined wrt to 
# just one component
Qo = model.fluxes.total.ocean.heat
Qi = model.fluxes.total.sea_ice.heat

fluxes = (; Qo, Qi, Ql, Qs, τx, τy, Fv)
ocean_outputs = merge(ocean.model.velocities, ocean.model.tracers)
sea_ice_outputs = (; h=sea_ice.model.ice_thickness)
outputs = merge(ocean_outputs, sea_ice_outputs, fluxes)

ow = JLD2OutputWriter(ocean.model, outputs,
                      schedule = TimeInterval(1hour),
                      filename = "idealized_atmosphere.jld2",
                      overwrite_existing = true)

simulation.output_writers[:ow] = ow

run!(simulation)

ut = FieldTimeSeries("idealized_atmosphere.jld2", "u")
vt = FieldTimeSeries("idealized_atmosphere.jld2", "v")
Tt = FieldTimeSeries("idealized_atmosphere.jld2", "T")
St = FieldTimeSeries("idealized_atmosphere.jld2", "S")
τx = FieldTimeSeries("idealized_atmosphere.jld2", "τx")
Qo = FieldTimeSeries("idealized_atmosphere.jld2", "Qo")
h = FieldTimeSeries("idealized_atmosphere.jld2", "h")
Nt = length(St)

using GLMakie

fig = Figure(size=(1000, 800))

axτ = Axis(fig[1, 1:3], xlabel="Time (days)", ylabel="Zonal momentum \n flux (N m⁻²)", xaxisposition=:top)
axQ = Axis(fig[2, 1:3], xlabel="Time (days)", ylabel="Heat flux \n (W m⁻²)")
axh = Axis(fig[3, 1:3], xlabel="Time (days)", ylabel="Sea ice thickness \n (m)")
axT = Axis(fig[4, 1], xlabel="Temperature (ᵒC)", ylabel="z (m)")
axS = Axis(fig[4, 2], xlabel="Salinity (g kg⁻¹)", ylabel="z (m)")
axu = Axis(fig[4, 3], xlabel="Velocities (m s⁻¹)", ylabel="z (m)", yaxisposition=:right)

slider = Slider(fig[5, 1:3], startvalue=1, range=1:Nt)
n = slider.value

un = @lift ut[$n]
vn = @lift vt[$n]
Tn = @lift Tt[$n]
Sn = @lift St[$n]

lines!(axT, Tn)
lines!(axS, Sn)
lines!(axu, un, label="u")
lines!(axu, vn, label="v")
axislegend(axu)

t = τx.times ./ days
tn = @lift t[$n]

lines!(axτ, t, τx[:], label="τx")
vlines!(axτ, tn)

lines!(axQ, t, Qo[:])
vlines!(axQ, tn)

lines!(axh, t, h[:])
vlines!(axh, tn)

rowsize!(fig.layout, 1, Relative(0.2))
rowsize!(fig.layout, 2, Relative(0.2))
rowsize!(fig.layout, 3, Relative(0.2))

xlims!(axS, 34.9, 35.1)
xlims!(axu, -0.1, 4.0)
hideydecorations!(axS, grid=false)
hidespines!(axT, :t, :r)
hidespines!(axS, :t, :l, :r)
hidespines!(axu, :t, :l)
hidespines!(axτ, :b, :r)
hidespines!(axQ, :t, :r)

record(fig, "idealized_atmosphere.mp4", 1:Nt, framerate=24) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end

