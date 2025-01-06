using ClimaOcean
using ClimaOcean.ECCO: ECCO4Monthly
using ClimaSeaIce
using OrthogonalSphericalShellGrids
using Oceananigans
using Oceananigans.Units
using CFTime
using Dates
using Printf

arch = CPU()

#=
Nx = 120
Ny = 60
Nz = 50
z_faces = exponential_z_faces(; Nz, depth=6000, h=34)
underlying_grid = TripolarGrid(arch; size=(Nx, Ny, Nz), z=z_faces)
=#

latitude = (-80, -20)
longitude = (0, 360)

Nx = 120
Ny = 20
Nz = 30

z_faces = exponential_z_faces(; Nz, depth=6000, h=34)
underlying_grid = LatitudeLongitudeGrid(arch; size=(Nx, Ny, Nz), latitude, longitude, z=z_faces)
bottom_height = regrid_bathymetry(underlying_grid)
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_height))

gm = Oceananigans.TurbulenceClosures.IsopycnalSkewSymmetricDiffusivity(κ_skew=4000, κ_symmetric=4000)
catke = ClimaOcean.OceanSimulations.default_ocean_closure()
viscous_closure = Oceananigans.TurbulenceClosures.HorizontalScalarDiffusivity(ν=4000)

dates = DateTimeProlepticGregorian(1993, 1, 1) : Month(1) : DateTimeProlepticGregorian(1993, 12, 1)
temperature = ECCOMetadata(:temperature, dates, ECCO4Monthly())
salinity = ECCOMetadata(:salinity, dates, ECCO4Monthly())

restoring_rate  = 1/2days
mask = LinearlyTaperedPolarMask(southern=(-80, -70), northern=(70, 90))
FT = ECCORestoring(temperature, grid; mask, rate=restoring_rate)
FS = ECCORestoring(salinity, grid; mask, rate=restoring_rate)

ocean = ocean_simulation(grid;
                         momentum_advection =  VectorInvariant(),
                         tracer_advection = Centered(order=2),
                         closure = (gm, catke, viscous_closure),
                         # forcing = (T=FT, S=FT),
                         tracers = (:T, :S, :e))

set!(ocean.model, T=ECCOMetadata(:temperature; dates=first(dates)),
                  S=ECCOMetadata(:salinity; dates=first(dates)))

radiation  = Radiation(arch)
atmosphere = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(20))

sea_ice_grid = LatitudeLongitudeGrid(arch; size=(Nx, Ny), latitude, longitude, topology=(Periodic, Bounded, Flat))
land_mask = bottom_height .>= 0
sea_ice_grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBoundary(land_mask))
top_sea_ice_temperature = Field{Center, Center, Nothing}(sea_ice_grid)
top_heat_boundary_condition = PrescribedTemperature(top_sea_ice_temperature)
ice_thermodynamics = SlabSeaIceThermodynamics(sea_ice_grid; top_heat_boundary_condition)
                                              
top_sea_ice_heat_flux = Field{Center, Center, Nothing}(sea_ice_grid)
bottom_sea_ice_heat_flux = Field{Center, Center, Nothing}(sea_ice_grid)

sea_ice_model = SeaIceModel(sea_ice_grid;
                            top_heat_flux = top_sea_ice_heat_flux,
                            bottom_heat_flux = bottom_sea_ice_heat_flux,
                            ice_thermodynamics)

sea_ice = Simulation(sea_ice_model, Δt=20minutes)
set!(sea_ice_model, h=1, ℵ=1)

coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation) 
simulation = Simulation(coupled_model; Δt=20minutes, stop_time=3days)

wall_time = Ref(time_ns())

function progress(sim)
    sea_ice = sim.model.sea_ice
    h = sea_ice.model.ice_thickness
    T = top_sea_ice_temperature
    hmax = maximum(interior(h))
    Timax = maximum(interior(Ti))
    Timin = minimum(interior(Ti))

    ocean = sim.model.ocean
    u, v, w = ocean.model.velocities
    T = ocean.model.tracers.T
    Tmax = maximum(interior(T))
    Tmin = minimum(interior(T))
    umax = (maximum(abs, interior(u)), maximum(abs, interior(v)), maximum(abs, interior(w)))
    step_time = 1e-9 * (time_ns() - wall_time[])

    @info @sprintf("Time: %s, n: %d, Δt: %s, max(h): %.1f, extrema(Ti): (%.2f, %.2f) K, \
                   max|u|: (%.2e, %.2e, %.2e) m s⁻¹, \
                   extrema(T): (%.2f, %.2f) ᵒC, wall time: %s \n",
                   prettytime(sim), iteration(sim), prettytime(sim.Δt),
                   hmax, Timin, Timax, umax..., Tmin, Tmax, prettytime(step_time))

    wall_time[] = time_ns()

    return nothing
end

add_callback!(simulation, progress, IterationInterval(10))

Ql = coupled_model.fluxes.turbulent.fields.ocean.latent_heat
Qs = coupled_model.fluxes.turbulent.fields.ocean.sensible_heat
τx = coupled_model.fluxes.turbulent.fields.ocean.x_momentum
τy = coupled_model.fluxes.turbulent.fields.ocean.y_momentum
Fv = coupled_model.fluxes.turbulent.fields.ocean.water_vapor

# TODO: the total fluxes are defined on _interfaces_ between components:
# atmopshere_ocean, atmosphere_sea_ice, ocean_sea_ice. They aren't defined wrt to 
# just one component
Qo = coupled_model.fluxes.total.ocean.heat
Qi = coupled_model.fluxes.total.sea_ice.top_heat

fluxes = (; Qo, Qi, Ql, Qs, τx, τy, Fv)
ocean_outputs = merge(ocean.model.velocities, ocean.model.tracers)

h = sea_ice_model.ice_thickness
ℵ = sea_ice_model.ice_concentration
Ti = top_sea_ice_temperature
sea_ice_outputs = (; h, ℵ, Ti)

outputs = merge(ocean_outputs, sea_ice_outputs, fluxes)

ow = JLD2OutputWriter(ocean.model, outputs,
                      filename = "three_degree_simulation.jld2",
                      schedule = TimeInterval(1hours),
                      overwrite_existing = true)

simulation.output_writers[:jld2] = ow

run!(simulation)

using GLMakie

ht = FieldTimeSeries("three_degree_simulation.jld2", "h")
ℵt = FieldTimeSeries("three_degree_simulation.jld2", "ℵ")
Nt = length(ht)

fig = Figure()
axh = Axis(fig[1, 1])
axℵ = Axis(fig[2, 1])
slider = Slider(fig[3, 1], range=1:Nt, startvalue=1)
n = slider.value

hn = @lift ht[$n]
ℵn = @lift ℵt[$n]
heatmap!(axh, hn)
heatmap!(axℵ, ℵn)

display(fig)
