using ClimaOcean
using ClimaOcean.ECCO: ECCO4Monthly
using OrthogonalSphericalShellGrids
using Oceananigans
using Oceananigans.Units
using CFTime
using Dates
using Printf

arch = CPU()
Nx = 120
Ny = 60
Nz = 20
z = exponential_z_faces(; Nz, depth=6000)
grid = TripolarGrid(arch; z, size=(Nx, Ny, Nz))

bottom_height = regrid_bathymetry(grid;
                                  minimum_depth = 10,
                                  interpolation_passes = 5,
                                  major_basins = 3)

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

# Closure
gm = Oceananigans.TurbulenceClosures.IsopycnalSkewSymmetricDiffusivity(κ_skew=1000, κ_symmetric=1000)
catke = Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivity()
closure = (gm, catke)

# Polar restoring
@inline function restoring_mask(λ, φ, z, t=0)
    ϵN = (φ - 75) / 5
    ϵN = clamp(ϵN, zero(ϵN), one(ϵN))
    ϵS = - (φ + 75) / 5
    ϵS = clamp(ϵS, zero(ϵS), one(ϵS))
    return ϵN + ϵS
end

restoring_rate = 1 / 1days

restoring_mask_field = CenterField(grid)
set!(restoring_mask_field, restoring_mask)

@inline sponge_layer(λ, φ, z, t, c, ω) = - restoring_mask(λ, φ, z, t) * ω * c
Fu = Forcing(sponge_layer, field_dependencies=:u, parameters=restoring_rate)
Fv = Forcing(sponge_layer, field_dependencies=:v, parameters=restoring_rate)

dates = DateTimeProlepticGregorian(1993, 11, 1) : Month(1) : DateTimeProlepticGregorian(1994, 11, 1)
temperature = ECCOMetadata(:temperature, dates, ECCO4Monthly())
salinity = ECCOMetadata(:salinity, dates, ECCO4Monthly())

FT = ECCORestoring(arch, temperature; grid, mask=restoring_mask_field, rate=restoring_rate)
FS = ECCORestoring(arch, salinity;    grid, mask=restoring_mask_field, rate=restoring_rate)
forcing = (T=FT, S=FS, u=Fu, v=Fv)

momentum_advection = VectorInvariant()
tracer_advection = Centered(order=2)
ocean = ocean_simulation(grid; momentum_advection, tracer_advection,
                         closure, forcing,
                         tracers = (:T, :S, :e))

set!(ocean.model,
     T = ECCOMetadata(:temperature; dates=first(dates)),
     S = ECCOMetadata(:salinity; dates=first(dates)))

radiation = Radiation(arch)
atmosphere = JRA55_prescribed_atmosphere(arch; backend=JRA55NetCDFBackend(41))
sea_ice = ClimaOcean.OceanSeaIceModels.MinimumTemperatureSeaIce()
coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

simulation = Simulation(coupled_model; Δt=15minutes, stop_time=2*365days)

wall_time = Ref(time_ns())

function progress(sim)
    ocean = sim.model.ocean
    u, v, w = ocean.model.velocities
    T = ocean.model.tracers.T
    Tmax = maximum(interior(T))
    Tmin = minimum(interior(T))
    umax = (maximum(abs, interior(u)),
            maximum(abs, interior(v)),
            maximum(abs, interior(w)))

    step_time = 1e-9 * (time_ns() - wall_time[])

    @info @sprintf("Time: %s, n: %d, Δt: %s, max|u|: (%.2e, %.2e, %.2e) m s⁻¹, extrema(T): (%.2f, %.2f) ᵒC, wall time: %s \n",
                   prettytime(sim), iteration(sim), prettytime(sim.Δt),
                   umax..., Tmax, Tmin, prettytime(step_time))

     wall_time[] = time_ns()

     return nothing
end

add_callback!(simulation, progress, IterationInterval(10))

run!(simulation)
