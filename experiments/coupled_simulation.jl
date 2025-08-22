#####
##### A one degree coupled Atmosphere - Ocean - Sea Ice model
#####

using Oceananigans, ClimaSeaIce, SpeedyWeather, ClimaOcean
using Oceananigans, Oceananigans.Units
using Printf, Statistics, Dates

#####
##### The Ocean!!!
#####

Nx = 360 # Longitudinal direction 
Ny = 180 # Meridional direction 
Nz = 60  # Vertical levels

r_faces = ExponentialCoordinate(Nz, -6000, 0)
grid    = TripolarGrid(CPU(); size=(Nx, Ny, Nz), z=r_faces, halo=(5, 5, 4))

# Regridding the bathymetry...
bottom_height = regrid_bathymetry(grid; minimum_depth=15, major_basins=1, interpolation_passes=15)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)

# Advection
momentum_advection = WENOVectorInvariant(order=5)
tracer_advection   = WENO(order=5)

# Free Surface
free_surface = SplitExplicitFreeSurface(grid; substeps=80)

# Parameterizations
catke_closure = RiBasedVerticalDiffusivity()
eddy_closure  = Oceananigans.TurbulenceClosures.IsopycnalSkewSymmetricDiffusivity(κ_skew=1e3, κ_symmetric=1e3)
closures      = (catke_closure, eddy_closure, VerticalScalarDiffusivity(κ=1e-5, ν=1e-4))

# The ocean simulation
ocean = ocean_simulation(grid; 
                         momentum_advection,
                         tracer_advection,
                         free_surface,
                         timestepper = :SplitRungeKutta3,
                         closure = closures,
                         radiative_forcing = nothing)

Oceananigans.set!(ocean.model, T=Metadatum(:temperature, dataset=ECCO4Monthly()), 
                               S=Metadatum(:salinity,    dataset=ECCO4Monthly()))

#####
##### The Atmosphere!!!
#####

spectral_grid = SpectralGrid(trunc=31, nlayers=8, Grid=FullClenshawGrid)

humidity_flux_ocean = PrescribedOceanHumidityFlux(spectral_grid)
humidity_flux_land = SurfaceLandHumidityFlux(spectral_grid)
surface_humidity_flux = SurfaceHumidityFlux(ocean=humidity_flux_ocean, land=humidity_flux_land)

ocean_heat_flux = PrescribedOceanHeatFlux(spectral_grid)
land_heat_flux = SurfaceLandHeatFlux(spectral_grid)
surface_heat_flux = SurfaceHeatFlux(ocean=ocean_heat_flux, land=land_heat_flux)

atmosphere_model = PrimitiveWetModel(spectral_grid;
                                     surface_heat_flux,
                                     surface_humidity_flux,
                                     sea_ice=NoSeaIce()) # This is provided by ClimaSeaIce

atmosphere = initialize!(atmosphere_model)
initialize!(atmosphere)

function initialize_atmospheric_state!(simulation::SpeedyWeather.Simulation)
    progn, diagn, model  = SpeedyWeather.unpack(simulation)

    (; time) = progn.clock                           # current time

    # set the tendencies back to zero for accumulation
    fill!(diagn.tendencies, 0, typeof(model))

    if model.physics                   
        # calculate all parameterizations
        SpeedyWeather.parameterization_tendencies!(diagn, progn, time, model)
    end
    
    return nothing
end

initialize_atmospheric_state!(atmosphere)
# atmosphere.model.feedback.verbose = false

#####
##### The Sea-Ice!!!
#####

dynamics = ClimaOcean.SeaIceSimulations.sea_ice_dynamics(grid, ocean)
sea_ice = sea_ice_simulation(grid, ocean; dynamics, advection=WENO(order=7))

Oceananigans.set!(sea_ice.model, h=Metadatum(:sea_ice_thickness, dataset=ECCO4Monthly()), 
                                 ℵ=Metadatum(:sea_ice_concentration, dataset=ECCO4Monthly()))

#####
##### Coupled model
#####

Δt = convert(eltype(grid), atmosphere_model.time_stepping.Δt_sec)

# Remember in the future that reflected radiation is computed independently by speedy 
# so we need to communicate albedo in some way if this reflected radiation is to be
# absorbed by clouds 
radiation = Radiation(ocean_emissivity=0, sea_ice_emissivity=0)
earth_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
earth = Oceananigans.Simulation(earth_model; Δt, stop_time=20days)

wall_time = Ref(time_ns())

function progress(sim)
    sea_ice = sim.model.sea_ice
    ocean   = sim.model.ocean
    hmax  = maximum(sea_ice.model.ice_thickness)
    ℵmax  = maximum(sea_ice.model.ice_concentration)
    uimax = maximum(abs, sea_ice.model.velocities.u)
    vimax = maximum(abs, sea_ice.model.velocities.v)
    uomax = maximum(abs, ocean.model.velocities.u)
    vomax = maximum(abs, ocean.model.velocities.v)

    step_time = 1e-9 * (time_ns() - wall_time[])

    msg1 = @sprintf("time: %s, iteration: %d, Δt: %s, ", prettytime(sim), iteration(sim), prettytime(sim.Δt))
    msg2 = @sprintf("max(h): %.2e m, max(ℵ): %.2e ", hmax, ℵmax)
    msg3 = @sprintf("max uᵢ: (%.2f, %.2f) m s⁻¹, ", uimax, vimax)
    msg4 = @sprintf("max uₒ: (%.2f, %.2f) m s⁻¹, ", uomax, vomax)
    msg5 = @sprintf("wall time: %s \n", prettytime(step_time))

    @info msg1 * msg2 * msg3 * msg4 * msg5

    wall_time[] = time_ns()

     return nothing
end

# And add it as a callback to the simulation.
# add_callback!(earth, progress, IterationInterval(10))

# Oceananigans.run!(earth)
time_step!(earth)