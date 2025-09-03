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
Nz = 40  # Vertical levels

r_faces = ExponentialCoordinate(Nz, -6000, 0)
grid    = TripolarGrid(CPU(); size=(Nx, Ny, Nz), z=r_faces, halo=(7, 7, 6))

# Regridding the bathymetry...
bottom_height = regrid_bathymetry(grid; major_basins=1, interpolation_passes=15)
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
                         closure = closures)

Oceananigans.set!(ocean.model, T=Metadatum(:temperature, dataset=ECCO4Monthly()), 
                               S=Metadatum(:salinity,    dataset=ECCO4Monthly()))

#####
##### The Atmosphere!!!
#####

spectral_grid = SpectralGrid(trunc=63, nlayers=8, Grid=FullClenshawGrid)
atmosphere = atmosphere_simulation(spectral_grid; output=true)
atmosphere.model.output.output_dt = Hour(3)

#####
##### The Sea-Ice!!!
#####

dynamics = ClimaOcean.SeaIceSimulations.sea_ice_dynamics(grid, ocean; solver=ClimaSeaIce.SeaIceDynamics.SplitExplicitSolver(100))
sea_ice = sea_ice_simulation(grid, ocean; dynamics, advection=WENO(order=7))

Oceananigans.set!(sea_ice.model, h=Metadatum(:sea_ice_thickness, dataset=ECCO4Monthly()), 
                                 ℵ=Metadatum(:sea_ice_concentration, dataset=ECCO4Monthly()))

#####
##### Coupled model
#####

Δt = convert(eltype(grid), atmosphere_model.time_stepping.Δt_sec)

# Remember in the future that reflected radiation is computed independently by speedy 
# so we need to communicate albedo in some way if this reflected radiation is to be
# absorbed by clouds. At the moment, speedy does not provide longwave down, until this
# is provided we just estimate a lower emissivity 
# *** THIS IS A HACK TO BE REMOVED ***
radiation = Radiation(ocean_emissivity=0.0, sea_ice_emissivity=0.0)
earth_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
earth = Oceananigans.Simulation(earth_model; Δt, stop_time=360days)

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

outputs = merge(ocean.model.velocities, ocean.model.tracers)
sea_ice_fields = merge(sea_ice.model.velocities, sea_ice.model.dynamics.auxiliaries.fields, 
                       (; h=sea_ice.model.ice_thickness, ℵ=sea_ice.model.ice_concentration))

ocean.output_writers[:free_surf] = JLD2Writer(ocean.model, (; η=ocean.model.free_surface.η);
                                              overwrite_existing=true,
                                              schedule=TimeInterval(3600 * 3),
                                              filename="ocean_free_surface.jld2")

ocean.output_writers[:surface] = JLD2Writer(ocean.model, outputs;
                                            overwrite_existing=true,
                                            schedule=TimeInterval(3600 * 3),
                                            filename="ocean_surface_fields.jld2",
                                            indices=(:, :, grid.Nz))

sea_ice.output_writers[:fields] = JLD2Writer(sea_ice.model, sea_ice_fields;
                                             overwrite_existing=true,
                                             schedule=TimeInterval(3600 * 3),
                                             filename="sea_ice_fields.jld2")

Qcao = earth.model.interfaces.atmosphere_ocean_interface.fluxes.sensible_heat
Qvao = earth.model.interfaces.atmosphere_ocean_interface.fluxes.latent_heat
τxao = earth.model.interfaces.atmosphere_ocean_interface.fluxes.x_momentum
τyao = earth.model.interfaces.atmosphere_ocean_interface.fluxes.y_momentum
Qcai = earth.model.interfaces.atmosphere_sea_ice_interface.fluxes.sensible_heat
Qvai = earth.model.interfaces.atmosphere_sea_ice_interface.fluxes.latent_heat
τxai = earth.model.interfaces.atmosphere_sea_ice_interface.fluxes.x_momentum
τyai = earth.model.interfaces.atmosphere_sea_ice_interface.fluxes.y_momentum
Qoi = earth.model.interfaces.net_fluxes.sea_ice_bottom.heat
Soi = earth.model.interfaces.sea_ice_ocean_interface.fluxes.salt
fluxes = (; Qcao, Qvao, τxao, τyao, Qcai, Qvai, τxai, τyai, Qoi, Soi)

earth.output_writers[:fluxes] = JLD2Writer(earth.model.ocean.model, fluxes;
                                           overwrite_existing=true,
                                           schedule=TimeInterval(3600 * 3),
                                           filename="intercomponent_fluxes.jld2")

# And add it as a callback to the simulation.
add_callback!(earth, progress, IterationInterval(100))

Oceananigans.run!(earth)
