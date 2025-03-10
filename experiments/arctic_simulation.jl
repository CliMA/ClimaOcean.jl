using ClimaOcean
using ClimaSeaIce
using Oceananigans
using Oceananigans.Grids
using Oceananigans.Units
using Oceananigans.OrthogonalSphericalShellGrids
using ClimaOcean.OceanSimulations
using ClimaOcean.ECCO
using ClimaOcean.ECCO: all_ECCO_dates
using ClimaSeaIce.SeaIceThermodynamics: IceWaterThermalEquilibrium
using Printf

using CUDA
CUDA.device!(1)

r_faces = ClimaOcean.exponential_z_faces(; Nz=30, h=10, depth=2000)
z_faces = MutableVerticalDiscretization(r_faces)

Nx = 180 # longitudinal direction -> 250 points is about 1.5ᵒ resolution
Ny = 180 # meridional direction -> same thing, 48 points is about 1.5ᵒ resolution
Nz = length(r_faces) - 1

grid = RotatedLatitudeLongitudeGrid(GPU(), size = (Nx, Ny, Nz), 
                                           latitude = (-45, 45),
                                           longitude = (-45, 45),
                                           z = r_faces,
                                           north_pole = (180, 0),
                                           halo = (5, 5, 4),
                                           topology = (Bounded, Bounded, Bounded))

bottom_height = regrid_bathymetry(grid; minimum_depth=15, major_basins=1)

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

#####
##### A Propgnostic Ocean model
#####

# A very diffusive ocean
momentum_advection = WENOVectorInvariant(order=3) 
tracer_advection   = WENO(order=3)

free_surface = SplitExplicitFreeSurface(grid; cfl=0.7) 
closure = ClimaOcean.OceanSimulations.default_ocean_closure()

ocean = ocean_simulation(grid; 
                         momentum_advection, 
                         tracer_advection, 
                         free_surface,
                         closure)

set!(ocean.model, T=ECCOMetadata(:temperature),
                  S=ECCOMetadata(:salinity))

#####
##### A Prognostic Sea-ice model
#####

# Remember to pass the SSS as a bottom bc to the sea ice!
SSS = view(ocean.model.tracers.S, :, :, grid.Nz)
bottom_heat_boundary_condition = IceWaterThermalEquilibrium(SSS)

sea_ice = sea_ice_simulation(grid; bottom_heat_boundary_condition, dynamics=nothing, advection=nothing) 

set!(sea_ice.model, h=ECCOMetadata(:sea_ice_thickness), 
                    ℵ=ECCOMetadata(:sea_ice_concentration))

##### 
##### A Prescribed Atmosphere model
#####

atmosphere = JRA55PrescribedAtmosphere(GPU(); backend=JRA55NetCDFBackend(40))
radiation  = Radiation()

#####
##### Arctic coupled model
#####

arctic = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
arctic = Simulation(arctic, Δt=5minutes, stop_time=365days)

# Sea-ice variables
h  = sea_ice.model.ice_thickness
ℵ  = sea_ice.model.ice_concentration
Gh = sea_ice.model.timestepper.Gⁿ.h
Gℵ = sea_ice.model.timestepper.Gⁿ.ℵ

# Fluxes
Tu = arctic.model.interfaces.atmosphere_sea_ice_interface.temperature
Qˡ = arctic.model.interfaces.atmosphere_sea_ice_interface.fluxes.latent_heat
Qˢ = arctic.model.interfaces.atmosphere_sea_ice_interface.fluxes.sensible_heat
Qⁱ = arctic.model.interfaces.sea_ice_ocean_interface.fluxes.interface_heat
Qᶠ = arctic.model.interfaces.sea_ice_ocean_interface.fluxes.frazil_heat
Qᵗ = arctic.model.interfaces.net_fluxes.sea_ice_top.heat
Qᴮ = arctic.model.interfaces.net_fluxes.sea_ice_bottom.heat

# Output writers
arctic.output_writers[:vars] = JLD2OutputWriter(sea_ice.model, (; h, ℵ, Gh, Gℵ, Tu, Qˡ, Qˢ, Qⁱ, Qᶠ, Qᵗ, Qᴮ),
                                                 filename = "sea_ice_quantities.jld2",
                                                 schedule = IterationInterval(12),
                                                 overwrite_existing=true)

arctic.output_writers[:avrages] = JLD2OutputWriter(sea_ice.model, (; h, ℵ, Tu, Qˡ, Qˢ, Qⁱ, Qᶠ, Qᵗ, Qᴮ),
                                                    filename = "averaged_sea_ice_quantities.jld2",
                                                    schedule = AveragedTimeInterval(1days),
                                                    overwrite_existing=true)

wall_time = Ref(time_ns())

using Statistics

function progress(sim)
    sea_ice = sim.model.sea_ice
    hmax = maximum(sea_ice.model.ice_thickness)
    ℵmax = maximum(sea_ice.model.ice_concentration)
    hmean = mean(sea_ice.model.ice_thickness)
    ℵmean = mean(sea_ice.model.ice_concentration)
    Tmax = maximum(sim.model.interfaces.atmosphere_sea_ice_interface.temperature)
    Tmin = minimum(sim.model.interfaces.atmosphere_sea_ice_interface.temperature)

    step_time = 1e-9 * (time_ns() - wall_time[])

    msg1 = @sprintf("time: %s, iteration: %d, Δt: %s, ", prettytime(sim), iteration(sim), prettytime(sim.Δt))
    msg2 = @sprintf("max(h): %.2e m, max(ℵ): %.2e ", hmax, ℵmax)
    msg3 = @sprintf("mean(h): %.2e m, mean(ℵ): %.2e ", hmean, ℵmean)
    msg4 = @sprintf("extrema(T): (%.2f, %.2f) ᵒC, ", Tmax, Tmin)
    msg5 = @sprintf("wall time: %s \n", prettytime(step_time))

    @info msg1 * msg2 * msg3 * msg4 * msg5

     wall_time[] = time_ns()

     return nothing
end

# And add it as a callback to the simulation.
add_callback!(arctic, progress, IterationInterval(10))

run!(arctic)

#####
##### Comparison to ECCO Climatology
#####

version = ECCO4Monthly()
dates   = all_ECCO_dates(version)[1:12]

h_metadata = ECCOMetadata(:sea_ice_thickness;     version, dates)
ℵ_metadata = ECCOMetadata(:sea_ice_concentration; version, dates)

# Montly averaged ECCO data
hE = ECCOFieldTimeSeries(h_metadata, grid; time_indices_in_memory=12)
ℵE = ECCOFieldTimeSeries(ℵ_metadata, grid; time_indices_in_memory=12)

# Daily averaged Model output
h = FieldTimeSeries("averaged_sea_ice_quantities.jld2", "h")
ℵ = FieldTimeSeries("averaged_sea_ice_quantities.jld2", "ℵ")

# Montly average the model output
hm = FieldTimeSeries{Center, Center, Nothing}(grid, hE.times; backend=InMemory())
ℵm = FieldTimeSeries{Center, Center, Nothing}(grid, hE.times; backend=InMemory())

for (i, time) in enumerate(hm.times)
    counter = 0
    for (j, t) in enumerate(h.times)
        if t ≤ time
            hm[i] .+= h[j]
            ℵm[i] .+= ℵ[j]
            counter += 1
        end
    end
    @show counter
    hm[i] ./= counter
    ℵm[i] ./= counter
end
