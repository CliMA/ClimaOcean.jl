using ClimaOcean
using ClimaSeaIce
using Oceananigans
using Oceananigans.Grids
using Oceananigans.Units
using Oceananigans.OrthogonalSphericalShellGrids
using ClimaOcean.OceanSimulations
using ClimaOcean.ECCO
using ClimaOcean.ECCO: all_ECCO_dates
using Printf

r_faces = ClimaOcean.exponential_z_faces(; Nz=30, h=10, depth=3000)
z_faces = MutableVerticalDiscretization(r_faces)

Nx = 180 # longitudinal direction -> 250 points is about 1.5ᵒ resolution
Ny = 180 # meridional direction -> same thing, 48 points is about 1.5ᵒ resolution
Nz = length(r_faces) - 1

grid = RotatedLatitudeLongitudeGrid(size = (Nx, Ny, Nz), 
                                    latitude = (-45, 45),
                                    longitude = (-45, 45),
                                    z = z_faces,
                                    north_pole = (180, 0),
                                    halo = (5, 5, 4),
                                    topology = (Bounded, Bounded, Bounded))

bottom_height = regrid_bathymetry(grid; minimum_depth=15, major_basins=1)

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

#####
##### A Prescribed Ocean model
#####

# ...with prescribed velocity and tracer fields
version = ECCO2Daily()
dates   = all_ECCO_dates(version)[1:200]

u = ECCOFieldTimeSeries(:u_velocity,  version; dates, time_indices_in_memory=2)
v = ECCOFieldTimeSeries(:v_velocity,  version; dates, time_indices_in_memory=2)
T = ECCOFieldTimeSeries(:temperature, version; dates, time_indices_in_memory=2)
S = ECCOFieldTimeSeries(:salinity,    version; dates, time_indices_in_memory=2)

ocean_model = PrescribedOcean((; u, v, T, S); grid)
ocean       = Simulation(ocean_model, Δt=10minutes, verbose=false)

#####
##### A Prognostic Sea-ice model
#####

sea_ice = sea_ice_simulation(grid; dynamics=nothing, advection=nothing) 

set!(sea_ice.model, h=ECCOMetadata(:sea_ice_thickness), 
                    ℵ=ECCOMetadata(:sea_ice_concentration))

#####
##### Atmosphere model
#####

atmosphere = JRA55PrescribedAtmosphere(; backend=JRA55NetCDFBackend(40))
radiation  = Radiation()

#####
##### Arctic coupled model
#####

arctic = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
arctic = Simulation(arctic, Δt=600, stop_time=365days)

h = sea_ice.model.ice_thickness
ℵ = sea_ice.model.ice_concentration
Gh = sea_ice.model.timestepper.Gⁿ.h
Gℵ = sea_ice.model.timestepper.Gⁿ.ℵ
Tu = arctic.model.interfaces.atmosphere_sea_ice_interface.temperature

sea_ice.output_writers[:vars] = JLD2OutputWriter(sea_ice.model, (; h, ℵ, Gh, Gℵ, Tu),
                                                 filename = "sea_ice_quantities.jld2",
                                                 schedule = IterationInterval(12))

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