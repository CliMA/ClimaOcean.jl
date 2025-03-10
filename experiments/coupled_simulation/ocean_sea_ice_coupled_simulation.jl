using ClimaOcean
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids
using Oceananiagans.OrthogonalSphericalShellGrids
using CFTime
using Dates
using Printf

using Oceananigans.Utils: launch!
using ClimaOcean.DataWrangling: NearestNeighborInpainting
using ClimaSeaIce.SeaIceThermodynamics: melting_temperature
using Oceananigans.Grids: architecture
using KernelAbstractions: @kernel, @index


arch = CPU()

depth = 1000meters
Nz    = 10
h     = 3 

r_faces = ClimaOcean.exponential_z_faces(; Nz, h, depth)
z_faces = MutableVerticalDiscretization(r_faces)

Nx = 256 # longitudinal direction -> 250 points is about 1.5ᵒ resolution
Ny = 128 # meridional direction -> same thing, 48 points is about 1.5ᵒ resolution
Nz   = length(r_faces) - 1
grid = TripolarGrid(arch, Float64; size=(Nx, Ny, Nz), z=z_faces)
sea_ice_grid = TripolarGrid(arch, Float64; size=(Nx, Ny, 1), z = (-10, 0))

# ## Adding a bathymetry to the grid
url = "https://www.dropbox.com/scl/fi/zy1cu64ybn93l67rjgiq0/Downsampled_ETOPO_2022.nc?rlkey=5upqvoxrnljj205amqf663vcw&st=ou8b32tt&dl=0"
filename = isfile("Downsampled_ETOPO_2022.nc") ? "Downsampled_ETOPO_2022.nc" : download(url, "Downsampled_ETOPO_2022.nc")
bottom_height = regrid_bathymetry(grid; minimum_depth=15, major_basins=1, filename, dir="./")

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)
sea_ice_grid = ImmersedBoundaryGrid(sea_ice_grid, GridFittedBottom(bottom_height))

#####
##### Ocean model
#####

momentum_advection = WENOVectorInvariant(order=3) 
tracer_advection   = Centered()

free_surface = SplitExplicitFreeSurface(grid; substeps=30) 

using Oceananigans.TurbulenceClosures: IsopycnalSkewSymmetricDiffusivity, 
                                       DiffusiveFormulation

using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEVerticalDiffusivity

eddy_closure = IsopycnalSkewSymmetricDiffusivity(κ_skew=1e3, κ_symmetric=1e3, skew_flux_formulation=DiffusiveFormulation())
vertical_mixing = ClimaOcean.OceanSimulations.default_ocean_closure() 

closure = (eddy_closure, vertical_mixing) 

ocean = ocean_simulation(grid; 
                         momentum_advection, 
                         tracer_advection, 
                         closure, 
                         free_surface)

#####
##### Sea-ice model
#####

# Define the model!
# TODO: Fix the melting temperature to -1.96 degrees Celsius
# Verify that the minimum temperature is -1.96 degrees Celsius
sea_ice = sea_ice_simulation(sea_ice_grid; dynamics=nothing, advection=nothing) 

#####
##### Initialize Ocean and Sea ice models
#####

temperature = ECCOMetadata(:temperature; dir="./")
salinity    = ECCOMetadata(:salinity;    dir="./")

ice_thickness     = ECCOMetadata(:sea_ice_thickness; dir="./")
ice_concentration = ECCOMetadata(:sea_ice_concentration; dir="./")

atmosphere  = JRA55PrescribedAtmosphere(arch, backend=JRA55NetCDFBackend(20))
radiation   = Radiation(ocean_albedo = LatitudeDependentAlbedo(), sea_ice_albedo=0.6)

set!(ocean.model, T=temperature, S=salinity) 
set!(sea_ice.model.ice_thickness,     ice_thickness,     inpainting=NearestNeighborInpainting(1))
set!(sea_ice.model.ice_concentration, ice_concentration, inpainting=NearestNeighborInpainting(1))

earth_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
earth = Simulation(earth_model; Δt=30minutes, stop_iteration=10, stop_time=30days)

u, v, _ = ocean.model.velocities
T = ocean.model.tracers.T
S = ocean.model.tracers.S
s = sqrt(u^2 + v^2)

h = sea_ice.model.ice_thickness
ℵ = sea_ice.model.ice_concentration
η = ocean.model.free_surface.η 

earth.output_writers[:surface_tracers] = JLD2OutputWriter(ocean.model, (; T, S, s),
                                                          schedule = TimeInterval(12hours),
                                                          indices = (:, :, grid.Nz),
                                                          overwrite_existing = true,
                                                          filename = "surface_fields.jld2")


earth.output_writers[:sea_ice_variables] = JLD2OutputWriter(sea_ice.model, (; h, ℵ),
                                                            schedule = TimeInterval(12hours),
                                                            overwrite_existing = true,
                                                            filename = "sea_ice_fields.jld2")


earth.output_writers[:free_surface] = JLD2OutputWriter(ocean.model, (; η),
                                                       schedule = TimeInterval(12hours),
                                                       overwrite_existing = true,
                                                       filename = "free_surface.jld2")

Q  = earth.model.interfaces.net_fluxes.ocean_surface.T
τx = earth.model.interfaces.net_fluxes.ocean_surface.u
τy = earth.model.interfaces.net_fluxes.ocean_surface.v
PE = earth.model.interfaces.net_fluxes.ocean_surface.S

earth.output_writers[:fluxes] = JLD2OutputWriter(ocean.model, (; Q, τx, τy, PE),
                                                 schedule = TimeInterval(12hours),
                                                 overwrite_existing = true,
                                                 filename = "surface_fluxes.jld2")

# Also, we add a callback to print a message about how the simulation is going

wall_time = [time_ns()]

function progress(earth)
    clock   = earth.model.clock

    maxu = maximum(abs, u)
    maxv = maximum(abs, v)
    maxT = maximum(T)
    minS = minimum(S)
    
    @info @sprintf("Iteration: %d, time: %s, wall_time: %s, max(|u|, |v|): %.2e %.2e max(T): %.2e, min(S): %.2e\n",
                   clock.iteration, prettytime(clock.time), prettytime(1e-9 * (time_ns() - wall_time[1])), maxu, maxv, maxT, minS)

    wall_time[1] = time_ns()
end

add_callback!(earth, progress, IterationInterval(10))


run!(earth)

# ## Visualizing the results
#
# We can visualize the results using CairoMakie. We record a video of surface variables and fluxes.
# To load the data we can use Oceananigans' `FieldTimeSeries` object.