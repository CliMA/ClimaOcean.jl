using ClimaOcean
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids
using OrthogonalSphericalShellGrids
using CFTime
using Dates
using Printf
using CairoMakie

using ClimaOcean.DataWrangling: NearestNeighborInpainting

arch = CPU()

depth = 5000meters
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
sea_ice_grid = ImmersedBoundaryGrid(sea_ice_grid, GridFittedBottom(bottom_height); active_cells_map=true)

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

# ocean_velocities = (u = interior(ocean.model.velocities.u, :, :, grid.Nz:grid.Nz), 
#                     v = interior(ocean.model.velocities.u, :, :, grid.Nz:grid.Nz))


# τₒ = SemiImplicitStress(uₑ = ocean_velocities.u, vₑ = ocean_velocities.v)
# τᵤₐ
# # We use an elasto-visco-plastic rheology and WENO seventh order 
# # for advection of h and ℵ
# momentum_equations = SeaIceMomentumEquation(grid; 
#                                             top_momentum_stress = (u = τᵤₐ, v = τᵥₐ),
#                                             bottom_momentum_stress = τₒ,
#                                             coriolis = FPlane(f=1e-4),
#                                             ocean_velocities,
#                                             rheology = ElastoViscoPlasticRheology(min_substeps=50, 
#                                                                                   max_substeps=100,
#                                                                                   minimum_plastic_stress=1e-10),
#                                             solver = SplitExplicitSolver(substeps=150))

# Define the model!
sea_ice = sea_ice_simulation(sea_ice_grid; dynamics=nothing, advection=nothing) 

#####
##### Initialize Ocean and Sea ice models
#####

temperature = ECCOMetadata(:temperature; dir="./")
salinity    = ECCOMetadata(:salinity;    dir="./")

ice_thickness     = ECCOMetadata(:sea_ice_thickness; dir="./")
ice_concentration = ECCOMetadata(:sea_ice_concentration; dir="./")

set!(ocean.model, T=temperature, S=salinity) 
set!(sea_ice.model.ice_thickness,     ice_thickness,     inpainting=NearestNeighborInpainting(1))
set!(sea_ice.model.ice_concentration, ice_concentration, inpainting=NearestNeighborInpainting(1))

atmosphere  = JRA55PrescribedAtmosphere(arch, backend=JRA55NetCDFBackend(20))
radiation   = Radiation(ocean_albedo = LatitudeDependentAlbedo(), sea_ice_albedo=0.6)
earth_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

earth = Simulation(earth_model; Δt=30minutes, stop_time=30days)

u, v, _ = ocean.model.velocities
T = ocean.model.tracers.T
S = ocean.model.tracers.S
s = sqrt(u^2 + v^2)

η = ocean.model.free_surface.η 

earth.output_writers[:surface_tracers] = JLD2OutputWriter(ocean.model, (; T, S, s),
                                                          schedule = TimeInterval(12hours),
                                                          indices = (:, :, grid.Nz),
                                                          overwrite_existing = true,
                                                          filename = "surface_fields.jld2")


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

using JLD2
using Oceananigans
using Oceananigans.Grids: halo_size
using CairoMakie 
using Statistics: mean

file  = jldopen("free_surface.jld2")
iters = keys(file["timeseries/t"]) 

Hx, Hy, _ = halo_size(η.grid)
T  = FieldTimeSeries("surface_fields.jld2", "T")
S  = FieldTimeSeries("surface_fields.jld2", "S")
s  = FieldTimeSeries("surface_fields.jld2", "s")

n  = Observable(1)
Tn = @lift(interior(T[$n], :, :, 1))
Sn = @lift(interior(S[$n], :, :, 1))
sn = @lift(interior(s[$n], :, :, 1))
ηn = @lift(file["timeseries/η/" * iters[$n]][Hx+1:end-Hx, Hy+1:end-Hy, 1])

fig = Figure(size = (1800, 800))
axT = Axis(fig[1, 1], title="Surface temperature ᵒC")
axS = Axis(fig[1, 2], title="Surface salinity psu")
axs = Axis(fig[2, 1], title="Surface speed ms⁻¹")
axη = Axis(fig[2, 2], title="Sea surface height m")

λ, φ, z = nodes(T[1])

hmT = heatmap!(axT, Tn, colormap=:magma,  colorrange=(-1, 30))
hmS = heatmap!(axS, Sn, colormap=:haline, colorrange=(25, 40))
hms = heatmap!(axs, sn, colormap=:deep,   colorrange=( 0, 0.8))
hmη = heatmap!(axη, ηn, colormap=:bwr,    colorrange=(-1, 1))

CairoMakie.record(fig, "surface_fields.mp4", 1:length(T.times); framerate=5) do i 
    @info "doing $i of $(length(T.times))"
    n[] = i
end

# let's also visualize the surface fluxes that force the model

Q  = FieldTimeSeries("surface_fluxes.jld2", "Q")
τx = FieldTimeSeries("surface_fluxes.jld2", "τx")
τy = FieldTimeSeries("surface_fluxes.jld2", "τy")
PE = FieldTimeSeries("surface_fluxes.jld2", "PE")

Qn  = @lift(interior(Q[$n],  :, :, 1))
τxn = @lift(interior(τx[$n], :, :, 1))
τyn = @lift(interior(τy[$n], :, :, 1))
PEn = @lift(interior(PE[$n], :, :, 1))

fig  = Figure(size = (1800, 800))
axQ  = Axis(fig[1, 1], title="Net heat flux Wm⁻²")
axPE = Axis(fig[1, 2], title="Net salt flux psu m s⁻¹")
axτx = Axis(fig[2, 1], title="Zonal wind stress Nm⁻²")
axτy = Axis(fig[2, 2], title="Meridional wind stress Nm⁻²")

hmQ  = heatmap!(axQ,  Qn,  colormap=:magma,   colorrange=(-800,  800))
hmPE = heatmap!(axPE, PEn, colormap=:haline,  colorrange=(-1e-5, 5e-5))
hmτx = heatmap!(axτx, τxn, colormap=:balance, colorrange=(-5e-4, 5e-4))
hmτy = heatmap!(axτy, τyn, colormap=:balance, colorrange=(-5e-4, 5e-4))

CairoMakie.record(fig, "surface_fluxes.mp4", 1:length(Q.times); framerate=5) do i 
    @info "doing $i of $(length(Q.times))"
    n[] = i
end
