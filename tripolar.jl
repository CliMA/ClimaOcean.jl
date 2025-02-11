using ClimaOcean
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids
using OrthogonalSphericalShellGrids
using CFTime
using Dates
using Printf

arch = Distributed(CPU(), partition=Partition(4, 2), synchronized_communication=true)
child_arch = CPU()

depth = 6000meters
Nz    = 100

r_faces = ClimaOcean.exponential_z_faces(; Nz, depth)
z_faces = MutableVerticalDiscretization(r_faces)

Nx = 4320 # longitudinal direction -> 250 points is about 1.5ᵒ resolution
Ny = 1800 # meridional direction -> same thing, 48 points is about 1.5ᵒ resolution
Nz   = length(r_faces) - 1
grid = TripolarGrid(arch, Float64; size=(Nx, Ny, Nz), z=z_faces)

bottom_height = regrid_bathymetry(grid; minimum_depth=15, major_basins=1, filename, dir="./")

fig = Figure(size = (800, 400))
ax  = Axis(fig[1, 1])
heatmap!(ax, interior(bottom_height, :, :, 1), colormap=:deep)
display(fig)

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)

momentum_advection = WENOVectorInvariant() 
tracer_advection   = WENO(order=7)

free_surface = SplitExplicitFreeSurface(grid; substeps=70) 

ocean = ocean_simulation(grid; 
                         momentum_advection, 
                         tracer_advection, 
                         free_surface)

temperature = ECCOMetadata(:temperature; dir="./")
salinity    = ECCOMetadata(:salinity;    dir="./")

set!(ocean.model, T=temperature, S=salinity) 

atmosphere = JRA55PrescribedAtmosphere(child_arch; backend=JRA55NetCDFBackend(10), dir="./")
radiation  = Radiation()

similarity_theory = SimilarityTheoryTurbulentFluxes(grid; maxiter=10)
earth_model = OceanSeaIceModel(ocean; atmosphere, radiation, similarity_theory)

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
                                                          filename = "surface_fields_$(arch.local_rank).jld2")


earth.output_writers[:free_surface] = JLD2OutputWriter(ocean.model, (; η),
                                                       schedule = TimeInterval(12hours),
                                                       overwrite_existing = true,
                                                       filename = "free_surface_$(arch.local_rank).jld2")

Q  = earth.model.fluxes.total.ocean.heat
τx = earth.model.fluxes.total.ocean.momentum.u
τy = earth.model.fluxes.total.ocean.momentum.v
PE = earth.model.fluxes.total.ocean.tracers.S

earth.output_writers[:fluxes] = JLD2OutputWriter(ocean.model, (; Q, τx, τy, PE),
                                                 schedule = TimeInterval(12hours),
                                                 overwrite_existing = true,
                                                 filename = "surface_fluxes_$(arch.local_rank).jld2")

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