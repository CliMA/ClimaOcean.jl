using Pkg
Pkg.update()
using ClimaOcean
using ClimaSeaIce
using Oceananigans
using Oceananigans.Grids
using Oceananigans.Grids: AbstractGrid
using Oceananigans.Units
using Oceananigans.OrthogonalSphericalShellGrids
using ClimaOcean.OceanSimulations
using ClimaOcean.ECCO
using ClimaOcean.JRA55
using ClimaOcean.DataWrangling
using ClimaSeaIce.SeaIceThermodynamics: IceWaterThermalEquilibrium
using Printf
using Dates
using CUDA
using Oceananigans.Fields: ConstantField
using Oceananigans.Operators

import Oceananigans.OutputWriters: checkpointer_address
import Oceananigans.OutputWriters: saveproperty!

saveproperty!(file, address, p::ConstantField) = file[address] = p

checkpointer_address(::SeaIceModel) = "SeaIceModel"

function synch!(clock1::Clock, clock2)
    # Synchronize the clocks
    clock1.time = clock2.time
    clock1.iteration = clock2.iteration
    clock1.last_Δt = clock2.last_Δt
end

synch!(model1, model2) = synch!(model1.clock, model2.clock)

arch    = GPU()
r_faces = ClimaOcean.exponential_z_faces(; Nz=60, depth=6200)
z_faces = MutableVerticalDiscretization(r_faces)

Nx = 2160 # longitudinal direction 
Ny = 1080 # meridional direction 
Nz = length(r_faces) - 1

grid = TripolarGrid(arch;
                    size = (Nx, Ny, Nz),
                    z = z_faces,
                    halo = (7, 7, 7))

bottom_height = regrid_bathymetry(grid; minimum_depth=15, major_basins=1, interpolation_passes=10)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)

#####
##### A Propgnostic Ocean model
#####

using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEVerticalDiffusivity, CATKEMixingLength, CATKEEquation

momentum_advection = WENOVectorInvariant()
tracer_advection   = WENO(order=7)

free_surface = SplitExplicitFreeSurface(grid; cfl=0.7, fixed_Δt=20minutes)

# closure = (catke_closure, VerticalScalarDiffusivity(κ=1e-5, ν=1e-5))
# closure = (ClimaOcean.OceanSimulations.default_ocean_closure(), VerticalScalarDiffusivity(κ=1e-5, ν=1e-5))

closure = RiBasedVerticalDiffusivity()

ocean = ocean_simulation(grid; Δt=1minutes,
                         momentum_advection,
                         tracer_advection,
                         timestepper = :SplitRungeKutta3,
                         free_surface,
                         closure)

dataset = ECCO4Monthly()

set!(ocean.model, T=Metadatum(:temperature; dataset),
                  S=Metadatum(:salinity;    dataset))

#####
##### A Prognostic Sea-ice model
#####

sea_ice = sea_ice_simulation(grid, ocean; dynamics = nothing) #advection=WENO(order=7))

set!(sea_ice.model, h=Metadatum(:sea_ice_thickness;     dataset),
                    ℵ=Metadatum(:sea_ice_concentration; dataset))

# using JLD2
# file = jldopen("sea_ice_checkpoint_iteration10000.jld2")
# parent(sea_ice.model.ice_thickness)     .= CuArray(file["SeaIceModel/h/data"])
# parent(sea_ice.model.ice_concentration) .= CuArray(file["SeaIceModel/ℵ/data"])
# parent(sea_ice.model.velocities.u)      .= CuArray(file["SeaIceModel/u/data"])
# parent(sea_ice.model.velocities.v)      .= CuArray(file["SeaIceModel/v/data"])

#####
##### A Prescribed Atmosphere model
#####

dir = "./forcing_data"
dataset = MultiYearJRA55()
backend = JRA55NetCDFBackend(40)

atmosphere = JRA55PrescribedAtmosphere(arch; dir, dataset, backend, include_rivers_and_icebergs=true)
radiation  = Radiation(sea_ice_albedo=0.7)

#####
##### An ocean-sea ice coupled model
#####
 
omip = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
omip = Simulation(omip, Δt=10, stop_time=Inf)

synch!(sea_ice.model, ocean.model)
synch!(omip.model, ocean.model)

# Figure out the outputs....
ocean.output_writers[:checkpointer] = Checkpointer(ocean.model,
                                                  schedule = IterationInterval(10000),
                                                  prefix = "ocean_checkpoint",
                                                  overwrite_existing = true)

sea_ice.output_writers[:checkpointer] = Checkpointer(sea_ice.model,
                                                     schedule = IterationInterval(10000),
                                                     prefix = "sea_ice_checkpoint",
                                                     overwrite_existing = true)

uo = ocean.model.velocities.u
vo = ocean.model.velocities.v
wo = ocean.model.velocities.w
T, S = ocean.model.tracers
η  = ocean.model.free_surface.η

ui = sea_ice.model.velocities.u
vi = sea_ice.model.velocities.v
h  = sea_ice.model.ice_thickness
ℵ  = sea_ice.model.ice_concentration

omip.output_writers[:ocean_free_surface] = JLD2Writer(ocean.model, (; η),
                                                        schedule = TimeInterval(1days),
                                                        filename = "ocean_free_surface",
                                                        overwrite_existing = true)

omip.output_writers[:ocean_surface_fields] = JLD2Writer(ocean.model, (; uo, vo, wo, T, S),
                                                        schedule = TimeInterval(1days),
                                                        filename = "ocean_surface_fields",
                                                        indices = (:, :, grid.Nz),
                                                        overwrite_existing = true)

omip.output_writers[:sea_ice_surface_fields] = JLD2Writer(sea_ice.model, (; ui, vi, h, ℵ),
                                                          schedule = TimeInterval(1days),
                                                          filename = "sea_ice_surface_fields",
                                                          overwrite_existing = true)

wall_time = Ref(time_ns())

using Statistics

function progress(sim)
    sea_ice = sim.model.sea_ice
    ocean   = sim.model.ocean
    hmax = maximum(sea_ice.model.ice_thickness)
    ℵmax = maximum(sea_ice.model.ice_concentration)
    Tmax = maximum(sim.model.interfaces.atmosphere_sea_ice_interface.temperature)
    Tmin = minimum(sim.model.interfaces.atmosphere_sea_ice_interface.temperature)
    umax = maximum(abs, ocean.model.velocities.u)
    vmax = maximum(abs, ocean.model.velocities.v)
    uimax = maximum(abs, sea_ice.model.velocities.u)
    vimax = maximum(abs, sea_ice.model.velocities.v)
    wmax = maximum(ocean.model.velocities.w)

    step_time = 1e-9 * (time_ns() - wall_time[])

    msg1 = @sprintf("time: %s, iteration: %d, Δt: %s, ", prettytime(sim), iteration(sim), prettytime(sim.Δt))
    msg2 = @sprintf("max(h): %.2e m, max(ℵ): %.2e ", hmax, ℵmax)
    msg4 = @sprintf("extrema(T): (%.2f, %.2f) ᵒC, ", Tmax, Tmin)
    msg5 = @sprintf("maximum(u): (%.2f, %.2f, %.2f) m/s, ", umax, vmax, wmax)
    msg6 = @sprintf("wall time: %s \n", prettytime(step_time))

    @info msg1 * msg2 * msg4 * msg5 * msg6

     wall_time[] = time_ns()

     return nothing
end

# And add it as a callback to the simulation.
add_callback!(omip, progress, IterationInterval(50))
wizard = TimeStepWizard(cfl=0.7, max_change=1.1, max_Δt=20minutes)

function sea_ice_cell_advection_timescale(grid, velocities)
    u, v = velocities
    τ = KernelFunctionOperation{Center, Center, Center}(cell_advection_timescaleᶜᶜ, grid, u, v)
    return minimum(τ)
end

@inline _inverse_timescale(i, j, k, Δ⁻¹, U, topo) = @inbounds abs(U[i, j, k]) * Δ⁻¹
@inline _inverse_timescale(i, j, k, Δ⁻¹, U, topo::Flat) = 0

@inline function cell_advection_timescaleᶜᶜ(i, j, k, grid::AbstractGrid{FT, TX, TY}, u, v) where {FT, TX, TY}
    Δx⁻¹ = Δx⁻¹ᶠᶜᶜ(i, j, k, grid)
    Δy⁻¹ = Δy⁻¹ᶜᶠᶜ(i, j, k, grid)

    inverse_timescale_x = _inverse_timescale(i, j, k, Δx⁻¹, u, TX())
    inverse_timescale_y = _inverse_timescale(i, j, k, Δy⁻¹, v, TY())

    inverse_timescale = inverse_timescale_x + inverse_timescale_y 

    return 1 / inverse_timescale
end

function add_wizard!(sim)
   wizard(sim.model.ocean)
   sea_ice = sim.model.sea_ice
   Δti = 0.15 * sea_ice_cell_advection_timescale(sea_ice.model.grid, sea_ice.model.velocities)
   @info "Wizard says: ocean Δt: $(ocean.Δt), sea ice Δt: $(Δti)"
   iter = sea_ice.model.clock.iteration
   Δtw = min(ocean.Δt, Δti)

   if iter < 5000
       sim.Δt = 60
   else
       sim.Δt = Δtw
   end

   @info "Final Δt is $(sim.Δt)"
end

omip.callbacks[:wizard] = Callback(add_wizard!, IterationInterval(10))

run!(omip)
