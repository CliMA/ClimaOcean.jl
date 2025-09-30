using ClimaOcean
using ClimaSeaIce
using Oceananigans
using Oceananigans.Grids
using Oceananigans.Units
using Oceananigans.OrthogonalSphericalShellGrids
using ClimaOcean.OceanSimulations
using ClimaOcean.JRA55
using ClimaOcean.DataWrangling
using ClimaOcean.DataWrangling: NearestNeighborInpainting
using ClimaSeaIce.SeaIceThermodynamics: IceWaterThermalEquilibrium
using Printf
using Dates
using CUDA
using PythonCall
using Oceananigans.BuoyancyFormulations: buoyancy, buoyancy_frequency

import Oceananigans.OutputWriters: checkpointer_address

arch = GPU()
# arch = Distributed(GPU(), partition=Partition(1, 4), synchronized_communication=true)

Nx = 2880 # longitudinal direction 
Ny = 1440 # meridional direction 
Nz = 80

z_faces = ExponentialCoordinate(Nz, -6000, 0)
# z_faces = MutableVerticalDiscretization(z_faces)

const z_surf = z_faces(Nz)

@info "Building grid..."
grid = TripolarGrid(arch;
                    size = (Nx, Ny, Nz),
                    z = z_faces,
                    halo = (7, 7, 7))

@info "Regridding bathymetry..."
bottom_height = regrid_bathymetry(grid; minimum_depth=15, major_basins=1, interpolation_passes=10)

@info "Building immersed boundary grid..."
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)

#####
##### A Propgnostic Ocean model
#####

using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEVerticalDiffusivity, CATKEMixingLength, CATKEEquation
using Oceananigans.TurbulenceClosures: RiBasedVerticalDiffusivity

momentum_advection = WENOVectorInvariant()
tracer_advection   = WENO(order=7)

free_surface = SplitExplicitFreeSurface(grid; cfl=0.8, fixed_Δt=12minutes)

obl_closure = ClimaOcean.OceanSimulations.default_ocean_closure() # CATKE
# obl_closure = RiBasedVerticalDiffusivity()
closure = (obl_closure, VerticalScalarDiffusivity(κ=1e-5, ν=1e-4))

if obl_closure isa RiBasedVerticalDiffusivity
    prefix = "RiBased"
else
    prefix = "CATKE"
end

prefix *= "oneeighth_degree"

glorys_dir = joinpath(homedir(), "GLORYS_data")
mkpath(glorys_dir)

glorys_dataset = GLORYSMonthly()

# Turn off salinity restoring for now
# @inline mask(x, y, z, t) = z ≥ z_surf - 1
# Smetadata = Metadata(:salinity; dataset=glorys_dataset, dir=glorys_dir)
# FS = DatasetRestoring(Smetadata, grid; rate = 1/18days, mask, time_indices_in_memory = 3)

@info "Building ocean component..."
ocean = ocean_simulation(grid; Δt=1minutes,
                         momentum_advection,
                         tracer_advection,
                         timestepper = :SplitRungeKutta3,
                         free_surface,
                        #  forcing = (; S = FS),
                         closure)

start_date = DateTime(1993, 1, 1)
end_date   = DateTime(2003, 4, 1)
simulation_period = Dates.value(Second(end_date - start_date))

inpainting = NearestNeighborInpainting(50)
@info "Setting initial conditions..."
set!(ocean.model, T=Metadatum(:temperature; dataset=glorys_dataset, date=start_date, dir=glorys_dir, inpainting),
                  S=Metadatum(:salinity;    dataset=glorys_dataset, date=start_date, dir=glorys_dir, inpainting))

@info ocean.model.clock

#####
##### A Prognostic Sea-ice model
#####

@info "Building sea-ice component..."
# Default sea-ice dynamics and salinity coupling are included in the defaults
# sea_ice = sea_ice_simulation(grid, ocean; advection=WENO(order=7))
sea_ice = sea_ice_simulation(grid, ocean; dynamics=nothing)

@info "Setting sea-ice initial conditions..."
set!(sea_ice.model, h=Metadatum(:sea_ice_thickness;     dataset=glorys_dataset, dir=glorys_dir, inpainting=nothing),
                    ℵ=Metadatum(:sea_ice_concentration; dataset=glorys_dataset, dir=glorys_dir, inpainting=nothing))

#####
##### A Prescribed Atmosphere model
#####

jra55_dir = joinpath(homedir(), "JRA55_data")
mkpath(jra55_dir)
dataset = MultiYearJRA55()
# jra55_dataset = RepeatYearJRA55()
jra55_backend = JRA55NetCDFBackend(10)

@info "Building atmospheric forcing..."
atmosphere = JRA55PrescribedAtmosphere(arch; dir=jra55_dir, dataset=jra55_backend, backend=jra55_backend, include_rivers_and_icebergs=true, start_date)
radiation  = Radiation()

#####
##### An ocean-sea ice coupled model
#####

@info "Building coupled ocean-sea ice model..."
omip = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
# omip = Simulation(omip, Δt=10minutes, stop_time=60days)
omip = Simulation(omip, Δt=10minutes, stop_time=simulation_period)

# Figure out the outputs....
checkpointer_address(::SeaIceModel) = "SeaIceModel"

FILE_DIR = "./Data/$(prefix)/"
mkpath(FILE_DIR)

@info "Setting up output writers..."
ocean.output_writers[:checkpointer] = Checkpointer(ocean.model,
                                                  schedule = SpecifiedTimes([simulation_period]),
                                                  prefix = "$(FILE_DIR)/ocean_checkpoint_oneeighthdegree",
                                                  overwrite_existing = true)

sea_ice.output_writers[:checkpointer] = Checkpointer(sea_ice.model,
                                                     schedule = TimeInterval(60days),
                                                     prefix = "$(FILE_DIR)/sea_ice_checkpoint_oneeighthdegree",
                                                     overwrite_existing = true)

sea_ice.output_writers[:checkpointer] = Checkpointer(sea_ice.model,
                                                     schedule = SpecifiedTimes([simulation_period]),
                                                     prefix = "$(FILE_DIR)/sea_ice_checkpoint_oneeighthdegree_final",
                                                     overwrite_existing = true)
                                                     
u, v, w = ocean.model.velocities
T, S = ocean.model.tracers
b = Field(buoyancy(ocean.model))
N² = Field(buoyancy_frequency(ocean.model))

ocean_outputs = merge(ocean.model.tracers, ocean.model.velocities, (; b, N²))

sea_ice_outputs = merge((h = sea_ice.model.ice_thickness,
                         ℵ = sea_ice.model.ice_concentration,
                         T = sea_ice.model.ice_thermodynamics.top_surface_temperature),
                         sea_ice.model.velocities)

ocean.output_writers[:surface] = JLD2Writer(ocean.model, ocean_outputs;
                                            schedule = TimeInterval(30days),
                                            filename = "$(FILE_DIR)/ocean_surface_fields",
                                            indices = (:, :, grid.Nz),
                                            overwrite_existing = true)

sea_ice.output_writers[:surface] = JLD2Writer(ocean.model, sea_ice_outputs;
                                              schedule = TimeInterval(30days),
                                              filename = "$(FILE_DIR)/sea_ice_surface_fields",
                                              overwrite_existing = true)

ocean.output_writers[:full] = JLD2Writer(ocean.model, ocean_outputs;
                                         schedule = SpecifiedTimes([simulation_period]),
                                         filename = "$(FILE_DIR)/ocean_complete_fields",
                                         overwrite_existing = true)

sea_ice.output_writers[:full] = JLD2Writer(sea_ice.model, sea_ice_outputs;
                                            schedule = SpecifiedTimes([simulation_period]),
                                            filename = "$(FILE_DIR)/sea_ice_complete_fields",
                                            overwrite_existing = true)

ocean.output_writers[:monthly_average] = JLD2Writer(ocean.model, ocean_outputs;
                                                    schedule = AveragedTimeInterval(30days, window=30days),
                                                    filename = "$(FILE_DIR)/ocean_complete_fields_monthly_average",
                                                    overwrite_existing = true)

sea_ice.output_writers[:monthly_average] = JLD2Writer(sea_ice.model, sea_ice_outputs;
                                                   schedule = AveragedTimeInterval(30days, window=30days),
                                                   filename = "$(FILE_DIR)/sea_ice_complete_fields_monthly_average",
                                                   overwrite_existing = true)

ocean.output_writers[:yearly_average] = JLD2Writer(ocean.model, ocean_outputs;
                                                   schedule = AveragedTimeInterval(365days, window=365days),
                                                   filename = "$(FILE_DIR)/ocean_complete_fields_yearly_average",
                                                   overwrite_existing = true)

sea_ice.output_writers[:yearly_average] = JLD2Writer(sea_ice.model, sea_ice_outputs;
                                                   schedule = AveragedTimeInterval(365days, window=365days),
                                                   filename = "$(FILE_DIR)/sea_ice_complete_fields_yearly_average",
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
    umax = maximum(ocean.model.velocities.u)
    vmax = maximum(ocean.model.velocities.v)
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
add_callback!(omip, progress, IterationInterval(100))

@info "Starting simulation..."
run!(omip)

# omip.Δt = 10minutes
# omip.stop_time = simulation_period

run!(omip)

#%%
@info "Plotting results..."
using CairoMakie

uo = FieldTimeSeries("$(FILE_DIR)/ocean_surface_fields.jld2",  "u"; backend = OnDisk())
vo = FieldTimeSeries("$(FILE_DIR)/ocean_surface_fields.jld2",  "v"; backend = OnDisk())
To = FieldTimeSeries("$(FILE_DIR)/ocean_surface_fields.jld2",  "T"; backend = OnDisk())

# and sea ice fields with "i":
ui = FieldTimeSeries("$(FILE_DIR)/sea_ice_surface_fields.jld2", "u"; backend = OnDisk())
vi = FieldTimeSeries("$(FILE_DIR)/sea_ice_surface_fields.jld2", "v"; backend = OnDisk())
hi = FieldTimeSeries("$(FILE_DIR)/sea_ice_surface_fields.jld2", "h"; backend = OnDisk())
ℵi = FieldTimeSeries("$(FILE_DIR)/sea_ice_surface_fields.jld2", "ℵ"; backend = OnDisk())
Ti = FieldTimeSeries("$(FILE_DIR)/sea_ice_surface_fields.jld2", "T"; backend = OnDisk())

times = uo.times
Nt = length(times)
n = Observable(Nt)

# We create a land mask and use it to fill land points with `NaN`s.
land = interior(To.grid.immersed_boundary.bottom_height) .≥ 0

Toₙ = @lift begin
    Tₙ = interior(To[$n])
    Tₙ[land] .= NaN
    view(Tₙ, :, :, 1)
end

heₙ = @lift begin
    hₙ = interior(hi[$n])
    ℵₙ = interior(ℵi[$n])
    hₙ[land] .= NaN
    view(hₙ, :, :, 1) .* view(ℵₙ, :, :, 1)
end

# We compute the surface speeds for the ocean and the sea ice.
uoₙ = Field{Face, Center, Nothing}(uo.grid)
voₙ = Field{Center, Face, Nothing}(vo.grid)

uiₙ = Field{Face, Center, Nothing}(ui.grid)
viₙ = Field{Center, Face, Nothing}(vi.grid)

so = Field(sqrt(uoₙ^2 + voₙ^2))
si = Field(sqrt(uiₙ^2 + viₙ^2))

soₙ = @lift begin
    parent(uoₙ) .= parent(uo[$n])
    parent(voₙ) .= parent(vo[$n])
    compute!(so)
    soₙ = interior(so)
    soₙ[land] .= NaN
    view(soₙ, :, :, 1)
end

siₙ = @lift begin
    parent(uiₙ) .= parent(ui[$n])
    parent(viₙ) .= parent(vi[$n])
    compute!(si)
    siₙ = interior(si)
    hₙ = interior(hi[$n])
    ℵₙ = interior(ℵi[$n])
    he = hₙ .* ℵₙ
    siₙ[he .< 1e-7] .= 0
    siₙ[land] .= NaN
    view(siₙ, :, :, 1)
end

# Finally, we plot a snapshot of the surface speed, temperature, and the turbulent
# eddy kinetic energy from the CATKE vertical mixing parameterization as well as the
# sea ice speed and the effective sea ice thickness.
fig = Figure(size=(2000, 1000))

title = @lift string("Global 1/8ᵒ ocean simulation after ", prettytime(times[$n] - times[1]))

axso = Axis(fig[1, 1])
axsi = Axis(fig[1, 3])
axTo = Axis(fig[2, 1])
axhi = Axis(fig[2, 3])

hmo = heatmap!(axso, soₙ, colorrange = (0, 0.5), colormap = :deep,  nan_color=:lightgray)
hmi = heatmap!(axsi, siₙ, colorrange = (0, 0.5), colormap = :greys, nan_color=:lightgray)
Colorbar(fig[1, 2], hmo, label = "Ocean Surface speed (m s⁻¹)")
Colorbar(fig[1, 4], hmi, label = "Sea ice speed (m s⁻¹)")

hmo = heatmap!(axTo, Toₙ, colorrange = (-1, 32), colormap = :magma, nan_color=:lightgray)
hmi = heatmap!(axhi, heₙ, colorrange =  (0, 4),  colormap = :blues, nan_color=:lightgray)
Colorbar(fig[2, 2], hmo, label = "Surface Temperature (ᵒC)")
Colorbar(fig[2, 4], hmi, label = "Effective ice thickness (m)")

for ax in (axso, axsi, axTo, axhi)
    hidedecorations!(ax)
end

Label(fig[0, :], title)

save("$(FILE_DIR)/global_snapshot.png", fig)
nothing #hide

# ![](global_snapshot.png)

# And now a movie:

CairoMakie.record(fig, "$(FILE_DIR)/oneeighth_degree_global_ocean_surface.mp4", 1:Nt, framerate = 8) do nn
    n[] = nn
end
nothing
#%%