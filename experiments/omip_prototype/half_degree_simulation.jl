using ClimaOcean
using ClimaSeaIce
using Oceananigans
using Oceananigans.Grids
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
using Oceananigans.BuoyancyFormulations: buoyancy, buoyancy_frequency
using ArgParse

import Oceananigans.OutputWriters: checkpointer_address

function parse_commandline()
    s = ArgParseSettings()
  
    @add_arg_table! s begin
      "--kappa_skew"
        help = "Isopycnal skew diffusivity (m^2/s)"
        arg_type = Float64
        default = 5e2
      "--kappa_symmetric"
        help = "Isopycnal skew diffusivity (m^2/s)"
        arg_type = Float64
        default = 5e2
    end
    return parse_args(s)
end

args = parse_commandline()
κ_skew = args["kappa_skew"]
κ_symmetric = args["kappa_symmetric"]

arch = GPU()

Nx = 720 # longitudinal direction 
Ny = 360 # meridional direction 
Nz = 100

z_faces = ExponentialDiscretization(Nz, -6000, 0)

const z_surf = z_faces(Nz)

grid = TripolarGrid(arch;
                    size = (Nx, Ny, Nz),
                    z = z_faces,
                    halo = (7, 7, 7))

bottom_height = regrid_bathymetry(grid; minimum_depth=15, major_basins=1, interpolation_passes=40)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)

@info "Built grid $(grid)"

#####
##### A Prognostic Ocean model
#####

using Oceananigans.TurbulenceClosures: RiBasedVerticalDiffusivity

momentum_advection = WENOVectorInvariant()
tracer_advection   = WENO(order=7)

free_surface = SplitExplicitFreeSurface(grid; cfl=0.8, fixed_Δt=65minutes)

gm_closure = Oceananigans.TurbulenceClosures.EddyAdvectiveClosure(; κ_skew)
redi_closure = Oceananigans.TurbulenceClosures.IsopycnalDiffusivity(; κ_symmetric)

# obl_closure = ClimaOcean.OceanSimulations.default_ocean_closure()  
obl_closure = RiBasedVerticalDiffusivity()
closure = (obl_closure, VerticalScalarDiffusivity(κ=1e-5, ν=1e-4), gm_closure, redi_closure)

prefix = "halfdegree"
if obl_closure isa RiBasedVerticalDiffusivity
    prefix *= "_RiBased"
else
    prefix *= "_CATKE"
end

prefix *= "_$(κ_skew)_$(κ_symmetric)"
prefix *= "_newgm_multiyearjra55"

dir = joinpath(homedir(), "forcing_data_half_degree")
mkpath(dir)

start_date = DateTime(1958, 1, 1)
end_date = start_date + Year(40)
simulation_period = Dates.value(Second(end_date - start_date))
yearly_times = cumsum(vcat([0.], Dates.value.(Second.(diff(start_date:Year(1):end_date)))))
decadal_times = cumsum(vcat([0.], Dates.value.(Second.(diff(start_date:Year(10):end_date)))))

@info "Settting up salinity restoring..."
@inline mask(x, y, z, t) = z ≥ z_surf - 1
Smetadata = Metadata(:salinity; dataset=EN4Monthly(), dir, start_date, end_date)
FS = DatasetRestoring(Smetadata, grid; rate = 1/18days, mask, time_indices_in_memory = 10)

ocean = ocean_simulation(grid; Δt=1minutes,
                         momentum_advection,
                         tracer_advection,
                         timestepper = :SplitRungeKutta3,
                         free_surface,
                         forcing = (; S = FS),
                         closure)

@info "Built ocean model $(ocean)"

set!(ocean.model, T=Metadatum(:temperature; dataset=EN4Monthly(), date=start_date, dir),
                  S=Metadatum(:salinity;    dataset=EN4Monthly(), date=start_date, dir))
@info "Initialized T and S"

#####
##### A Prognostic Sea-ice model
#####

# Default sea-ice dynamics and salinity coupling are included in the defaults
# sea_ice = sea_ice_simulation(grid, ocean; advection=WENO(order=7))
sea_ice = sea_ice_simulation(grid, ocean; dynamics=nothing)
@info "Built sea ice model $(sea_ice)"

set!(sea_ice.model, h=Metadatum(:sea_ice_thickness;     dataset=ECCO4Monthly(), dir),
                    ℵ=Metadatum(:sea_ice_concentration; dataset=ECCO4Monthly(), dir))

@info "Initialized sea ice fields"

#####
##### A Prescribed Atmosphere model
#####
jra55_dir = joinpath(homedir(), "JRA55_data")
dataset = MultiYearJRA55()
backend = JRA55NetCDFBackend(20)

@info "Setting up presctibed atmosphere $(dataset)"
atmosphere = JRA55PrescribedAtmosphere(arch; dir=jra55_dir, dataset, backend, include_rivers_and_icebergs=true, start_date, end_date)
radiation  = Radiation()

@info "Built atmosphere model $(atmosphere)"

#####
##### An ocean-sea ice coupled model
#####

omip = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

@info "Built coupled model $(omip)"

omip = Simulation(omip, Δt=20minutes, stop_time=60days) 
@info "Built simulation $(omip)"

# Figure out the outputs....
checkpointer_address(::SeaIceModel) = "SeaIceModel"

FILE_DIR = "./Data/$(prefix)/"
mkpath(FILE_DIR)

ocean.output_writers[:checkpointer] = Checkpointer(ocean.model,
                                                  schedule = SpecifiedTimes(decadal_times),
                                                  prefix = "$(FILE_DIR)/ocean_checkpoint",
                                                  overwrite_existing = true)

sea_ice.output_writers[:checkpointer] = Checkpointer(sea_ice.model,
                                                     schedule = SpecifiedTimes(decadal_times),
                                                     prefix = "$(FILE_DIR)/sea_ice_checkpoint",
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
                                            schedule = TimeInterval(180days),
                                            filename = "$(FILE_DIR)/ocean_surface_fields",
                                            indices = (:, :, grid.Nz),
                                            overwrite_existing = true)

sea_ice.output_writers[:surface] = JLD2Writer(ocean.model, sea_ice_outputs;
                                              schedule = TimeInterval(180days),
                                              filename = "$(FILE_DIR)/sea_ice_surface_fields",
                                              overwrite_existing = true)

ocean.output_writers[:full] = JLD2Writer(ocean.model, ocean_outputs;
                                         schedule = TimeInterval(1825days),
                                         filename = "$(FILE_DIR)/ocean_complete_fields",
                                         overwrite_existing = true)

ocean.output_writers[:time_average] = JLD2Writer(ocean.model, ocean_outputs;
                                         schedule = AveragedTimeInterval(1825days, window=1825days),
                                         filename = "$(FILE_DIR)/ocean_complete_fields_10year_average",
                                         overwrite_existing = true)

sea_ice.output_writers[:time_average] = JLD2Writer(sea_ice.model, sea_ice_outputs;
                                                   schedule = AveragedTimeInterval(1825days, window=1825days),
                                                   filename = "$(FILE_DIR)/sea_ice_complete_fields_10year_average",
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

run!(omip)

omip.Δt = 60minutes
omip.stop_time = simulation_period

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
fig = Figure(size=(1500, 1000))

title = @lift string("Global 1ᵒ ocean simulation after ", prettytime(times[$n] - times[1]))

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

CairoMakie.record(fig, "$(FILE_DIR)/global_ocean_surface.mp4", 1:Nt, framerate = 8) do nn
    n[] = nn
end
nothing

#%%
uoac = FieldTimeSeries("$(FILE_DIR)/ocean_complete_fields_10year_average.jld2", "u"; backend = OnDisk())
voac = FieldTimeSeries("$(FILE_DIR)/ocean_complete_fields_10year_average.jld2", "v"; backend = OnDisk())
woac = FieldTimeSeries("$(FILE_DIR)/ocean_complete_fields_10year_average.jld2", "w"; backend = OnDisk())
Toac = FieldTimeSeries("$(FILE_DIR)/ocean_complete_fields_10year_average.jld2", "T"; backend = OnDisk())
Soac = FieldTimeSeries("$(FILE_DIR)/ocean_complete_fields_10year_average.jld2", "S"; backend = OnDisk())

times = uoac.times
Nt = length(times)
n = Observable(Nt)

# We create a land mask and use it to fill land points with `NaN`s.
land = interior(Toac.grid.immersed_boundary.bottom_height) .≥ 0

Toacₙ = @lift begin
    Tₙ = interior(Toac[$n])
    Tₙ[land[:, :, 1], :] .= NaN
    view(Tₙ, :, :, :)
end

Soacₙ = @lift begin
    Sₙ = interior(Soac[$n])
    Sₙ[land[:, :, 1], :] .= NaN
    view(Sₙ, :, :, :)
end

surface = Toac.grid.underlying_grid.Nz
middle = findlast(x -> x <= -500, Toac.grid.underlying_grid.z.cᵃᵃᶜ)
bottom = findlast(x -> x <= -2000, Toac.grid.underlying_grid.z.cᵃᵃᶜ)

T_surfaceₙ = @lift $Toacₙ[:, :, surface]
T_middleₙ  = @lift $Toacₙ[:, :, middle]
T_bottomₙ  = @lift $Toacₙ[:, :, bottom]

S_surfaceₙ = @lift $Soacₙ[:, :, surface]
S_middleₙ  = @lift $Soacₙ[:, :, middle]
S_bottomₙ  = @lift $Soacₙ[:, :, bottom]

Tlim_surface = extrema(interior(Toac[Nt], :, :, surface))
Tlim_middle  = extrema(interior(Toac[Nt], :, :, middle))
Tlim_bottom  = extrema(interior(Toac[Nt], :, :, bottom))

Slim_surface = extrema(interior(Soac[Nt], :, :, surface)[interior(Soac[Nt], :, :, surface) .!= 0])
Slim_middle  = extrema(interior(Soac[Nt], :, :, middle)[interior(Soac[Nt], :, :, middle) .!= 0])
Slim_bottom  = extrema(interior(Soac[Nt], :, :, bottom)[interior(Soac[Nt], :, :, bottom) .!= 0])

# Finally, we plot a snapshot of the surface speed, temperature, and the turbulent
# eddy kinetic energy from the CATKE vertical mixing parameterization as well as the
# sea ice speed and the effective sea ice thickness.
fig = Figure(size=(2400, 1000))

title = @lift string("10 year average starting from year ", round((times[$n] - times[1]) / 365 / 24 / 60^2))

axTs = Axis(fig[1, 1], title = "z = $(round(Toac.grid.underlying_grid.z.cᵃᵃᶜ[surface])) m")
axTm = Axis(fig[1, 3], title = "z = $(round(Toac.grid.underlying_grid.z.cᵃᵃᶜ[middle])) m")
axTb = Axis(fig[1, 5], title = "z = $(round(Toac.grid.underlying_grid.z.cᵃᵃᶜ[bottom])) m")

axSs = Axis(fig[2, 1])
axSm = Axis(fig[2, 3])
axSb = Axis(fig[2, 5])

T_colorscheme = :turbo
S_colorscheme = :turbo

hmTs = heatmap!(axTs, T_surfaceₙ, colormap = T_colorscheme,  nan_color=:lightgray, colorrange = Tlim_surface)
hmTm = heatmap!(axTm, T_middleₙ,  colormap = T_colorscheme,  nan_color=:lightgray, colorrange = Tlim_middle)
hmTb = heatmap!(axTb, T_bottomₙ,  colormap = T_colorscheme,  nan_color=:lightgray, colorrange = Tlim_bottom)
Colorbar(fig[1, 2], hmTs)
Colorbar(fig[1, 4], hmTm)
Colorbar(fig[1, 6], hmTb, label = "Temperature (ᵒC)")

hmSs = heatmap!(axSs, S_surfaceₙ, colormap = S_colorscheme,  nan_color=:lightgray, colorrange = Slim_surface)
hmSm = heatmap!(axSm, S_middleₙ,  colormap = S_colorscheme,  nan_color=:lightgray, colorrange = Slim_middle)
hmSb = heatmap!(axSb, S_bottomₙ,  colormap = S_colorscheme,  nan_color=:lightgray, colorrange = Slim_bottom)
Colorbar(fig[2, 2], hmSs)
Colorbar(fig[2, 4], hmSm)
Colorbar(fig[2, 6], hmSb, label = "Salinity (psu)")

for ax in (axTs, axTm, axTb, axSs, axSm, axSb)
    hidedecorations!(ax)
end

Label(fig[0, :], title)

save("$(FILE_DIR)/10year_time_average_level_snapshot.png", fig)
nothing #hide

# ![](global_snapshot.png)

# And now a movie:

CairoMakie.record(fig, "$(FILE_DIR)/10year_time_average_level.mp4", 1:Nt, framerate = 1) do nn
    n[] = nn
end
nothing
#%%