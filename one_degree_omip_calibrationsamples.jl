using ClimaOcean
using ClimaSeaIce
using Oceananigans
using Oceananigans.Grids
using Oceananigans.Units
using Oceananigans.OrthogonalSphericalShellGrids
using Oceananigans.BuoyancyFormulations: buoyancy, buoyancy_frequency
using ClimaOcean.OceanSimulations
using ClimaOcean.ECCO
using ClimaOcean.JRA55
using ClimaOcean.DataWrangling
using ClimaSeaIce.SeaIceThermodynamics: IceWaterThermalEquilibrium
using Printf
using Dates
using CUDA
using JLD2
using ArgParse
using Oceananigans.OutputWriters: AveragedSpecifiedTimes

import Oceananigans.OutputWriters: checkpointer_address

using Libdl
ucx_libs = filter(lib -> occursin("ucx", lowercase(lib)), Libdl.dllist())
if isempty(ucx_libs)
    @info "✓ No UCX - safe to run!"
else
    @warn "✗ UCX libraries detected! This can cause issues with MPI+CUDA. Detected libs:\n$(join(ucx_libs, "\n"))"
end

function parse_commandline()
    s = ArgParseSettings()
  
    @add_arg_table! s begin
      "--kappa_skew"
        help = "Isopycnal skew diffusivity (m^2/s)"
        arg_type = Float64
        default = 1e3
      "--kappa_symmetric"
        help = "Isopycnal skew diffusivity (m^2/s)"
        arg_type = Float64
        default = 1e3
      "--start_year"
        help = "Start year of the simulation"
        arg_type = Int
        default = 1992
      "--simulation_length"
        help = "Length of the simulation in years"
        arg_type = Int
        default = 25
      "--sampling_length"
        help = "Length of the sampling window in years"
        arg_type = Int
        default = 10
    end
    return parse_args(s)
end

args = parse_commandline()
κ_skew = args["kappa_skew"]
κ_symmetric = args["kappa_symmetric"]
start_year = args["start_year"]
simulation_length = args["simulation_length"]
sampling_length = args["sampling_length"]

@info "1-degree omip"
@info "Using κ_skew = $(κ_skew) m²/s and κ_symmetric = $(κ_symmetric) m²/s, starting in year $(start_year) for a length of $(simulation_length) years with a $(sampling_length)-year sample."

function synch!(clock1::Clock, clock2)
    # Synchronize the clocks
    clock1.time = clock2.time
    clock1.iteration = clock2.iteration
    clock1.last_Δt = clock2.last_Δt
end

synch!(model1, model2) = synch!(model1.clock, model2.clock)
# restart_iteration = "446000" 
restart_iteration = nothing

arch = GPU()

Nx = 360 # longitudinal direction 
Ny = 180 # meridional direction 
Nz = 100

z_faces = ExponentialDiscretization(Nz, -6000, 0; scale=1800)
const z_surf = z_faces(Nz)

grid = TripolarGrid(arch;
                    size = (Nx, Ny, Nz),
                    z = z_faces,
                    halo = (7, 7, 7))

bottom_height = regrid_bathymetry(grid; minimum_depth=15, major_basins=1, interpolation_passes=75)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)

#####
##### A Propgnostic Ocean model
#####

using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization, AdvectiveFormulation, IsopycnalSkewSymmetricDiffusivity
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEVerticalDiffusivity, CATKEMixingLength, CATKEEquation

momentum_advection = WENOVectorInvariant(order=5)
tracer_advection   = WENO(order=5)
free_surface       = SplitExplicitFreeSurface(grid; cfl=0.8, fixed_Δt=50minutes)

using Oceananigans.Operators: Δx, Δy

eddy_closure  = IsopycnalSkewSymmetricDiffusivity(; κ_skew, κ_symmetric, skew_flux_formulation=AdvectiveFormulation())
# obl_closure = ClimaOcean.OceanSimulations.default_ocean_closure()
obl_closure = RiBasedVerticalDiffusivity()

closure = (obl_closure, VerticalScalarDiffusivity(κ=1e-5, ν=3e-4), eddy_closure)

prefix = "onedegree"
if obl_closure isa RiBasedVerticalDiffusivity
    prefix *= "_RiBased"
else
    prefix *= "_CATKE"
end

prefix *= "_$(κ_skew)_$(κ_symmetric)"
prefix *= "_$(start_year)"
prefix *= "_$(simulation_length)year_$(sampling_length)yearsample"
prefix *= "_advectiveGM_multiyearjra55_calibrationsamples"

dir = joinpath(homedir(), "EN4_data")
mkpath(dir)

start_date = DateTime(start_year, 1, 1)
end_date = start_date + Year(simulation_length)
simulation_period = Dates.value(Second(end_date - start_date))
yearly_times = cumsum(vcat([0.], Dates.value.(Second.(diff(start_date:Year(1):end_date)))))
decadal_times = cumsum(vcat([0.], Dates.value.(Second.(diff(start_date:Year(10):end_date)))))
# sampling_endtimes = decadal_times[3:end]
sampling_start_date = end_date - Year(sampling_length)
sampling_window = Dates.value(Second(end_date - sampling_start_date))

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

dir = joinpath(homedir(), "ECCO_data")

set!(sea_ice.model, h=Metadatum(:sea_ice_thickness;     dataset=ECCO4Monthly(), dir),
                    ℵ=Metadatum(:sea_ice_concentration; dataset=ECCO4Monthly(), dir))

@info "Initialized sea ice fields"

#####
##### A Prescribed Atmosphere model
#####
jra55_dir = joinpath(homedir(), "JRA55_data")
mkpath(jra55_dir)
dataset = MultiYearJRA55()
backend = JRA55NetCDFBackend(100)

@info "Setting up presctibed atmosphere $(dataset)"
atmosphere = JRA55PrescribedAtmosphere(arch; dir=jra55_dir, dataset, backend, include_rivers_and_icebergs=true, start_date, end_date)
radiation  = Radiation()

@info "Built atmosphere model $(atmosphere)"

#####
##### An ocean-sea ice coupled model
#####

omip = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

@info "Built coupled model $(omip)"

omip = Simulation(omip, Δt=40minutes, stop_time=simulation_period) 
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
                                         schedule = AveragedTimeInterval(3650days, window=3650days),
                                         filename = "$(FILE_DIR)/ocean_complete_fields_10year_average",
                                         overwrite_existing = true)

sea_ice.output_writers[:time_average] = JLD2Writer(sea_ice.model, sea_ice_outputs;
                                                   schedule = AveragedTimeInterval(3650days, window=3650days),
                                                   filename = "$(FILE_DIR)/sea_ice_complete_fields_10year_average",
                                                   overwrite_existing = true)

ocean.output_writers[:sample_decadal_average] = JLD2Writer(ocean.model, ocean_outputs;
                                         schedule = AveragedTimeInterval(simulation_period, window=sampling_window),
                                         filename = "$(FILE_DIR)/ocean_complete_fields_$(sampling_length)year_average_calibrationsample",
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