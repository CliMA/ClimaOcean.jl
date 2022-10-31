using Oceananigans
using Oceananigans.Architectures: arch_array
using Oceananigans.Units
using Oceananigans.Grids: on_architecture
using Oceananigans.Utils: WallTimeInterval
using Oceananigans.BuoyancyModels: buoyancy
using Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField
using ClimaOcean.NearGlobalSimulations: one_degree_near_global_simulation
using ParameterEstimocean
using ParameterEstimocean.Utils: map_gpus_to_ranks!
using ParameterEstimocean.Observations: FieldTimeSeriesCollector
using ParameterEstimocean.Parameters: closure_with_parameters
using DataDeps
using JLD2 

using MPI
using CUDA

MPI.Init()

using Random
Random.seed!(123)

#####
##### Setting up multi architecture infrastructure
#####

comm = MPI.COMM_WORLD

rank  = MPI.Comm_rank(comm)
nproc = MPI.Comm_size(comm)

arch = GPU()

if arch isa GPU
    map_gpus_to_ranks!()
end

#####
##### Simulation parameters
#####

prefix = "gm_one_degree_calibration_"
start_time = 0
stop_time  = 45days
slice_indices = (UnitRange(11, 11), :, :)

# Finished simulation, wait rank 0
MPI.Barrier(comm)

for r in 0:nproc-1
    rank == r && @info "rank $rank on device $(CUDA.device())"
    MPI.Barrier(comm)
end

#####
##### On all ranks
#####
 
initial_condition = datadep"near_global_one_degree/initial_conditions_month_01_360_150_48.jld2"
comparison_file   = datadep"near_global_one_degree/initial_conditions_month_02_360_150_48.jld2"

T₀ = jldopen(initial_condition)["T"]
S₀ = jldopen(initial_condition)["S"]
 
T₀[isnan.(T₀)] .= 0.0
S₀[isnan.(S₀)] .= 0.0

T₁ = jldopen(comparison_file)["T"]
S₁ = jldopen(comparison_file)["S"]

T₁[isnan.(T₁)] .= 0.0
T₁[isnan.(S₁)] .= 0.0

simulation_kw = (; start_time, stop_time, 
                   isopycnal_κ_skew = 900.0,
                   isopycnal_κ_symmetric = 900.0)

simulation = one_degree_near_global_simulation(arch; simulation_kw...) 

T, S = simulation.model.tracers

dir = "./" 

save_indices = Dict(
    :depth_5_meters    => (:,   :, 48),
    :depth_45_meters   => (:,   :, 44),
    :depth_257_meters  => (:,   :, 30),
    :depth_1007_meters => (:,   :, 20),
    :pacific_transect  => (11,  :, :),
    :atlantic_transect => (150, :, :),
    :southern_ocean    => (:,  15, :)) 

eki_iteration = 0

function initialize_output_writers!(sim, iteration, rank)
    model = sim.model
    T, S  = model.tracers
    for (name, idx) in save_indices
        delete!(sim.output_writers, name)

        output_prefix = prefix * string(name) * "_eki_iteration" * string(iteration) * "_particle$(rank)"
        sim.output_writers[name] = JLD2OutputWriter(model, (; T, S); dir,
                                                    schedule = TimeInterval(44days),
                                                    filename = output_prefix,
                                                    indices = idx,
                                                    overwrite_existing = true)
    end
end

priors = (κ_skew      = ScaledLogitNormal(bounds=(0.0, 4000.0)),
          κ_symmetric = ScaledLogitNormal(bounds=(0.0, 4000.0)))

free_parameters = FreeParameters(priors) 

times = [start_time, stop_time]

Tfield = FieldTimeSeries{Center, Center, Center}(on_architecture(CPU(), simulation.model.grid), times; indices = slice_indices)
Sfield = FieldTimeSeries{Center, Center, Center}(on_architecture(CPU(), simulation.model.grid), times; indices = slice_indices)

interior(Tfield[1]) .= T₀[slice_indices...]
interior(Sfield[1]) .= S₀[slice_indices...]
interior(Tfield[2]) .= T₁[slice_indices...]
interior(Sfield[2]) .= S₁[slice_indices...]

observations = SyntheticObservations(; field_names=(:T, :S), field_time_serieses = (; T = Tfield, S = Sfield))

function initialize_simulation!(sim, parameters)
    fill!(sim.model.velocities.u, 0.0)
    fill!(sim.model.velocities.v, 0.0)
    T, S = sim.model.tracers
    fill!(T, 0.0)
    fill!(S, 0.0)
    set!(T, T₀)
    set!(S, S₀)
    sim.model.clock.time = start_time
    global eki_iteration += 1
    initialize_output_writers!(sim, eki_iteration, rank)
    return nothing
end

function slice_collector(sim)
    T = sim.model.tracers.T
    S = sim.model.tracers.S
    T_slice = Field(T; indices=slice_indices)
    S_slice = Field(S; indices=slice_indices)
    return FieldTimeSeriesCollector((T=T_slice, S=S_slice), times, architecture = CPU(), averaging_window = 30days)
end

##### 
##### Building the distributed inverse problem
#####

time_series_collector = slice_collector(simulation) 

ip = InverseProblem(observations, simulation, free_parameters;
                    time_series_collector,
                    initialize_with_observations = false,
                    initialize_simulation = initialize_simulation!)

dip = DistributedInverseProblem(ip)

eki = EnsembleKalmanInversion(dip; pseudo_stepping=ConstantConvergence(0.2))
@info "finished setting up eki"

##### 
##### Let's run!
#####

iterate!(eki, iterations=10)

using JLD2
jldsave("calibration-$rank.jld2", eki = eki)

@show eki.iteration_summaries[end]
