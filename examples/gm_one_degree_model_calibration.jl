using Oceananigans
using Oceananigans.Architectures: arch_array
using Oceananigans.Units
using Oceananigans.Utils: WallTimeInterval
using Oceananigans.BuoyancyModels: buoyancy
using Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField
using ClimaOcean.NearGlobalSimulations: one_degree_near_global_simulation
using ParameterEstimocean
using ParameterEstimocean.Utils: map_gpus_to_ranks!
using ParameterEstimocean.Observations: FieldTimeSeriesCollector
using ParameterEstimocean.Parameters: closure_with_parameters

using MPI
using CUDA

MPI.Init()

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

prefix = "gm_one_degree_calibration"
start_time = 0
stop_time  = 45days
slice_indices = (11, :, :)

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
 
T₁ = jldopen(comparison_file)["T"]
S₁ = jldopen(comparison_file)["S"]

simulation_kw = (; start_time, stop_time, 
                   isopycnal_κ_skew = 900.0,
                   isopycnal_κ_symmetric = 900.0,
                   initial_condition_fields = (T = T₀, S = S₀))

simulation = one_degree_near_global_simulation(arch; simulation_kw...) 

priors = (κ_skew      = ScaledLogitNormal(bounds=(0.0, 2000.0)),
          κ_symmetric = ScaledLogitNormal(bounds=(0.0, 2000.0)))

free_parameters = FreeParameters(priors) 

times = [start_time, stop_time]

Tfield = FieldTimeSeries(on_architecture(CPU(), simulation.model.grid), times; indices = slice_indices)
Sfield = FieldTimeSeries(on_architecture(CPU(), simulation.model.grid), times; indices = slice_indices)

set!(Tfield[1], T₀)
set!(Sfield[1], S₀)
set!(Tfield[2], T₁)
set!(Sfield[2], S₁)

observations = SyntheticObservations(; field_names=(:T, :S), field_time_serieses = (; T = Tfield, S = Sfield))

function initialize_simulation!(sim, parameters)
    fill!(sim.model.velocities.u, 0.0)
    fill!(sim.model.velocities.v, 0.0)
    T, S = sim.model.tracers
    fill!(T, 0.0)
    fill!(S, 0.0)
    set!(T, T₀)
    set!(S, S₀)
    return nothing
end

for sim in simulation_ensemble
    initialize_simulation!(sim, nothing)
    run!(sim)
end

function slice_collector(sim)
    T = sim.model.tracers.T
    S = sim.model.tracers.S
    T_slice = Field(T; indices=slice_indices)
    S_slice = Field(S; indices=slice_indices)
    return FieldTimeSeriesCollector((T=T_slice, S=S_slice), times, architecture = CPU(), averaging_window=30days)
end

##### 
##### Building the distributed inverse problem
#####

time_series_collector = slice_collector(simulation) 

ip = InverseProblem(observations, simulation, free_parameters;
                    time_series_collector,
                    initialize_with_observations = false,
                    initialize_simulation = initialize_simulation!)

dip = DistributedInverseProblems(ip)

eki = EnsembleKalmanInversion(dip; pseudo_stepping=ConstantConvergence(0.2))
@info "finished setting up eki"

##### 
##### Let's run!
#####

iterate!(eki, iterations=10)

@info "final parameters: $(eki.iteration_summaries[end].ensemble_mean)"
