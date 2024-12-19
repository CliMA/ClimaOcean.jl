using Oceananigans
using Oceananigans.Architectures: arch_array
using Oceananigans.Units
using Oceananigans.Utils: WallTimeInterval
using Oceananigans.BuoyancyFormulations: buoyancy
using Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField
using ClimaOcean.NearGlobalSimulations: one_degree_near_global_simulation
using ParameterEstimocean
using ParameterEstimocean.Utils: map_gpus_to_ranks!
using ParameterEstimocean.Observations: FieldTimeSeriesCollector
using ParameterEstimocean.Parameters: closure_with_parameters

using MPI
using CUDA
using JLD2

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

initial_conditions_path = datadep"near_global_one_degree/initial_conditions_month_01_360_150_48.jld2",

prefix = "perfect_one_degree_calibration"
start_time = 345days
stop_iteration = 100

simulation_kw = (; start_time, stop_iteration,
                 isopycnal_κ_skew = 900.0,
                 isopycnal_κ_symmetric = 900.0,
                 initial_conditions_path)

slice_indices = (11, :, :)

#####
##### Benchmark simulation only on rank (and GPU) 0
#####

# Reinitialize the problem
initialization = false

if rank == 0 && initialization
    test_simulation = one_degree_near_global_simulation(arch; simulation_kw...)

    model = test_simulation.model

    test_simulation.output_writers[:d3] = JLD2OutputWriter(model, model.tracers,
                                                           schedule = IterationInterval(100),
                                                           filename = prefix * "_fields",
                                                           overwrite_existing = true)

    slice_indices = (11, :, :)
    test_simulation.output_writers[:d2] = JLD2OutputWriter(model, model.tracers,
                                                           schedule = IterationInterval(100),
                                                           filename = prefix * "_slices",
                                                           indices = slice_indices,
                                                           overwrite_existing = true)


    @info "Running simulation..."; timer = time_ns()

    run!(test_simulation)

    @info "... finished. (" * prettytime(1e-9 * (time_ns() - timer)) * ")"
    stop_time = [time(test_simulation)]
else 
    stop_time = stop_iteration * 20minutes + start_time
end

# Finished simulation, wait rank 0
MPI.Barrier(comm)

for r in 0:nproc-1
    rank == r && @info "rank $rank on device $(CUDA.device())"
    MPI.Barrier(comm)
end

#####
##### On all ranks
#####
          
simulation = one_degree_near_global_simulation(arch; simulation_kw...) 

priors = (κ_skew      = ScaledLogitNormal(bounds=(0.0, 2000.0)),
          κ_symmetric = ScaledLogitNormal(bounds=(0.0, 2000.0)))

free_parameters = FreeParameters(priors) 

obspath = prefix * "_slices.jld2"

times = [start_time, stop_time]
observations = SyntheticObservations(obspath; field_names=(:T, :S), times)

# Initial conditions
T_init = jldopen(initial_conditions_path)["T"]
S_init = jldopen(initial_conditions_path)["S"]

function initialize_simulation!(sim, parameters)
    @info "initializing model on rank $(MPI.Comm_rank(MPI.COMM_WORLD))"
    fill!(sim.model.velocities.u, 0.0)
    fill!(sim.model.velocities.v, 0.0)
    T, S = sim.model.tracers
    fill!(T, 0.0)
    fill!(S, 0.0)
    set!(T, T_init)
    set!(S, S_init)
    return nothing
end

function slice_collector(sim)
    T = sim.model.tracers.T
    S = sim.model.tracers.S
    T_slice = Field(T; indices=slice_indices)
    S_slice = Field(S; indices=slice_indices)
    return FieldTimeSeriesCollector((T=T_slice, S=S_slice), times, architecture = CPU())
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

@show eki.iteration_summaries[end]
