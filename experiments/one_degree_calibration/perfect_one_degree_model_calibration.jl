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
using JLD2

arch = GPU()

#####
##### Simulation parameters
#####

prefix = "perfect_one_degree_calibration"
start_time = 345days
stop_iteration = 1000

simulation_kw = (; start_time, stop_iteration,
                 isopycnal_κ_skew = 900.0,
                 isopycnal_κ_symmetric = 900.0)

slice_indices = (11, :, :)

#####
##### Benchmark simulation
#####

initialization = false

if initialization
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
end

simulation_ensemble = [one_degree_near_global_simulation(arch; simulation_kw...) for _ in 1:4]

priors = (κ_skew      = ScaledLogitNormal(bounds=(0.0, 2000.0)),
          κ_symmetric = ScaledLogitNormal(bounds=(0.0, 2000.0)))

free_parameters = FreeParameters(priors) 

obspath = prefix * "_slices.jld2"
 
T₀ = FieldTimeSeries(prefix * "_fields.jld2", "T")
S₀ = FieldTimeSeries(prefix * "_fields.jld2", "S")

times = T₀.times
observations = SyntheticObservations(obspath; field_names=(:T, :S), times)

# Initial conditions
T_init = Array(interior(T₀[1]))
S_init = Array(interior(S₀[1]))

function initialize_simulation!(sim, parameters)
    fill!(sim.model.velocities.u, 0.0)
    fill!(sim.model.velocities.v, 0.0)
    T, S = sim.model.tracers
    fill!(T, 0.0)
    fill!(S, 0.0)
    set!(T, T_init)
    set!(S, S_init)
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
    return FieldTimeSeriesCollector((T=T_slice, S=S_slice), times, architecture = CPU())
end

##### 
##### Building the inverse problem
#####

time_series_collector_ensemble = [slice_collector(sim) for sim in simulation_ensemble]

ip = InverseProblem(observations, simulation_ensemble, free_parameters;
                    time_series_collector = time_series_collector_ensemble,
                    initialize_with_observations = false,
                    initialize_simulation = initialize_simulation!)

eki = EnsembleKalmanInversion(ip; pseudo_stepping=ConstantConvergence(0.2))

##### 
##### Let's run!
#####

iterate!(eki, iterations=10)

@info "final parameters: $(eki.iteration_summaries[end].ensemble_mean)"
