using Oceananigans
using Oceananigans.Architectures: arch_array
using Oceananigans.Units
using Oceananigans.Utils: WallTimeInterval
using Oceananigans.BuoyancyModels: buoyancy
using Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField
using ClimaOcean.NearGlobalSimulations: one_degree_near_global_simulation
using ParameterEstimocean
using ParameterEstimocean.Observations: FieldTimeSeriesCollector
using ParameterEstimocean.Parameters: closure_with_parameters

prefix = "perfect_one_degree_calibration"
start_time = 345days
stop_iteration = 10

simulation_kw = (; start_time, stop_iteration,
                 isopycnal_κ_skew = 900.0,
                 isopycnal_κ_symmetric = 900.0)

test_simulation = one_degree_near_global_simulation(; simulation_kw...)

# Test...
# new_parameters = (κ_skew = 900.0, κ_symmetric = 900.0, max_slope=1e-3)
# new_closure = closure_with_parameters(test_simulation.model.closure, new_parameters)
# @show new_closure
# @show new_closure[4]

model = test_simulation.model

test_simulation.output_writers[:d3] = JLD2OutputWriter(model, model.tracers,
                                                       schedule = IterationInterval(stop_iteration),
                                                       filename = prefix * "_fields",
                                                       overwrite_existing = true)

slice_indices = (11, :, :)
test_simulation.output_writers[:d2] = JLD2OutputWriter(model, model.tracers,
                                                       schedule = IterationInterval(stop_iteration),
                                                       filename = prefix * "_slices",
                                                       indices = slice_indices,
                                                       overwrite_existing = true)


@info "Running simulation..."; timer = time_ns()

run!(test_simulation)

@info "... finished. (" * prettytime(1e-9 * (time_ns() - timer)) * ")"

#=
@show test_simulation.model.tracers.T
@show test_simulation.model.tracers.S

T_slices = FieldTimeSeries(prefix * "_slices.jld2", "T")
S_slices = FieldTimeSeries(prefix * "_slices.jld2", "S")
=#

simulation_ensemble = [one_degree_near_global_simulation(; simulation_kw...) for i = 1:5]

priors = (κ_skew = ScaledLogitNormal(bounds=(0.0, 2000.0)),
          κ_symmetric = ScaledLogitNormal(bounds=(0.0, 2000.0)))

free_parameters = FreeParameters(priors) 

obspath = prefix * "_slices.jld2"
times = [start_time, time(test_simulation)]
observations = SyntheticObservations(obspath; field_names=(:T, :S), times)
 
# Initial condition
T₀ = FieldTimeSeries(prefix * "_fields.jld2", "T")[1]
S₀ = FieldTimeSeries(prefix * "_fields.jld2", "S")[1]

T₀_GPU = arch_array(GPU(), parent(T₀))
S₀_GPU = arch_array(GPU(), parent(S₀))

function initialize_simulation!(sim, parameters)
    model = sim.model
    parent(sim.model.velocities.u) .= 0 
    parent(sim.model.velocities.v) .= 0 
    T, S = sim.model.tracers
    T .= T₀_GPU
    S .= S₀_GPU
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

time_series_collector_ensemble = [slice_collector(sim) for sim in simulation_ensemble]

ip = InverseProblem(observations, simulation_ensemble, free_parameters;
                    time_series_collector = time_series_collector_ensemble,
                    initialize_with_observations = false,
                    initialize_simulation = initialize_simulation!)

#θ = [(κ_skew=900.0 + randn(), κ_symmetric=900.0+randn()) for sim in simulation_ensemble]
#forward_run!(ip, θ)

eki = EnsembleKalmanInversion(ip; pseudo_stepping=ConstantConvergence(0.2))
iterate!(eki, iterations=10)

#=
fig = Figure()
ax = Axis(fig[1, 1])

for iter in [0, 1, 4, 10]
    summary = eki.iteration_summaries[iter]
    κh = map(θ -> θ.κh, summary.parameters)
    κz = map(θ -> θ.κz, summary.parameters)
    scatter!(ax, κh, κz, label="Iteration $iter")
end

axislegend(ax)

display(fig)
=#
