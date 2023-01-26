using ClimaOcean.NearGlobalSimulations: one_degree_near_global_simulation

using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: WallTimeInterval
using Oceananigans.BuoyancyModels: buoyancy
using Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField
using JLD2

# Choose closure
boundary_layer_turbulence_closure = ConvectiveAdjustmentVerticalDiffusivity(convective_κz = 1.0)

@show boundary_layer_turbulence_closure

#####
##### Build the simulation
#####

start_time = 345days
stop_time  = 435days 
with_isopycnal_skew_symmetric_diffusivity = false

initial_conditions               = "input/initial_conditions.jld2"
bathymetry_path                  = "input/bathy_t4.jld2"
surface_boundary_conditions_path = "input/surface_boundary_conditions.jld2"

using ClimaOcean.NearGlobalSimulations: u_bottom_linear_drag, v_bottom_linear_drag
using ClimaOcean.NearGlobalSimulations: u_immersed_bottom_linear_drag, v_immersed_bottom_linear_drag

simulation = one_degree_near_global_simulation(CPU(); 
    start_time, 
    stop_time,
    with_isopycnal_skew_symmetric_diffusivity,
    boundary_layer_turbulence_closure,
    surface_boundary_conditions_path,
    initial_conditions,
    u_drag = u_bottom_linear_drag,
    v_drag = u_bottom_linear_drag,
    u_immersed_drag = u_immersed_bottom_linear_drag,
    v_immersed_drag = u_immersed_bottom_linear_drag,
    bottom_drag_coefficient = 1e-3,
    horizontal_diffusivity = 300.0,
    tracer_advection = Centered,
    bathymetry_path,
    equation_of_state = LinearEquationOfState()
)

output_prefix = "linear_test"

fields_save_interval = 5days

simulation.output_writers[:checkpointer] = Checkpointer(simulation.model; 
                                                        prefix = output_prefix * "_checkpointer",
                                                        schedule = WallTimeInterval(10minutes),
                                                        cleanup = true,
                                                        overwrite_existing = true)

model = simulation.model

simulation.output_writers[:fields] = JLD2OutputWriter(model, merge(model.velocities, model.tracers);
                                                      schedule = TimeInterval(fields_save_interval),
                                                      filename = output_prefix * "_fields",
                                                      with_halos = true,
                                                      overwrite_existing = true)

run!(simulation)
