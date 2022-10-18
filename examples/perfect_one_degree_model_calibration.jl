using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: WallTimeInterval
using Oceananigans.BuoyancyModels: buoyancy
using Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField
using ClimaOcean.NearGlobalSimulations: one_degree_near_global_simulation

simulation = one_degree_near_global_simulation(CPU(),
                                               start_time = 0.0,
                                               with_isopycnal_skew_symmetric_diffusivity = true, 
                                               isopycnal_κ_skew = 2000.0,
                                               isopycnal_κ_symmetric = 2000.0,
                                               stop_iteration = 2)

@info "Running simulation..."; start_time = time_ns()

run!(simulation)

@info "... finished. (" * prettytime(1e-9 * (time_ns() - start_time)) * ")"

