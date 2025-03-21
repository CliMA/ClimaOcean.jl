# A stress tests for the atmospheric sea ice flux solver

using Thermodynamics
using ClimaSeaIce
using ClimaOcean.OceanSeaIceModels: PrescribedAtmosphereThermodynamicsParameters
using ClimaOcean.OceanSeaIceModels.InterfaceComputations: InterfaceState, 
                                                          InterfaceProperties, 
                                                          SpecificHumidityFormulation, 
                                                          CoefficientBasedFluxes,
                                                          RelativeVelocity, 
                                                          compute_interface_state,
                                                          atmosphere_sea_ice_stability_functions,
                                                          SimilarityTheoryFluxes

# Building states (surface, interior, atmosphere)
u★ = 1e-4
Ψs = InterfaceState(u★, u★, u★, 0.0, 0.0, 268.97, 0.0, 0.00268863)
Ψi = (; u = 0.0, v = 0.0, T = 271.53, S = 30.0, h = 0.2, ℵ = 1.0)
Ψa = (; z = 10.0, 
        u = 10.0, 
        v = 0.0, 
        𝒬 = Thermodynamics.PhaseEquil{Float64}(1.2910431086728036, 101325.0, -71648.65771694925, 0.0028001374147162287, 272.9929607369637), 
        h_bℓ = 600.0)

# Downwelling radiation
Rd = (; Qs = 100.0, Qℓ = 250.0)

# Building propertities (surface, interior, atmosphere)
Qc = ConductiveFlux(conductivity=2)
Ts = SkinTemperature(Qc)
ub = RelativeVelocity()
Rp = (; α = 0.6, σ = 5.67e-8, ϵ = 1.0)
qs = SpecificHumidityFormulation(Thermodynamics.Ice())
Ps = InterfaceProperties(Rp, qs, Ts, ub)

TP = PrescribedAtmosphereThermodynamicsParameters(Float64)
Pa = (; thermodynamics_parameters = TP, surface_layer_height = 10.0)

lq = ClimaSeaIce.SeaIceThermodynamics.LinearLiquidus()
tu = ClimaOcean.OceanSeaIceModels.InterfaceComputations.DegreesCelsius()
Pi = (; reference_density = 900.0, heat_capacity = 2100.0, freshwater_density = 1000.0, liquidus = lq, temperature_units = tu)

# Flux solver
solver_tolerance = 1e-8
solver_maxiter = 1000000
stability_functions = atmosphere_sea_ice_stability_functions()
solver = SimilarityTheoryFluxes(; solver_tolerance, solver_maxiter, stability_functions)

# Compute the interface state and check for convergence
compute_interface_state(solver, Ψs, Ψa, Ψi, Rd, Ps, Pa, Pi)