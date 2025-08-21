# Make sure the atmospheric parameters from SpeedyWeather can be used in the compute fluxes function
import ClimaOcean.OceanSeaIceModels.PrescribedAtmospheres: 
    thermodynamics_parameters, 
    boundary_layer_height, 
    surface_layer_height

const SpeedySimulation = SpeedyWeather.Simulation
const SpeedyCoupledModel = ClimaOcean.OceanSeaIceModel{<:Any, <:SpeedySimulation}
Base.summary(::SpeedySimulation) = "SpeedyWeather.Simulation"

# Take one time-step
Oceananigans.TimeSteppers.time_step!(atmos::SpeedySimulation, Δt) = SpeedyWeather.timestep!(atmos)

# The height of near-surface variables used in the turbulent flux solver
function surface_layer_height(s::SpeedySimulation)
    T = s.model.atmosphere.temp_ref
    g = s.model.planet.gravity
    Φ = s.model.geopotential.Δp_geopot_full
    return Φ[end] * T / g
end

# This is a parameter that is used in the computation of the fluxes,
# It probably should not be here but in the similarity theory type.
boundary_layer_height(atmos::SpeedySimulation) = 600

# This is a _hack_!! The parameters should be consistent with what is specified in SpeedyWeather
thermodynamics_parameters(atmos::SpeedyWeather.Simulation) = 
    ClimaOcean.OceanSeaIceModels.AtmosphereThermodynamicsParameters(Float32)



