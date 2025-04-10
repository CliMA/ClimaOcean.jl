# Make sure the atmospheric parameters from SpeedyWeather can be used in the compute fluxes function
import ClimaOcean.OceanSeaIceModels.PrescribedAtmospheres: 
    thermodynamics_parameters, 
    boundary_layer_height, 
    surface_layer_height

import ClimaOcean: atmosphere_simulation

# This can be left blank:
update_model_field_time_series!(::SpeedySimulation, time) = nothing

# Take one time-step
time_step!(atmos::SpeedySimulation, Δt) = SpeedyWeather.timestep!(atmos)

function atmosphere_simulation(grid::SpeedyWeather.SpectralGrid)
                               # orography = zeros(grid.Grid, grid.nlat_half))

    surface_heat_flux = SpeedyWeather.PrescribedOceanHeatFlux()
    surface_evaporation = SpeedyWeather.PrescribedOceanEvaporation()
    atmos_model = SpeedyWeather.PrimitiveWetModel(grid;
                                                  # orography,
                                                  surface_heat_flux,
                                                  surface_evaporation)

    atmos_sim = SpeedyWeather.initialize!(atmos_model)
    return atmos_sim
end

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

# Base.eltype(::EarthAtmosphere{FT}) where FT = FT

# This is a _hack_!! The parameters should be consistent with what is specified in SpeedyWeather
thermodynamics_parameters(atmos::SpeedyWeather.Simulation) = 
    ClimaOcean.OceanSeaIceModels.PrescribedAtmosphereThermodynamicsParameters(Float32)



