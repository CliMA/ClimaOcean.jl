import ClimaOcean: atmosphere_simulation

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
thermodynamics_parameters(atmos::SpeedySimulation) = 
    ClimaOcean.OceanSeaIceModels.AtmosphereThermodynamicsParameters(Float32)

function initialize_atmospheric_state!(simulation::SpeedyWeather.Simulation)
    progn, diagn, model  = SpeedyWeather.unpack(simulation)
    (; time) = progn.clock  # current time

    # set the tendencies back to zero for accumulation
    fill!(diagn.tendencies, 0, typeof(model))

    if model.physics                   
        SpeedyWeather.parameterization_tendencies!(diagn, progn, time, model)
    end
    
    return nothing
end

function atmosphere_simulation(spectral_grid::SpectralGrid; output=false)
    # Surface fluxes
    humidity_flux_ocean = PrescribedOceanHumidityFlux(spectral_grid)
    humidity_flux_land = SurfaceLandHumidityFlux(spectral_grid)
    surface_humidity_flux = SurfaceHumidityFlux(ocean=humidity_flux_ocean, land=humidity_flux_land)

    ocean_heat_flux = PrescribedOceanHeatFlux(spectral_grid)
    land_heat_flux = SurfaceLandHeatFlux(spectral_grid)
    surface_heat_flux = SurfaceHeatFlux(ocean=ocean_heat_flux, land=land_heat_flux)

    # The atmospheric model
    atmosphere_model = PrimitiveWetModel(spectral_grid;
                                        surface_heat_flux,
                                        surface_humidity_flux,
                                        ocean = nothing,
                                        sea_ice = nothing) # This is provided by ClimaSeaIce

    # Construct the simulation
    atmosphere = initialize!(atmosphere_model)

    # Initialize the simulation
    initialize!(atmosphere; output)

    # Fill in prognostic fields
    initialize_atmospheric_state!(atmosphere)

    return atmosphere
end
