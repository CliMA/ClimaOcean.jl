using StaticArrays
using Thermodynamics
using OffsetArrays

using ..OceanSeaIceModels: reference_density,
                           heat_capacity,
                           sea_ice_concentration,
                           sea_ice_thickness,
                           thermodynamics_parameters

using ClimaSeaIce: SeaIceModel

using Oceananigans: HydrostaticFreeSurfaceModel, architecture
using Oceananigans.Units: Time
using Oceananigans.Grids: inactive_node, node, topology
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Fields: ConstantField, interpolate, FractionalIndices
using Oceananigans.Utils: launch!, KernelParameters
using Oceananigans.Operators: ℑxᶜᵃᵃ, ℑyᵃᶜᵃ, ℑxᶠᵃᵃ, ℑyᵃᶠᵃ

using KernelAbstractions: @kernel, @index

import Oceananigans.Simulations: initialize!

#####
##### Container for organizing information related to fluxes
#####

mutable struct AtmosphereInterface{J, F, ST, P}
    fluxes :: J
    flux_formulation :: F
    temperature :: ST
    properties :: P
end

mutable struct SeaIceOceanInterface{J, P}
    fluxes :: J
    properties :: P
end

# Utilities to get the computed fluxes
@inline computed_fluxes(interface::AtmosphereInterface)  = interface.fluxes
@inline computed_fluxes(interface::SeaIceOceanInterface) = interface.fluxes
@inline computed_fluxes(::Nothing) = nothing

@inline get_possibly_zero_flux(fluxes, name)    = getfield(fluxes, name)
@inline get_possibly_zero_flux(::Nothing, name) = ZeroField()

mutable struct ComponentInterfaces{AO, ASI, SIO, C, AP, OP, SIP, EX, P}
    atmosphere_ocean_interface :: AO
    atmosphere_sea_ice_interface :: ASI
    sea_ice_ocean_interface :: SIO
    atmosphere_properties :: AP
    ocean_properties :: OP
    sea_ice_properties :: SIP
    exchanger :: EX
    net_fluxes :: C
    properties :: P
end

# Possible units for temperature and salinity
struct DegreesCelsius end
struct DegreesKelvin end

const celsius_to_kelvin = 273.15
@inline convert_to_kelvin(::DegreesCelsius, T::FT) where FT = T + convert(FT, celsius_to_kelvin)
@inline convert_to_kelvin(::DegreesKelvin, T) = T

@inline convert_from_kelvin(::DegreesCelsius, T::FT) where FT = T - convert(FT, celsius_to_kelvin)
@inline convert_from_kelvin(::DegreesKelvin, T) = T

Base.summary(crf::ComponentInterfaces) = "ComponentInterfaces"
Base.show(io::IO, crf::ComponentInterfaces) = print(io, summary(crf))

#####
##### Atmosphere-Ocean Interface
#####

atmosphere_ocean_interface(grid, ::Nothing,   ocean,    args...) = nothing
atmosphere_ocean_interface(grid, ::Nothing,  ::Nothing, args...) = nothing
atmosphere_ocean_interface(grid, atmosphere, ::Nothing, args...) = nothing

function atmosphere_ocean_interface(grid, 
                                    atmosphere,
                                    ocean,
                                    radiation,
                                    ao_flux_formulation,
                                    temperature_formulation,
                                    velocity_formulation,
                                    specific_humidity_formulation)

    water_vapor           = Field{Center, Center, Nothing}(grid)
    latent_heat           = Field{Center, Center, Nothing}(grid)
    sensible_heat         = Field{Center, Center, Nothing}(grid)
    x_momentum            = Field{Center, Center, Nothing}(grid)
    y_momentum            = Field{Center, Center, Nothing}(grid)
    friction_velocity     = Field{Center, Center, Nothing}(grid)
    temperature_scale     = Field{Center, Center, Nothing}(grid)
    water_vapor_scale     = Field{Center, Center, Nothing}(grid)
    upwelling_longwave    = Field{Center, Center, Nothing}(grid)
    downwelling_longwave  = Field{Center, Center, Nothing}(grid)
    downwelling_shortwave = Field{Center, Center, Nothing}(grid)

    ao_fluxes = (; latent_heat,
                   sensible_heat,
                   water_vapor,
                   x_momentum,
                   y_momentum,
                   friction_velocity,
                   temperature_scale,
                   water_vapor_scale,
                   upwelling_longwave,
                   downwelling_longwave,
                   downwelling_shortwave)

    σ = radiation.stefan_boltzmann_constant
    αₐₒ = radiation.reflection.ocean
    ϵₐₒ = radiation.emission.ocean
    radiation = (σ=σ, α=αₐₒ, ϵ=ϵₐₒ)

    ao_properties = InterfaceProperties(radiation,
                                        specific_humidity_formulation,
                                        temperature_formulation,
                                        velocity_formulation)

    interface_temperature = Field{Center, Center, Nothing}(grid)

    return AtmosphereInterface(ao_fluxes, ao_flux_formulation, interface_temperature, ao_properties)
end

#####
##### Atmosphere-Sea Ice Interface
#####

atmosphere_sea_ice_interface(grid, atmos, ::Nothing,     args...) = nothing
atmosphere_sea_ice_interface(grid, ::Nothing, sea_ice,   args...) = nothing
atmosphere_sea_ice_interface(grid, ::Nothing, ::Nothing, args...) = nothing

function atmosphere_sea_ice_interface(grid, 
                                      atmosphere,
                                      sea_ice,
                                      radiation,
                                      ai_flux_formulation,
                                      temperature_formulation,
                                      velocity_formulation)

    water_vapor   = Field{Center, Center, Nothing}(grid)
    latent_heat   = Field{Center, Center, Nothing}(grid)
    sensible_heat = Field{Center, Center, Nothing}(grid)
    x_momentum    = Field{Center, Center, Nothing}(grid)
    y_momentum    = Field{Center, Center, Nothing}(grid)
    fluxes = (; latent_heat, sensible_heat, water_vapor, x_momentum, y_momentum)

    σ   = radiation.stefan_boltzmann_constant
    αₐᵢ = radiation.reflection.sea_ice
    ϵₐᵢ = radiation.emission.sea_ice
    radiation = (σ=σ, α=αₐᵢ, ϵ=ϵₐᵢ)

    phase = AtmosphericThermodynamics.Ice()
    specific_humidity_formulation = ImpureSaturationSpecificHumidity(phase)

    properties = InterfaceProperties(radiation,
                                     specific_humidity_formulation,
                                     temperature_formulation,
                                     velocity_formulation)

    interface_temperature = sea_ice.model.ice_thermodynamics.top_surface_temperature

    return AtmosphereInterface(fluxes, ai_flux_formulation, interface_temperature, properties)
end

#####
##### Sea Ice-Ocean Interface
#####

sea_ice_ocean_interface(grid, ::Nothing, ocean;     kwargs...) = nothing
sea_ice_ocean_interface(grid, ::Nothing, ::Nothing; kwargs...) = nothing
sea_ice_ocean_interface(grid, sea_ice,   ::Nothing; kwargs...) = nothing

function sea_ice_ocean_interface(grid, sea_ice, ocean; characteristic_melting_speed = 1e-5)

    io_bottom_heat_flux = Field{Center, Center, Nothing}(grid)
    io_frazil_heat_flux = Field{Center, Center, Nothing}(grid)
    io_salt_flux = Field{Center, Center, Nothing}(grid)
    x_momentum = Field{Face, Center, Nothing}(grid)
    y_momentum = Field{Center, Face, Nothing}(grid)

    io_fluxes = (interface_heat=io_bottom_heat_flux,
                 frazil_heat=io_frazil_heat_flux,
                 salt=io_salt_flux,
                 x_momentum=x_momentum,
                 y_momentum=y_momentum)

    io_properties = (; characteristic_melting_speed)

    return SeaIceOceanInterface(io_fluxes, io_properties)
end

#####
##### Component Interfaces
#####

default_ai_temperature(::Nothing) = nothing

function default_ai_temperature(sea_ice)
    conductive_flux = sea_ice.model.ice_thermodynamics.internal_heat_flux.parameters.flux
    return SkinTemperature(conductive_flux)
end

function default_ao_specific_humidity(ocean)
    FT    = eltype(ocean)
    phase = AtmosphericThermodynamics.Liquid()
    x_H₂O = convert(FT, 0.98)
    return ImpureSaturationSpecificHumidity(phase, x_H₂O)
end

default_exchange_grid(atmosphere, ocean, sea_ice) = ocean.model.grid

"""
    ComponentInterfaces(atmosphere, ocean, sea_ice=nothing;
                        radiation = Radiation(),
                        freshwater_density = default_freshwater_density,
                        atmosphere_ocean_fluxes = SimilarityTheoryFluxes(),
                        atmosphere_sea_ice_fluxes = SimilarityTheoryFluxes(eltype(ocean.model.grid)),
                        atmosphere_ocean_interface_temperature = BulkTemperature(),
                        atmosphere_ocean_interface_specific_humidity = default_ao_specific_humidity(ocean),
                        atmosphere_sea_ice_interface_temperature = default_ai_temperature(sea_ice),
                        ocean_reference_density = reference_density(ocean),
                        ocean_heat_capacity = heat_capacity(ocean),
                        ocean_temperature_units = DegreesCelsius(),
                        sea_ice_temperature_units = DegreesCelsius(),
                        sea_ice_reference_density = reference_density(sea_ice),
                        sea_ice_heat_capacity = heat_capacity(sea_ice),
                        gravitational_acceleration = default_gravitational_acceleration)
"""
function ComponentInterfaces(atmosphere, ocean, sea_ice=nothing;
                             exchange_grid = default_exchange_grid(atmosphere, ocean, sea_ice),
                             radiation = Radiation(),
                             freshwater_density = default_freshwater_density,
                             atmosphere_ocean_fluxes = SimilarityTheoryFluxes(eltype(exchange_grid)),
                             atmosphere_sea_ice_fluxes = atmosphere_sea_ice_similarity_theory(eltype(exchange_grid)),
                             atmosphere_ocean_interface_temperature = BulkTemperature(),
                             atmosphere_ocean_velocity_difference = RelativeVelocity(),
                             atmosphere_ocean_interface_specific_humidity = default_ao_specific_humidity(ocean),
                             atmosphere_sea_ice_interface_temperature = default_ai_temperature(sea_ice),
                             atmosphere_sea_ice_velocity_difference = RelativeVelocity(),
                             ocean_reference_density = reference_density(ocean),
                             ocean_heat_capacity = heat_capacity(ocean),
                             ocean_temperature_units = DegreesCelsius(),
                             sea_ice_temperature_units = DegreesCelsius(),
                             sea_ice_reference_density = reference_density(sea_ice),
                             sea_ice_heat_capacity = heat_capacity(sea_ice),
                             gravitational_acceleration = default_gravitational_acceleration)

    FT = eltype(exchange_grid)

    ocean_reference_density    = convert(FT, ocean_reference_density)
    ocean_heat_capacity        = convert(FT, ocean_heat_capacity)
    sea_ice_reference_density  = convert(FT, sea_ice_reference_density)
    sea_ice_heat_capacity      = convert(FT, sea_ice_heat_capacity)
    freshwater_density         = convert(FT, freshwater_density)
    gravitational_acceleration = convert(FT, gravitational_acceleration)

    # Component properties
    atmosphere_properties = thermodynamics_parameters(atmosphere)

    ocean_properties = (reference_density  = ocean_reference_density,
                        heat_capacity      = ocean_heat_capacity,
                        freshwater_density = freshwater_density,
                        temperature_units  = ocean_temperature_units)

    # Only build sea_ice_properties if sea_ice is an actual Simulation with a model
    if sea_ice isa Simulation
        sea_ice_properties = (reference_density  = sea_ice_reference_density,
                              heat_capacity      = sea_ice_heat_capacity,
                              freshwater_density = freshwater_density,
                              liquidus           = sea_ice.model.ice_thermodynamics.phase_transitions.liquidus,
                              temperature_units  = sea_ice_temperature_units)
    else
        sea_ice_properties = nothing
    end

    # Component interfaces
    ao_interface = atmosphere_ocean_interface(exchange_grid,
                                              atmosphere,
                                              ocean,
                                              radiation,
                                              atmosphere_ocean_fluxes,
                                              atmosphere_ocean_interface_temperature,
                                              atmosphere_ocean_velocity_difference,
                                              atmosphere_ocean_interface_specific_humidity)

    io_interface = sea_ice_ocean_interface(exchange_grid, sea_ice, ocean)

    ai_interface = atmosphere_sea_ice_interface(exchange_grid, 
                                                atmosphere,
                                                sea_ice,
                                                radiation,
                                                atmosphere_sea_ice_fluxes,
                                                atmosphere_sea_ice_interface_temperature,
                                                atmosphere_sea_ice_velocity_difference)
    # Total interface fluxes
    total_fluxes = (ocean      = net_fluxes(ocean),
                    sea_ice    = net_fluxes(sea_ice),
                    atmosphere = net_fluxes(atmosphere))

    exchanger = StateExchanger(exchange_grid, atmosphere, ocean, sea_ice)

    properties = (; gravitational_acceleration)

    return ComponentInterfaces(ao_interface,
                               ai_interface,
                               io_interface,
                               atmosphere_properties,
                               ocean_properties,
                               sea_ice_properties,
                               exchanger,
                               total_fluxes,
                               properties)
end
