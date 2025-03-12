using StaticArrays
using Thermodynamics
using SurfaceFluxes
using OffsetArrays

using ..OceanSeaIceModels: reference_density,
                           heat_capacity,
                           sea_ice_concentration,
                           sea_ice_thickness,
                           downwelling_radiation,
                           freshwater_flux,
                           SeaIceSimulation

using ..OceanSeaIceModels.PrescribedAtmospheres:
    PrescribedAtmosphere,
    thermodynamics_parameters

using ClimaSeaIce: SeaIceModel

using Oceananigans: HydrostaticFreeSurfaceModel, architecture
using Oceananigans.Grids: inactive_node, node, topology
using Oceananigans.BoundaryConditions: fill_halo_regions!

using Oceananigans.Fields: ConstantField, interpolate, FractionalIndices
using Oceananigans.Utils: launch!, Time, KernelParameters

using Oceananigans.Operators: ℑxᶜᵃᵃ, ℑyᵃᶜᵃ, ℑxᶠᵃᵃ, ℑyᵃᶠᵃ

using KernelAbstractions: @kernel, @index

#####
##### Container for organizing information related to fluxes
#####

struct AtmosphereInterface{J, F, ST, P}
    fluxes :: J
    flux_formulation :: F
    temperature :: ST
    properties :: P
end

struct SeaIceOceanInterface{J, P, H, A}
    fluxes :: J
    properties :: P
    previous_ice_thickness :: H
    previous_ice_concentration :: A
end

struct ComponentInterfaces{AO, ASI, SIO, C, AP, OP, SIP, EX}
    atmosphere_ocean_interface :: AO
    atmosphere_sea_ice_interface :: ASI
    sea_ice_ocean_interface :: SIO
    atmosphere_properties :: AP
    ocean_properties :: OP
    sea_ice_properties :: SIP
    exchanger :: EX
    net_fluxes :: C
end

struct StateExchanger{G, AST, AEX}
    exchange_grid :: G
    exchange_atmosphere_state :: AST
    atmosphere_exchanger :: AEX
end

# Note that Field location can also affect fractional index type.
# Here we assume that we know the location of Fields that will be interpolated.
fractional_index_type(FT, Topo) = FT
fractional_index_type(FT, ::Flat) = Nothing

function StateExchanger(ocean::Simulation, atmosphere)
    # TODO: generalize this
    exchange_grid = ocean.model.grid

    exchange_atmosphere_state = (u  = Field{Center, Center, Nothing}(exchange_grid),
                                 v  = Field{Center, Center, Nothing}(exchange_grid),
                                 T  = Field{Center, Center, Nothing}(exchange_grid),
                                 q  = Field{Center, Center, Nothing}(exchange_grid),
                                 p  = Field{Center, Center, Nothing}(exchange_grid),
                                 Qs = Field{Center, Center, Nothing}(exchange_grid),
                                 Qℓ = Field{Center, Center, Nothing}(exchange_grid),
                                 Mp = Field{Center, Center, Nothing}(exchange_grid))

    exchanger = atmosphere_exchanger(atmosphere, exchange_grid)

    return StateExchanger(ocean.model.grid, exchange_atmosphere_state, exchanger)
end

function atmosphere_exchanger(atmosphere::PrescribedAtmosphere, exchange_grid)
    atmos_grid = atmosphere.grid
    arch = architecture(exchange_grid)
    Nx, Ny, Nz = size(exchange_grid)

    # Make a NamedTuple of fractional indices
    # Note: we could use an array of FractionalIndices. Instead, for compatbility
    # with Reactant we construct FractionalIndices on the fly in `interpolate_atmospheric_state`.
    FT = eltype(atmos_grid)
    TX, TY, TZ = topology(exchange_grid)
    fi = TX() isa Flat ? nothing : Field{Center, Center, Nothing}(exchange_grid, FT)
    fj = TY() isa Flat ? nothing : Field{Center, Center, Nothing}(exchange_grid, FT)
    frac_indices = (i=fi, j=fj) # no k needed, only horizontal interpolation

    kernel_parameters = interface_kernel_parameters(exchange_grid)
    launch!(arch, exchange_grid, kernel_parameters,
            _compute_fractional_indices!, frac_indices, exchange_grid, atmos_grid)

    return frac_indices
end

@kernel function _compute_fractional_indices!(indices_tuple, exchange_grid, atmos_grid)
    i, j = @index(Global, NTuple)
    kᴺ = size(exchange_grid, 3) # index of the top ocean cell
    X = _node(i, j, kᴺ + 1, exchange_grid, c, c, f)
    fractional_indices_ij = FractionalIndices(X, atmos_grid, c, c, nothing)
    fi = indices_tuple.i
    fj = indices_tuple.j
    @inbounds begin
        if !isnothing(fi)
            fi[i, j, 1] = fractional_indices_ij.i
        end

        if !isnothing(fj)
            fj[i, j, 1] = fractional_indices_ij.j
        end
    end
end

const PATP = PrescribedAtmosphereThermodynamicsParameters

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

function atmosphere_ocean_interface(ocean,
                                    radiation,
                                    ao_flux_formulation,
                                    temperature_formulation, 
                                    velocity_formulation,
                                    specific_humidity_formulation)

    water_vapor   = Field{Center, Center, Nothing}(ocean.model.grid)
    latent_heat   = Field{Center, Center, Nothing}(ocean.model.grid)
    sensible_heat = Field{Center, Center, Nothing}(ocean.model.grid)
    x_momentum    = Field{Center, Center, Nothing}(ocean.model.grid)
    y_momentum    = Field{Center, Center, Nothing}(ocean.model.grid)
    ao_fluxes = (; latent_heat, sensible_heat, water_vapor, x_momentum, y_momentum)

    σ = radiation.stefan_boltzmann_constant
    αₐₒ = radiation.reflection.ocean
    ϵₐₒ = radiation.emission.ocean
    radiation = (σ=σ, α=αₐₒ, ϵ=ϵₐₒ)

    ao_properties = InterfaceProperties(radiation,
                                        specific_humidity_formulation,
                                        temperature_formulation,
                                        velocity_formulation)

    interface_temperature = Field{Center, Center, Nothing}(ocean.model.grid)

    return AtmosphereInterface(ao_fluxes, ao_flux_formulation, interface_temperature, ao_properties)
end

atmosphere_sea_ice_interface(sea_ice, args...) = nothing

function atmosphere_sea_ice_interface(sea_ice::SeaIceSimulation,
                                      radiation,
                                      ai_flux_formulation,
                                      temperature_formulation,
                                      velocity_formulation)

    water_vapor   = Field{Center, Center, Nothing}(sea_ice.model.grid)
    latent_heat   = Field{Center, Center, Nothing}(sea_ice.model.grid)
    sensible_heat = Field{Center, Center, Nothing}(sea_ice.model.grid)
    x_momentum    = Field{Center, Center, Nothing}(sea_ice.model.grid)
    y_momentum    = Field{Center, Center, Nothing}(sea_ice.model.grid)
    fluxes = (; latent_heat, sensible_heat, water_vapor, x_momentum, y_momentum)

    σ = radiation.stefan_boltzmann_constant
    αₐᵢ = radiation.reflection.sea_ice
    ϵₐᵢ = radiation.emission.sea_ice
    radiation = (σ=σ, α=αₐᵢ, ϵ=ϵₐᵢ)

    phase = AtmosphericThermodynamics.Ice()
    specific_humidity_formulation = SpecificHumidityFormulation(phase)

    properties = InterfaceProperties(radiation,
                                     specific_humidity_formulation,
                                     temperature_formulation,
                                     velocity_formulation)

    interface_temperature = sea_ice.model.ice_thermodynamics.top_surface_temperature

    return AtmosphereInterface(fluxes, ai_flux_formulation, interface_temperature, properties)
end

sea_ice_ocean_interface(sea_ice, ocean) = nothing

function sea_ice_ocean_interface(sea_ice::SeaIceSimulation, ocean;
                                 characteristic_melting_speed = 1e-5)

    previous_ice_thickness = deepcopy(sea_ice.model.ice_thickness)
    previous_ice_concentration = deepcopy(sea_ice.model.ice_concentration)
    io_bottom_heat_flux = Field{Center, Center, Nothing}(ocean.model.grid)
    io_frazil_heat_flux = Field{Center, Center, Nothing}(ocean.model.grid)
    io_salt_flux = Field{Center, Center, Nothing}(ocean.model.grid)

    @assert io_frazil_heat_flux isa Field{Center, Center, Nothing}
    @assert io_bottom_heat_flux isa Field{Center, Center, Nothing}
    @assert io_salt_flux isa Field{Center, Center, Nothing}

    io_fluxes = (interface_heat=io_bottom_heat_flux,
                 frazil_heat=io_frazil_heat_flux,
                 salt=io_salt_flux)

    io_properties = (; characteristic_melting_speed)

    return SeaIceOceanInterface(io_fluxes,
                                io_properties,
                                previous_ice_thickness,
                                previous_ice_concentration)
end

default_ai_temperature(sea_ice) = nothing

function default_ai_temperature(sea_ice::SeaIceSimulation)
    conductive_flux = sea_ice.model.ice_thermodynamics.internal_heat_flux.parameters.flux
    return SkinTemperature(conductive_flux)
end

function default_ao_specific_humidity(ocean)
    FT = eltype(ocean.model.grid)
    phase = AtmosphericThermodynamics.Liquid()
    x_H₂O = convert(FT, 0.98)
    return SpecificHumidityFormulation(phase, x_H₂O)
end

"""
    ComponentInterfaces(atmosphere, ocean, sea_ice=nothing;
                        radiation = Radiation(),
                        freshwater_density = 1000,
                        atmosphere_ocean_flux_formulation = SimilarityTheoryFluxes(),
                        atmosphere_sea_ice_flux_formulation = CoefficientBasedFluxes(drag_coefficient=2e-3,
                                                                                     heat_transfer_coefficient=1e-4,
                                                                                     vapor_flux_coefficient=1e-4),
                        atmosphere_ocean_interface_temperature = BulkTemperature(),
                        atmosphere_ocean_interface_specific_humidity = default_ao_specific_humidity(ocean),
                        atmosphere_sea_ice_interface_temperature = default_ai_temperature(sea_ice),
                        ocean_reference_density = reference_density(ocean),
                        ocean_heat_capacity = heat_capacity(ocean),
                        ocean_temperature_units = DegreesCelsius(),
                        sea_ice_temperature_units = DegreesCelsius(),
                        sea_ice_reference_density = reference_density(sea_ice),
                        sea_ice_heat_capacity = heat_capacity(sea_ice))
"""
function ComponentInterfaces(atmosphere, ocean, sea_ice=nothing;
                             radiation = Radiation(),
                             freshwater_density = 1000,
                             atmosphere_ocean_flux_formulation = SimilarityTheoryFluxes(eltype(ocean.model.grid)),
                             atmosphere_sea_ice_flux_formulation = SimilarityTheoryFluxes(eltype(ocean.model.grid)),
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
                             sea_ice_heat_capacity = heat_capacity(sea_ice))

    ocean_grid = ocean.model.grid
    FT = eltype(ocean_grid)

    ocean_reference_density   = convert(FT, ocean_reference_density)
    ocean_heat_capacity       = convert(FT, ocean_heat_capacity)
    sea_ice_reference_density = convert(FT, sea_ice_reference_density)
    sea_ice_heat_capacity     = convert(FT, sea_ice_heat_capacity)
    freshwater_density        = convert(FT, freshwater_density)

    atmosphere_properties = thermodynamics_parameters(atmosphere)

    ocean_properties = (reference_density  = ocean_reference_density,
                        heat_capacity      = ocean_heat_capacity,
                        freshwater_density = freshwater_density,
                        temperature_units  = ocean_temperature_units)

    ao_interface = atmosphere_ocean_interface(ocean,
                                              radiation,
                                              atmosphere_ocean_flux_formulation,
                                              atmosphere_ocean_interface_temperature,
                                              atmosphere_ocean_velocity_difference,
                                              atmosphere_ocean_interface_specific_humidity)

    io_interface = sea_ice_ocean_interface(sea_ice, ocean)

    ai_interface = atmosphere_sea_ice_interface(sea_ice,
                                                radiation,
                                                atmosphere_sea_ice_flux_formulation,
                                                atmosphere_sea_ice_interface_temperature,
                                                atmosphere_sea_ice_velocity_difference)

    if sea_ice isa SeaIceSimulation
        sea_ice_properties = (reference_density  = sea_ice_reference_density,
                              heat_capacity      = sea_ice_heat_capacity,
                              freshwater_density = freshwater_density,
                              liquidus           = sea_ice.model.ice_thermodynamics.phase_transitions.liquidus,
                              temperature_units  = sea_ice_temperature_units)

        net_momentum_fluxes = if sea_ice.model.dynamics isa Nothing 
            u = Field{Face, Center, Nothing}(sea_ice.model.grid)
            v = Field{Center, Face, Nothing}(sea_ice.model.grid)
            (; u, v) 
        else
            u = sea_ice.model.dynamics.external_momentum_stresses.top.u
            v = sea_ice.model.dynamics.external_momentum_stresses.top.v
            (; u, v)
        end

        net_top_sea_ice_fluxes = merge((; heat=sea_ice.model.external_heat_fluxes.top), net_momentum_fluxes)
        net_bottom_sea_ice_fluxes = (; heat=sea_ice.model.external_heat_fluxes.bottom)
    else
        sea_ice_properties = nothing
        net_top_sea_ice_fluxes = nothing
        net_bottom_sea_ice_fluxes = nothing
    end

    τx = surface_flux(ocean.model.velocities.u)
    τy = surface_flux(ocean.model.velocities.v)
    tracers = ocean.model.tracers
    ρₒ = ocean_reference_density
    cₒ = ocean_heat_capacity
    Qₒ = ρₒ * cₒ * surface_flux(ocean.model.tracers.T)
    net_ocean_surface_fluxes = (u=τx, v=τy, Q=Qₒ)

    ocean_surface_tracer_fluxes = NamedTuple(name => surface_flux(tracers[name]) for name in keys(tracers))
    net_ocean_surface_fluxes = merge(ocean_surface_tracer_fluxes, net_ocean_surface_fluxes)

    # Total interface fluxes
    net_fluxes = (ocean_surface  = net_ocean_surface_fluxes,
                  sea_ice_top    = net_top_sea_ice_fluxes,
                  sea_ice_bottom = net_bottom_sea_ice_fluxes)

    exchanger = StateExchanger(ocean, atmosphere)

    return ComponentInterfaces(ao_interface,
                               ai_interface,
                               io_interface,
                               atmosphere_properties,
                               ocean_properties,
                               sea_ice_properties,
                               exchanger,
                               net_fluxes)
end

sea_ice_similarity_theory(sea_ice) = nothing

function sea_ice_similarity_theory(sea_ice::SeaIceSimulation)
    # Here we need to make sure the interface temperature type is
    # SkinTemperature. Also we need to pass the sea ice internal flux
    # The thickness and salinity need to be passed as well,
    # but the can be passed as state variables once we refactor the `StateValues` struct.
    internal_flux = sea_ice.model.ice_thermodynamics.internal_heat_flux
    interface_temperature_type = SkinTemperature(internal_flux)
    return SimilarityTheoryFluxes(; interface_temperature_type)
end

