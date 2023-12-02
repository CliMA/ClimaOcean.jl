using Oceananigans
using Oceananigans.Models: update_model_field_time_series!
using Oceananigans.TimeSteppers: Clock
using Oceananigans.BuoyancyModels: SeawaterBuoyancy

using SeawaterPolynomials: TEOS10EquationOfState

struct OceanSeaIceModel{FT, I, A, O, S, F, PI, PC, C, G} <: AbstractModel{Nothing}
    clock :: C
    grid :: G # TODO: make it so Oceananigans.Ssimulation does not require this
    atmosphere :: A
    sea_ice :: I
    ocean :: O
    surfaces :: S
    fluxes :: F
    previous_ice_thickness :: PI
    previous_ice_concentration :: PC
    # The ocean is Boussinesq, so these are _only_ coupled properties:
    ocean_reference_density :: FT
    ocean_heat_capacity :: FT
end

const OSIM = OceanSeaIceModel

Base.summary(::OSIM) = "OceanSeaIceModel"
prettytime(model::OSIM) = prettytime(model.clock.time)
iteration(model::OSIM) = model.clock.iteration
timestepper(::OSIM) = nothing
reset!(::OSIM) = nothing
initialize!(::OSIM) = nothing
default_included_properties(::OSIM) = tuple()
prognostic_fields(cm::OSIM) = nothing
fields(::OSIM) = NamedTuple()
default_clock(TT) = Oceananigans.TimeSteppers.Clock{TT}(0, 0, 1)

reference_density(ocean::Simulation) = reference_density(ocean.model.buoyancy.model)
reference_density(buoyancy_model::SeawaterBuoyancy) = reference_density(buoyancy_model.equation_of_state)
#reference_density(unsupported) = throw(ArgumentError("Cannot extract reference density from $(typeof(unsupported))"))
#reference_density(eos::TEOS10EquationOfState) = eos.reference_density
reference_density(eos) = eos.reference_density

heat_capacity(ocean::Simulation) = heat_capacity(ocean.model.buoyancy.model)
heat_capacity(buoyancy_model::SeawaterBuoyancy) = heat_capacity(buoyancy_model.equation_of_state)
#heat_capacity(unsupported) = throw(ArgumentError("Cannot deduce the heat capacity from $(typeof(unsupported))"))
#heat_capacity(eos::TEOS10EquationOfState) = 3991 # get the right value in here eventually
heat_capacity(eos) = 3991 # get the right value in here eventually

function OceanSeaIceModel(ocean, sea_ice=nothing;
                          atmosphere = nothing,
                          surface_radiation = nothing,
                          ocean_reference_density = reference_density(ocean),
                          ocean_heat_capacity = heat_capacity(ocean),
                          clock = deepcopy(ocean.model.clock))
    
    previous_ice_thickness = deepcopy(sea_ice.model.ice_thickness)
    previous_ice_concentration = deepcopy(sea_ice.model.ice_concentration)

    grid = ocean.model.grid
    ice_ocean_heat_flux = Field{Center, Center, Nothing}(grid)
    ice_ocean_salt_flux = Field{Center, Center, Nothing}(grid)
    
    # Contains information about flux contributions: bulk formula, prescribed
    # fluxes, etc.
    fluxes = OceanSeaIceModelFluxes(eltype(grid); surface_radiation)

    # Contains a reference to the Fields holding net surface fluxes:
    # ocean top surface, and both top and bottom sea ice surfaces
    surfaces = OceanSeaIceSurfaces(ocean, sea_ice)

    FT = eltype(ocean.model.grid)

    return OceanSeaIceModel(clock,
                            ocean.model.grid,
                            atmosphere,
                            sea_ice,
                            ocean,
                            surfaces,
                            fluxes,
                            previous_ice_thickness,
                            previous_ice_concentration,
                            convert(FT, ocean_reference_density),
                            convert(FT, ocean_heat_capacity))
end

time(coupled_model::OceanSeaIceModel) = coupled_model.clock.time

function time_step!(coupled_model::OceanSeaIceModel, Δt; callbacks=nothing)
    ocean = coupled_model.ocean
    sea_ice = coupled_model.sea_ice

    # Eventually, split out into OceanOnlyModel
    if !isnothing(sea_ice)
        h = sea_ice.model.ice_thickness
        fill_halo_regions!(h)

        # Initialization
        if coupled_model.clock.iteration == 0
            @info "Initializing coupled model ice thickness..."
            h⁻ = coupled_model.previous_ice_thickness
            hⁿ = coupled_model.sea_ice.model.ice_thickness
            parent(h⁻) .= parent(hⁿ)
        end

        sea_ice.Δt = Δt
        time_step!(sea_ice)
    end

    ocean.Δt = Δt

    # TODO after ice time-step:
    #   - Adjust ocean heat flux if the ice completely melts?

    time_step!(ocean)

    # TODO:
    # - Store fractional ice-free / ice-covered _time_ for more
    #   accurate flux computation?
   
    tick!(coupled_model.clock, Δt)
    
    return nothing
end

function update_state!(coupled_model::OceanSeaIceModel, callbacks=nothing)
    # update_model_field_time_series!(coupled_model.atmosphere) 
    compute_atmosphere_ocean_fluxes!(coupled_model) 
    # compute_atmosphere_sea_ice_fluxes!(coupled_model)
    # compute_sea_ice_ocean_fluxes!(coupled_model)
    return nothing
end

