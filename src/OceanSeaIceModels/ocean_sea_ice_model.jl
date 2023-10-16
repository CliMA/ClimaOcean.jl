using Oceananigans.Models: update_model_field_time_series!

struct OceanSeaIceModel{FT, I, O, F, C, G, S, PI, PC} <: AbstractModel{Nothing}
    clock :: C
    grid :: G # TODO: make it so simulation does not require this
    ice :: I
    previous_ice_thickness :: PI
    previous_ice_concentration :: PC
    ocean :: O
    atmospheric_forcing :: F
    solar_insolation :: S
    ocean_density :: FT
    ocean_heat_capacity :: FT
    ocean_emissivity :: FT
    reference_temperature :: FT
end

const OSIM = OceanSeaIceModel

Base.summary(::OSIM) = "OceanSeaIceModel"
prettytime(model::OSIM) = prettytime(model.clock.time)
iteration(model::OSIM) = model.clock.iteration
timestepper(::OSIM) = nothing
reset!(::OSIM) = nothing
initialize!(::OSIM) = nothing
default_included_properties(::OSIM) = tuple()
update_state!(::OSIM) = nothing
prognostic_fields(cm::OSIM) = nothing
fields(::OSIM) = NamedTuple()
default_clock(TT) = Clock{TT}(0, 0, 1)

# "Ocean only"
const OceanOnlyModel    = OceanSeaIceModel{<:Any, Nothing}
const NoAtmosphereModel = OceanSeaIceModel{<:Any, <:Any, Nothing}

OceanOnlyModel(ocean; atmospheric_forcing = nothing, clock = default_clock(eltype(ocean.model))) = 
    OceanSeaIceModel(nothing, ocean; atmospheric_forcing, clock)

function OceanSeaIceModel(ice, ocean; 
                          atmospheric_forcing = nothing,
                          clock = default_clock(eltype(ocean.model)))
    
    previous_ice_thickness = deepcopy(ice.model.ice_thickness)
    previous_ice_concentration = deepcopy(ice.model.ice_concentration)

    grid = ocean.model.grid
    ice_ocean_thermal_flux = Field{Center, Center, Nothing}(grid)
    ice_ocean_salt_flux = Field{Center, Center, Nothing}(grid)
    solar_insolation = Field{Center, Center, Nothing}(grid)

    ocean_density = 1024
    ocean_heat_capacity = 3991
    ocean_emissivity = 1
    reference_temperature = 273.15

    # How would we ensure consistency?
    try
        if ice.model.external_thermal_fluxes.top isa RadiativeEmission
            radiation = ice.model.external_thermal_fluxes.top
        else
            radiation = filter(flux isa RadiativeEmission, ice.model.external_thermal_fluxes.top) |> first
        end

        reference_temperature = radiation.reference_temperature
    catch
    end

    FT = eltype(ocean.model.grid)

    return OceanSeaIceModel(clock,
                            ocean.model.grid,
                            ice,
                            previous_ice_thickness,
                            previous_ice_concentration,
                            ocean,
                            solar_insolation,
                            atmospheric_forcing,
                            convert(FT, ocean_density),
                            convert(FT, ocean_heat_capacity),
                            convert(FT, ocean_emissivity),
                            convert(FT, stefan_boltzmann_constant),
                            convert(FT, reference_temperature))
end

time(coupled_model::OceanSeaIceModel) = coupled_model.clock.time

function time_step!(coupled_model::OceanSeaIceModel, Δt; callbacks=nothing)
    ocean = coupled_model.ocean
    ice = coupled_model.ice
    ice.Δt   = Δt
    ocean.Δt = Δt

    fill_halo_regions!(h)

    # Initialization
    if coupled_model.clock.iteration == 0
        h⁻ = coupled_model.previous_ice_thickness
        hⁿ = coupled_model.ice.model.ice_thickness
        parent(h⁻) .= parent(hⁿ)
    end

    time_step!(ice)

    # TODO: put this in update_state!
    # Air-sea and Air-ice fluxes substitute the previous values
    # while ice-ocean fluxes are additive
    update_model_field_time_series!(coupled_model.atmospheric_forcing) 
    compute_air_sea_flux!(coupled_model) 
    compute_air_ice_flux!(coupled_model) # TODO: we need to implement this, not sure how
    compute_ice_ocean_flux!(coupled_model)

    time_step!(ocean)

    # TODO:
    # - Store fractional ice-free / ice-covered _time_ for more
    #   accurate flux computation?
    # - Or, input "excess heat flux" into ocean after the ice melts
    # - Currently, non-conservative for heat due bc we don't account for excess
        
    # TODO after ice time-step:
    #   - Adjust ocean temperature if the ice completely melts?
   
    tick!(coupled_model.clock, Δt)
    
    return nothing
end
