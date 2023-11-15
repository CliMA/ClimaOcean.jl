using Oceananigans.Models: update_model_field_time_series!
using Oceananigans.TimeSteppers: Clock
using Oceananigans

struct Radiation{FT, SW, LW}
    downwelling_shortwave_radiation :: SW
    downwelling_longwave_radiation :: LW
    ocean_emissivity :: FT
    stefan_boltzmann_constant :: FT
    reference_temperature :: FT
end

struct OceanSeaIceModel{FT, I, A, O, C, AO, AI, IO, G, R, PI, PC} <: AbstractModel{Nothing}
    clock :: C
    grid :: G # TODO: make it so simulation does not require this
    atmosphere :: A
    sea_ice :: I
    ocean :: O
    atmosphere_ocean_fluxes :: AO
    atmosphere_sea_ice_fluxes :: AI
    sea_ice_ocean_fluxes :: IO
    radiation :: R
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

function OceanSeaIceModel(ice, ocean, atmosphere=nothing; 
                          clock = default_clock(eltype(ocean.model)))
    
    previous_ice_thickness = deepcopy(ice.model.ice_thickness)
    previous_ice_concentration = deepcopy(ice.model.ice_concentration)

    grid = ocean.model.grid
    ice_ocean_heat_flux = Field{Center, Center, Nothing}(grid)
    ice_ocean_salt_flux = Field{Center, Center, Nothing}(grid)

    ocean_reference_density = 1024
    ocean_heat_capacity = 3991
    reference_temperature = 273.15

    # How would we ensure consistency?
    try
        if ice.model.external_heat_fluxes.top isa RadiativeEmission
            ice_radiation = ice.model.external_heat_fluxes.top
        else
            ice_radiation = filter(flux isa RadiativeEmission, ice.model.external_heat_fluxes.top) |> first
        end

        reference_temperature = ice_radiation.reference_temperature
    catch
    end

    FT = eltype(ocean.model.grid)
    ocean_emissivity = 1
    stefan_boltzmann_constant = 5.67e-8

    radiation = Radiation(nothing, nothing,
                          convert(FT, ocean_emissivity),
                          convert(FT, stefan_boltzmann_constant),
                          convert(FT, reference_temperature))

    # Set up fluxes
    momentum_fluxes = (u = surface_flux(ocean.model.velocities.u),
                       v = surface_flux(ocean.model.velocities.v))

    @show momentum_fluxes
    @show momentum_fluxes.u
    @show momentum_fluxes.v

    atmosphere_ocean_fluxes = AtmosphereOceanMomentumFlux()
    atmosphere_sea_ice_fluxes = nothing
    sea_ice_ocean_fluxes = nothing

    return OceanSeaIceModel(clock,
                            ocean.model.grid,
                            atmosphere,
                            ice,
                            ocean,
                            atmosphere_ocean_fluxes,
                            atmosphere_sea_ice_fluxes,
                            sea_ice_ocean_fluxes,
                            radiation,
                            previous_ice_thickness,
                            previous_ice_concentration,
                            convert(FT, ocean_reference_density),
                            convert(FT, ocean_heat_capacity))
end

time(coupled_model::OceanSeaIceModel) = coupled_model.clock.time

function time_step!(coupled_model::OceanSeaIceModel, Δt; callbacks=nothing)
    ocean = coupled_model.ocean
    sea_ice = coupled_model.sea_ice

    h = sea_ice.model.ice_thickness
    fill_halo_regions!(h)

    # Initialization
    if coupled_model.clock.iteration == 0
        @info "Initializing coupled model ice thickness..."
        h⁻ = coupled_model.previous_ice_thickness
        hⁿ = coupled_model.sea_ice.model.ice_thickness
        parent(h⁻) .= parent(hⁿ)
    end

    sea_ice.Δt   = Δt
    ocean.Δt = Δt

    time_step!(sea_ice)

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
    # update_model_field_time_series!(coupled_model.atmosphere.model) 
    compute_atmosphere_ocean_fluxes!(coupled_model) 
    # compute_atmosphere_sea_ice_fluxes!(coupled_model)
    compute_sea_ice_ocean_fluxes!(coupled_model)
    return nothing
end

