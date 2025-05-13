using Printf
using Oceananigans.Operators: ℑxᶠᵃᵃ, ℑyᵃᶠᵃ
using Oceananigans.Forcings: MultipleForcings

using ClimaOcean.OceanSeaIceModels: sea_ice_concentration

@inline computed_sea_ice_ocean_fluxes(interface) = interface.fluxes
@inline computed_sea_ice_ocean_fluxes(::Nothing) = (interface_heat = ZeroField(),
                                                    frazil_heat = ZeroField(),
                                                    salt = ZeroField(),
                                                    x_momentum = ZeroField(),
                                                    y_momentum = ZeroField())

@inline shortwave_radiative_forcing(i, j, grid, Fᵀ, Qts, ocean_properties) = Qts

@inline function shortwave_radiative_forcing(i, j, grid, tcr::TwoColorRadiation, Iˢʷ, ocean_properties)
    ρₒ = ocean_properties.reference_density
    cₒ = ocean_properties.heat_capacity
    J₀ = tcr.surface_flux
    @inbounds J₀[i, j,  1] = - Iˢʷ / (ρₒ * cₒ)
    return zero(Iˢʷ)
end

get_radiative_forcing(FT) = FT
function get_radiative_forcing(FT::MultipleForcings)
    for forcing in FT.forcings
        forcing isa TwoColorRadiation && return forcing
    end
    return nothing
end

function compute_net_ocean_fluxes!(coupled_model)
    ocean = coupled_model.ocean
    sea_ice = coupled_model.sea_ice
    grid = ocean.model.grid
    arch = architecture(grid)
    clock = coupled_model.clock

    net_ocean_fluxes = coupled_model.interfaces.net_fluxes.ocean_surface
    atmos_ocean_fluxes = coupled_model.interfaces.atmosphere_ocean_interface.fluxes
    sea_ice_ocean_fluxes = computed_sea_ice_ocean_fluxes(coupled_model.interfaces.sea_ice_ocean_interface)

    # We remove the heat flux since does not need to be assembled and bloats the parameter space.
    net_ocean_fluxes = (u = net_ocean_fluxes.u,
                        v = net_ocean_fluxes.v,
                        T = net_ocean_fluxes.T,
                        S = net_ocean_fluxes.S)

    # Simplify NamedTuple to reduce parameter space consumption.
    # See https://github.com/CliMA/ClimaOcean.jl/issues/116.
    atmosphere_fields = coupled_model.interfaces.exchanger.exchange_atmosphere_state

    downwelling_radiation = (Qs = atmosphere_fields.Qs.data,
                             Qℓ = atmosphere_fields.Qℓ.data)

    freshwater_flux = atmosphere_fields.Mp.data

    ice_concentration = sea_ice_concentration(sea_ice)
    ocean_salinity = ocean.model.tracers.S
    atmos_ocean_properties = coupled_model.interfaces.atmosphere_ocean_interface.properties
    ocean_properties = coupled_model.interfaces.ocean_properties
    kernel_parameters = interface_kernel_parameters(grid)

    ocean_surface_temperature = coupled_model.interfaces.atmosphere_ocean_interface.temperature
    penetrating_radiation = get_radiative_forcing(ocean.model.forcing.T)

    launch!(arch, grid, kernel_parameters,
            _assemble_net_ocean_fluxes!,
            net_ocean_fluxes,
            penetrating_radiation,
            grid,
            clock,
            atmos_ocean_fluxes,
            sea_ice_ocean_fluxes,
            ocean_salinity,
            ocean_surface_temperature,
            ice_concentration,
            downwelling_radiation,
            freshwater_flux,
            atmos_ocean_properties,
            ocean_properties)

    return nothing
end

@inline τᶜᶜᶜ(i, j, k, grid, ρₒ⁻¹, ℵ, ρτᶜᶜᶜ) = @inbounds ρₒ⁻¹ * (1 - ℵ[i, j, k]) * ρτᶜᶜᶜ[i, j, k]

@kernel function _assemble_net_ocean_fluxes!(net_ocean_fluxes,
                                             penetrating_radiation,
                                             grid,
                                             clock,
                                             atmos_ocean_fluxes,
                                             sea_ice_ocean_fluxes,
                                             ocean_salinity,
                                             ocean_surface_temperature,
                                             sea_ice_concentration,
                                             downwelling_radiation,
                                             freshwater_flux,
                                             atmos_ocean_properties,
                                             ocean_properties)

    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3)
    time = Time(clock.time)
    ρτxao = atmos_ocean_fluxes.x_momentum # atmosphere - ocean zonal momentum flux
    ρτyao = atmos_ocean_fluxes.y_momentum # atmosphere - ocean meridional momentum flux
    ρτxio = sea_ice_ocean_fluxes.x_momentum # sea_ice - ocean zonal momentum flux
    ρτyio = sea_ice_ocean_fluxes.y_momentum # sea_ice - ocean meridional momentum flux

    @inbounds begin
        Sₒ = ocean_salinity[i, j, kᴺ]
        Tₛ = ocean_surface_temperature[i, j, 1]
        Tₛ = convert_to_kelvin(ocean_properties.temperature_units, Tₛ)

        Mp  = freshwater_flux[i, j, 1] # Prescribed freshwater flux
        Qs  = downwelling_radiation.Qs[i, j, 1] # Downwelling shortwave radiation
        Qℓ  = downwelling_radiation.Qℓ[i, j, 1] # Downwelling longwave radiation
        Qc  = atmos_ocean_fluxes.sensible_heat[i, j, 1] # sensible or "conductive" heat flux
        Qv  = atmos_ocean_fluxes.latent_heat[i, j, 1]   # latent heat flux
        Mv  = atmos_ocean_fluxes.water_vapor[i, j, 1]   # mass flux of water vapor
    end

    # Compute radiation fluxes
    σ = atmos_ocean_properties.radiation.σ
    α = atmos_ocean_properties.radiation.α
    ϵ = atmos_ocean_properties.radiation.ϵ
    Qu = emitted_longwave_radiation(i, j, kᴺ, grid, time, Tₛ, σ, ϵ) 
    Qaℓ = absorbed_longwave_radiation(i, j, kᴺ, grid, time, ϵ, Qℓ)

    # Compute the interior + surface absorbed shortwave radiation
    Qts = transmitted_shortwave_radiation(i, j, kᴺ, grid, time, α, Qs)
    Qss = shortwave_radiative_forcing(i, j, grid, penetrating_radiation, Qts, ocean_properties)

    # Compute the total radiation
    ΣQao = Qu + Qc + Qv + Qaℓ + Qss

    @inbounds begin
        # Write radiative components of the heat flux for diagnostic purposes
        atmos_ocean_fluxes.upwelling_longwave[i, j, 1] = Qu
        atmos_ocean_fluxes.downwelling_longwave[i, j, 1] = - Qaℓ
        atmos_ocean_fluxes.downwelling_shortwave[i, j, 1] = - Qts
    end

    # Convert from a mass flux to a volume flux (aka velocity)
    # by dividing with the density of freshwater.
    # Also switch the sign, for some reason we are given freshwater flux as positive down.
    ρf⁻¹ = 1 / ocean_properties.freshwater_density
    ΣFao = - Mp * ρf⁻¹

    # Add the contribution from the turbulent water vapor flux, which has
    # a different sign convention as the prescribed water mass fluxes (positive upwards)
    Fv = Mv * ρf⁻¹
    ΣFao += Fv

    # Compute fluxes for u, v, T, S from momentum, heat, and freshwater fluxes
    τx = net_ocean_fluxes.u
    τy = net_ocean_fluxes.v
    Jᵀ = net_ocean_fluxes.T
    Jˢ = net_ocean_fluxes.S
    ℵ = sea_ice_concentration
    ρₒ⁻¹ = 1 / ocean_properties.reference_density
    cₒ   = ocean_properties.heat_capacity

    @inbounds begin
        ℵᵢ   = ℵ[i, j, 1]
        Qio  = sea_ice_ocean_fluxes.interface_heat[i, j, 1]

        Jᵀao = ΣQao  * ρₒ⁻¹ / cₒ
        Jˢao = - Sₒ * ΣFao
        Jᵀio = Qio * ρₒ⁻¹ / cₒ
        Jˢio = sea_ice_ocean_fluxes.salt[i, j, 1] * ℵᵢ

        τxao = ℑxᶠᵃᵃ(i, j, 1, grid, τᶜᶜᶜ, ρₒ⁻¹, ℵ, ρτxao)
        τyao = ℑyᵃᶠᵃ(i, j, 1, grid, τᶜᶜᶜ, ρₒ⁻¹, ℵ, ρτyao)
        τxio = ρτxio[i, j, 1] * ρₒ⁻¹ * ℑxᶠᵃᵃ(i, j, 1, grid, ℵ)
        τyio = ρτyio[i, j, 1] * ρₒ⁻¹ * ℑyᵃᶠᵃ(i, j, 1, grid, ℵ)

        # Stresses
        τx[i, j, 1] = τxao + τxio
        τy[i, j, 1] = τyao + τyio

        # Tracer fluxes
        Jᵀ[i, j, 1] = (1 - ℵᵢ) * Jᵀao + Jᵀio
        Jˢ[i, j, 1] = (1 - ℵᵢ) * Jˢao + Jˢio
    end
end

function compute_net_sea_ice_fluxes!(coupled_model)
    sea_ice = coupled_model.sea_ice

    if !(sea_ice isa SeaIceSimulation)
        return nothing
    end

    ocean = coupled_model.ocean
    grid  = ocean.model.grid
    arch  = architecture(grid)
    clock = coupled_model.clock

    top_fluxes = coupled_model.interfaces.net_fluxes.sea_ice_top
    bottom_heat_flux = coupled_model.interfaces.net_fluxes.sea_ice_bottom.heat
    sea_ice_ocean_fluxes = coupled_model.interfaces.sea_ice_ocean_interface.fluxes
    atmosphere_sea_ice_fluxes = coupled_model.interfaces.atmosphere_sea_ice_interface.fluxes

    # Simplify NamedTuple to reduce parameter space consumption.
    # See https://github.com/CliMA/ClimaOcean.jl/issues/116.
    atmosphere_fields = coupled_model.interfaces.exchanger.exchange_atmosphere_state

    downwelling_radiation = (Qs = atmosphere_fields.Qs.data,
                             Qℓ = atmosphere_fields.Qℓ.data)

    freshwater_flux = atmosphere_fields.Mp.data

    atmos_sea_ice_properties = coupled_model.interfaces.atmosphere_sea_ice_interface.properties
    sea_ice_properties = coupled_model.interfaces.sea_ice_properties

    kernel_parameters = interface_kernel_parameters(grid)

    sea_ice_surface_temperature = coupled_model.interfaces.atmosphere_sea_ice_interface.temperature
    ice_concentration = sea_ice_concentration(sea_ice)
    
    launch!(arch, grid, kernel_parameters,
            _assemble_net_sea_ice_fluxes!,
            top_fluxes,
            bottom_heat_flux,
            grid,
            clock,
            atmosphere_sea_ice_fluxes,
            sea_ice_ocean_fluxes,
            freshwater_flux,
            ice_concentration,
            sea_ice_surface_temperature,
            downwelling_radiation,
            sea_ice_properties,
            atmos_sea_ice_properties)

    return nothing
end

@kernel function _assemble_net_sea_ice_fluxes!(top_fluxes,
                                               bottom_heat_flux,
                                               grid,
                                               clock,
                                               atmosphere_sea_ice_fluxes,
                                               sea_ice_ocean_fluxes,
                                               freshwater_flux, # Where do we add this one?
                                               ice_concentration,
                                               surface_temperature,
                                               downwelling_radiation,
                                               sea_ice_properties,
                                               atmos_sea_ice_properties)

    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3)
    time = Time(clock.time)

    @inbounds begin
        Ts = surface_temperature[i, j, kᴺ]
        Ts = convert_to_kelvin(sea_ice_properties.temperature_units, Ts)
        ℵi = ice_concentration[i, j, 1]
        
        Qs = downwelling_radiation.Qs[i, j, 1]
        Qℓ = downwelling_radiation.Qℓ[i, j, 1]
        Qc = atmosphere_sea_ice_fluxes.sensible_heat[i, j, 1] # sensible or "conductive" heat flux
        Qv = atmosphere_sea_ice_fluxes.latent_heat[i, j, 1]   # latent heat flux
        Qf = sea_ice_ocean_fluxes.frazil_heat[i, j, 1]        # frazil heat flux
        Qi = sea_ice_ocean_fluxes.interface_heat[i, j, 1]   # interfacial heat flux
    end

    ρτx = atmosphere_sea_ice_fluxes.x_momentum  # zonal momentum flux
    ρτy = atmosphere_sea_ice_fluxes.y_momentum  # meridional momentum flux

    # Compute radiation fluxes
    σ = atmos_sea_ice_properties.radiation.σ
    α = atmos_sea_ice_properties.radiation.α
    ϵ = atmos_sea_ice_properties.radiation.ϵ
    Qu = emitted_longwave_radiation(i, j, kᴺ, grid, time, Ts, σ, ϵ) 
    Qd = net_absorbed_interface_radiation(i, j, kᴺ, grid, time, α, ϵ, Qs, Qℓ)

    ΣQt = (Qd + Qu + Qc + Qv) * ℵi # If ℵi == 0 there is no heat flux from the top!
    ΣQb = Qf + Qi

    # Mask fluxes over land for convenience
    inactive = inactive_node(i, j, kᴺ, grid, c, c, c)

    @inbounds top_fluxes.heat[i, j, 1]  = ifelse(inactive, zero(grid), ΣQt)
    @inbounds top_fluxes.u[i, j, 1]     = ifelse(inactive, zero(grid), ℑxᶠᵃᵃ(i, j, 1, grid, ρτx))
    @inbounds top_fluxes.v[i, j, 1]     = ifelse(inactive, zero(grid), ℑyᵃᶠᵃ(i, j, 1, grid, ρτy))
    @inbounds bottom_heat_flux[i, j, 1] = ifelse(inactive, zero(grid), ΣQb)
end
