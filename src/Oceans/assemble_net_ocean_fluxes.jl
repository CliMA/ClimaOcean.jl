using Printf
using Oceananigans.Operators: ℑxᶠᵃᵃ, ℑyᵃᶠᵃ
using Oceananigans.Forcings: MultipleForcings
using ClimaOcean.OceanSeaIceModels: OceanSeaIceModel, NoOceanInterfaceModel, NoInterfaceModel

using ClimaOcean.OceanSeaIceModels.InterfaceComputations: interface_kernel_parameters, 
                                                          computed_fluxes, 
                                                          get_possibly_zero_flux, 
                                                          sea_ice_concentration,
                                                          convert_to_kelvin,
                                                          emitted_longwave_radiation,
                                                          absorbed_longwave_radiation,
                                                          transmitted_shortwave_radiation
                                                          

@inline τᶜᶜᶜ(i, j, k, grid, ρₒ⁻¹, ℵ, ρτᶜᶜᶜ) = @inbounds ρₒ⁻¹ * (1 - ℵ[i, j, k]) * ρτᶜᶜᶜ[i, j, k]

#####
##### Generic flux assembler
#####

# Fallback for an ocean-only model (it has no interfaces!)
update_net_fluxes!(coupled_model::Union{NoOceanInterfaceModel, NoInterfaceModel}, ocean::Simulation{<:HydrostaticFreeSurfaceModel}) = nothing

update_net_fluxes!(coupled_model, ocean::Simulation{<:HydrostaticFreeSurfaceModel}) = 
    update_net_ocean_fluxes!(coupled_model, ocean, ocean.model.grid)

# A generic ocean flux assembler for a coupled model with both an atmosphere and sea ice
function update_net_ocean_fluxes!(coupled_model, ocean_model, grid)
    sea_ice = coupled_model.sea_ice
    arch = architecture(grid)
    clock = coupled_model.clock

    net_ocean_fluxes = coupled_model.interfaces.net_fluxes.ocean
    atmos_ocean_fluxes = computed_fluxes(coupled_model.interfaces.atmosphere_ocean_interface)
    sea_ice_ocean_fluxes = computed_fluxes(coupled_model.interfaces.sea_ice_ocean_interface)

    # Simplify NamedTuple to reduce parameter space consumption.
    # See https://github.com/CliMA/ClimaOcean.jl/issues/116.
    atmosphere_fields = coupled_model.interfaces.exchanger.atmosphere.state

    downwelling_radiation = (Qs = atmosphere_fields.Qs.data,
                             Qℓ = atmosphere_fields.Qℓ.data)

    freshwater_flux = atmosphere_fields.Mp.data

    ice_concentration = sea_ice_concentration(sea_ice)
    ocean_salinity = OceanSeaIceModels.ocean_salinity(ocean_model)
    atmos_ocean_properties = coupled_model.interfaces.atmosphere_ocean_interface.properties
    ocean_properties = coupled_model.interfaces.ocean_properties
    kernel_parameters = interface_kernel_parameters(grid)

    ocean_surface_temperature = coupled_model.interfaces.atmosphere_ocean_interface.temperature
    penetrating_radiation = get_radiative_forcing(ocean_model)

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
    ρτxao = get_possibly_zero_flux(atmos_ocean_fluxes,   :x_momentum) # atmosphere - ocean zonal momentum flux
    ρτyao = get_possibly_zero_flux(atmos_ocean_fluxes,   :y_momentum) # atmosphere - ocean meridional momentum flux
    ρτxio = get_possibly_zero_flux(sea_ice_ocean_fluxes, :x_momentum) # sea_ice - ocean zonal momentum flux
    ρτyio = get_possibly_zero_flux(sea_ice_ocean_fluxes, :y_momentum) # sea_ice - ocean meridional momentum flux

    @inbounds begin
        ℵᵢ = sea_ice_concentration[i, j, 1]
        Sₒ = ocean_salinity[i, j, kᴺ]
        Tₛ = ocean_surface_temperature[i, j, 1]
        Tₛ = convert_to_kelvin(ocean_properties.temperature_units, Tₛ)

        Mp  = freshwater_flux[i, j, 1] # Prescribed freshwater flux
        Qs  = downwelling_radiation.Qs[i, j, 1] # Downwelling shortwave radiation
        Qℓ  = downwelling_radiation.Qℓ[i, j, 1] # Downwelling longwave radiation
        Qc  = get_possibly_zero_flux(atmos_ocean_fluxes, :sensible_heat)[i, j, 1] # sensible or "conductive" heat flux
        Qv  = get_possibly_zero_flux(atmos_ocean_fluxes, :latent_heat)[i, j, 1] # latent heat flux
        Mv  = get_possibly_zero_flux(atmos_ocean_fluxes, :water_vapor)[i, j, 1] # mass flux of water vapor
    end

    # Compute radiation fluxes (radiation is multiplied by the fraction of ocean, 1 - sea ice concentration)
    σ = atmos_ocean_properties.radiation.σ
    α = atmos_ocean_properties.radiation.α
    ϵ = atmos_ocean_properties.radiation.ϵ
    Qu  = emitted_longwave_radiation(i, j, kᴺ, grid, time, Tₛ, σ, ϵ) 
    Qaℓ = absorbed_longwave_radiation(i, j, kᴺ, grid, time, ϵ, Qℓ) 
  
    # Compute the interior + surface absorbed shortwave radiation
    Qts = transmitted_shortwave_radiation(i, j, kᴺ, grid, time, α, Qs)

    Qaℓ *= (1 - ℵᵢ)
    Qts *= (1 - ℵᵢ)
  
    Qss = shortwave_radiative_forcing(i, j, grid, penetrating_radiation, Qts, ocean_properties)

    # Compute the total heat flux
    ΣQao = (Qu + Qc + Qv) * (1 - ℵᵢ) + Qaℓ + Qss

    @inbounds begin
        # Write radiative components of the heat flux for diagnostic purposes
        atmos_ocean_fluxes.upwelling_longwave[i, j, 1] = Qu
        atmos_ocean_fluxes.downwelling_longwave[i, j, 1] = - Qaℓ
        atmos_ocean_fluxes.downwelling_shortwave[i, j, 1] = - Qts
    end

    # Convert from a mass flux to a volume flux (aka velocity)
    # by dividing with the ocean reference density.
    # Also switch the sign, for some reason we are given freshwater flux as positive down.
    ρₒ⁻¹ = 1 / ocean_properties.reference_density
    ΣFao = - Mp * ρₒ⁻¹

    # Add the contribution from the turbulent water vapor flux, which has
    # a different sign convention as the prescribed water mass fluxes (positive upwards)
    Fv = Mv * ρₒ⁻¹
    ΣFao += Fv

    # Compute fluxes for u, v, T, and S from momentum, heat, and freshwater fluxes
    τx = net_ocean_fluxes.u
    τy = net_ocean_fluxes.v
    Jᵀ = net_ocean_fluxes.T
    Jˢ = net_ocean_fluxes.S
    ℵ  = sea_ice_concentration
    cₒ = ocean_properties.heat_capacity

    @inbounds begin
        Qio  = get_possibly_zero_flux(sea_ice_ocean_fluxes, :interface_heat)[i, j, 1]
        Jˢio = get_possibly_zero_flux(sea_ice_ocean_fluxes, :salt)[i, j, 1]
        Jᵀao = ΣQao  * ρₒ⁻¹ / cₒ
        Jᵀio = Qio * ρₒ⁻¹ / cₒ
    
        # salinity flux > 0 extracts salinity from the ocean --- the opposite of a water vapor flux
        Jˢao = - Sₒ * ΣFao

        τxao = ℑxᶠᵃᵃ(i, j, 1, grid, τᶜᶜᶜ, ρₒ⁻¹, ℵ, ρτxao)
        τyao = ℑyᵃᶠᵃ(i, j, 1, grid, τᶜᶜᶜ, ρₒ⁻¹, ℵ, ρτyao)
        τxio = ρτxio[i, j, 1] * ρₒ⁻¹ * ℑxᶠᵃᵃ(i, j, 1, grid, ℵ)
        τyio = ρτyio[i, j, 1] * ρₒ⁻¹ * ℑyᵃᶠᵃ(i, j, 1, grid, ℵ)

        # Stresses
        τx[i, j, 1] = τxao + τxio
        τy[i, j, 1] = τyao + τyio

        # Tracer fluxes
        Jᵀ[i, j, 1] = Jᵀao + Jᵀio # Jᵀao is already multiplied by the sea ice concentration
        Jˢ[i, j, 1] = (1 - ℵᵢ) * Jˢao + Jˢio
    end
end
