using Printf
using Oceananigans.Operators: ℑxᶠᵃᵃ, ℑyᵃᶠᵃ

using ClimaOcean.OceanSeaIceModels: sea_ice_concentration

@inline computed_sea_ice_ocean_fluxes(interface) = interface.fluxes
@inline computed_sea_ice_ocean_fluxes(::Nothing) = nothing

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
    atmosphere_fields = coupled_model.interfaces.exchanger.near_surface_atmosphere_state

    downwelling_radiation = (Qs = atmosphere_fields.Qs.data,
                             Qℓ = atmosphere_fields.Qℓ.data)

    freshwater_flux = atmosphere_fields.Mp.data

    ice_concentration = sea_ice_concentration(sea_ice)
    ocean_salinity = ocean.model.tracers.S
    atmos_ocean_properties = coupled_model.interfaces.atmosphere_ocean_interface.properties
    ocean_properties = coupled_model.interfaces.ocean_properties
    kernel_parameters = interface_kernel_parameters(grid)

    ocean_surface_temperature = coupled_model.interfaces.atmosphere_ocean_interface.temperature

    launch!(arch, grid, kernel_parameters, 
            _assemble_net_ocean_fluxes!,
            net_ocean_fluxes,
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

@inline τᶜᶜᶜ(i, j, k, grid, ρₒ⁻¹, ℵ, ρτᶜᶜᶜ) = @inbounds ρₒ⁻¹ * (1 - ℵ[i, j, 1]) 

@kernel function _assemble_net_ocean_fluxes!(net_ocean_fluxes,
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

    assemble_net_atmosphere_ocean_fluxes!(i, j, grid, 
                                          net_ocean_fluxes,
                                          clock,
                                          atmos_ocean_fluxes,
                                          ocean_salinity,
                                          ocean_surface_temperature,
                                          downwelling_radiation,
                                          freshwater_flux,
                                          atmos_ocean_properties,
                                          ocean_properties)

    add_sea_ice_ocean_fluxes!(i, j, grid, 
                              net_ocean_fluxes,
                              sea_ice_ocean_fluxes,
                              sea_ice_concentration,
                              ocean_salinity,
                              ocean_surface_temperature)
end

@inline add_sea_ice_ocean_fluxes!(i, j, grid, net_ocean_fluxes, ::Nothing, args...) = nothing

@inline function add_sea_ice_ocean_fluxes!(i, j, grid, 
                                           net_ocean_fluxes,
                                           sea_ice_ocean_fluxes,
                                           sea_ice_concentration, args...)

    # Compute fluxes for u, v, T, S from momentum, heat, and freshwater fluxes
    τx = net_ocean_fluxes.u
    τy = net_ocean_fluxes.v
    Jᵀ = net_ocean_fluxes.T
    Jˢ = net_ocean_fluxes.S

    ρₒ⁻¹ = 1 / ocean_properties.reference_density
    cₒ   = ocean_properties.heat_capacity

    @inbounds begin
        ℵ    = sea_ice_concentration[i, j, 1]
        Qio  = sea_ice_ocean_fluxes.heat[i, j, 1]
        Jˢio = sea_ice_ocean_fluxes.salt[i, j, 1]
    end
    
    Jᵀio = Qio * ρₒ⁻¹ / cₒ

    @inbounds begin
        τx[i, j, 1] *= 1 - ℑxᶠᵃᵃ(i, j, 1, grid, ℵ)
        τy[i, j, 1] *= 1 - ℑyᵃᶠᵃ(i, j, 1, grid, ℵ)
        Jᵀ[i, j, 1] = (1 - ℵ) * Jᵀ[i, j, 1] + Jᵀio
        Jˢ[i, j, 1] = (1 - ℵ) * Jˢ[i, j, 1] + Jˢio
    end
end

@inline function assemble_net_atmosphere_ocean_fluxes!(i, j, grid, 
                                                       net_ocean_fluxes,
                                                       clock,
                                                       atmos_ocean_fluxes,
                                                       ocean_salinity,
                                                       ocean_surface_temperature,
                                                       downwelling_radiation,
                                                       freshwater_flux,
                                                       atmos_ocean_properties,
                                                       ocean_properties)
    kᴺ = size(grid, 3)
    time = Time(clock.time)
                                          
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
        ρτx = atmos_ocean_fluxes.x_momentum[i, j, 1]    # zonal momentum flux
        ρτy = atmos_ocean_fluxes.y_momentum[i, j, 1]    # meridional momentum flux
    end

    # Compute radiation fluxes
    σ = atmos_ocean_properties.radiation.σ
    α = stateindex(atmos_ocean_properties.radiation.α, i, j, kᴺ, grid, time)
    ϵ = stateindex(atmos_ocean_properties.radiation.ϵ, i, j, kᴺ, grid, time)
    Qu = upwelling_radiation(Tₛ, σ, ϵ)
    Qr = (; Qs, Qℓ)
    Qd = net_downwelling_radiation(Qr, α, ϵ)
    ΣQao = Qd + Qu + Qc + Qv

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

    ρₒ⁻¹ = 1 / ocean_properties.reference_density
    cₒ   = ocean_properties.heat_capacity

    τxao = ℑxᶠᵃᵃ(i, j, 1, grid, ρτx) * ρₒ⁻¹
    τyao = ℑyᵃᶠᵃ(i, j, 1, grid, ρτy) * ρₒ⁻¹
    Jᵀao = ΣQao  * ρₒ⁻¹ / cₒ
    Jˢao = - Sₒ * ΣFao

    @inbounds begin
        τx[i, j, 1] = τxao
        τy[i, j, 1] = τyao
        Jᵀ[i, j, 1] = Jᵀao
        Jˢ[i, j, 1] = Jˢao
    end
end

compute_net_sea_ice_fluxes!(coupled_model) = nothing

#=
function compute_net_sea_ice_fluxes!(coupled_model)
    ocean = coupled_model.ocean
    sea_ice = coupled_model.sea_ice
    grid = ocean.model.grid
    arch = architecture(grid)
    clock = coupled_model.clock

    top_heat_flux = coupled_model.interfaces.net_fluxes.sea_ice_top
    atmos_ocean_fluxes = coupled_model.interfaces.atmosphere_ocean_interface.fluxes
    sea_ice_ocean_fluxes = coupled_model.interfaces.sea_ice_ocean_interface.fluxes

    # Simplify NamedTuple to reduce parameter space consumption.
    # See https://github.com/CliMA/ClimaOcean.jl/issues/116.
    atmosphere_fields = coupled_model.interfaces.exchanger.near_surface_atmosphere_state

    downwelling_radiation = (Qs = atmosphere_fields.Qs.data,
                             Qℓ = atmosphere_fields.Qℓ.data)

    freshwater_flux = atmosphere_fields.Mp.data

    sea_ice_concentration = sea_ice.model.ice_concentration
    ocean_salinity = ocean.model.tracers.S
    atmos_ocean_properties = coupled_model.interfaces.atmosphere_ocean_interface.properties
    ocean_properties = coupled_model.interfaces.ocean_properties
    kernel_parameters = interface_kernel_parameters(grid)

    ocean_surface_temperature = coupled_model.interfaces.atmosphere_ocean_interface.temperature

    launch!(arch, grid, kernel_parameters, 
            _assemble_net_sea_ice_fluxes!,
            net_ocean_fluxes,
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

    return nothing
end

@kernel function _assemble_atmosphere_sea_ice_fluxes!(net_heat_flux,
                                                      grid,
                                                      clock,
                                                      surface_temperature,
                                                      surface_temperature_units,
                                                      turbulent_fluxes,
                                                      downwelling_radiation,
                                                      stefan_boltzmann_constant,
                                                      albedo,
                                                      emissivity)

    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3)
    time = Time(clock.time)

    @inbounds begin
        Ts = surface_temperature[i, j, kᴺ]
        Ts = convert_to_kelvin(surface_temperature_units, Ts)

        Qs  = downwelling_radiation.shortwave[i, j, 1]
        Qℓ  = downwelling_radiation.longwave[i, j, 1]
        Qc  = turbulent_fluxes.sensible_heat[i, j, 1] # sensible or "conductive" heat flux
        Qv  = turbulent_fluxes.latent_heat[i, j, 1]   # latent heat flux
    end

    # Compute radiation fluxes
    σ = stefan_boltzmann_constant
    α = stateindex(albedo, i, j, 1, grid, time)
    ϵ = stateindex(emissivity, i, j, 1, grid, time)
    Qu = upwelling_radiation(Ts, σ, ϵ)
    Qd = net_downwelling_radiation(i, j, grid, time, α, ϵ, Qs, Qℓ)

    ΣQ = Qd + Qu + Qc + Qv

    # Mask fluxes over land for convenience
    inactive = inactive_node(i, j, kᴺ, grid, c, c, c)

    @inbounds net_heat_flux[i, j, 1] = ifelse(inactive, zero(grid), ΣQ)
end
=#
