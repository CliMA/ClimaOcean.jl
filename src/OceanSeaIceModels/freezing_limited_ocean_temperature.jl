using ClimaSeaIce.SeaIceThermodynamics: LinearLiquidus

import ClimaOcean.OceanSeaIceModels.InterfaceComputations: add_sea_ice_ocean_fluxes!, 
                                                           computed_sea_ice_ocean_fluxes,
                                                           atmosphere_sea_ice_interface

#####
##### A workaround when you don't have a sea ice model
#####

struct FreezingLimitedOceanTemperature{L}
    liquidus :: L
end

"""
    FreezingLimitedOceanTemperature(FT=Float64)

The minimal possible sea ice representation, providing an "Insulating layer" on
the surface and clipping the temperature below to the freezing point. Not really
a "model"' per se, however, it is the most simple way to make sure that temperature
does not dip below freezing. All fluxes are shut down when the surface is below
the `T < Tₘ` except for heating to allow temperature to increase.

The melting temperature is a function of salinity and is controlled by the `liquidus`.
"""
FreezingLimitedOceanTemperature(FT::DataType=Float64) = FreezingLimitedOceanTemperature(LinearLiquidus(FT))

const FreezingLimitedCoupledModel = OceanSeaIceModel{<:FreezingLimitedOceanTemperature}

sea_ice_concentration(::FreezingLimitedOceanTemperature) = nothing
sea_ice_thickness(::FreezingLimitedOceanTemperature) = nothing

# does not matter
reference_density(::FreezingLimitedOceanTemperature) = 0
heat_capacity(::FreezingLimitedOceanTemperature) = 0

function compute_sea_ice_ocean_fluxes!(cm::FreezingLimitedCoupledModel)
    ocean = cm.ocean
    liquidus = cm.sea_ice.liquidus
    grid = ocean.model.grid
    arch = architecture(grid)
    Sₒ = ocean.model.tracers.S
    Tₒ = ocean.model.tracers.T

    launch!(arch, grid, :xyz, above_freezing_ocean_temperature!, Tₒ, Sₒ, liquidus)

    return nothing
end

@kernel function above_freezing_ocean_temperature!(Tₒ, Sₒ, liquidus)

    i, j, k = @index(Global, NTuple)

    @inbounds begin
        Sᵢ = Sₒ[i, j, k]
        Tᵢ = Tₒ[i, j, k]
    end

    Tₘ = melting_temperature(liquidus, Sᵢ)
    @inbounds Tₒ[i, j, k] = ifelse(Tᵢ < Tₘ, Tₘ, Tᵢ)
end

@kernel function _adjust_fluxes_over_sea_ice!(net_fluxes,
                                              grid,
                                              liquidus,
                                              ocean_temperature,
                                              ocean_salinity)

    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3)
    
    @inbounds begin
        Tₒ = ocean_temperature[i, j, kᴺ]
        Sₒ = ocean_salinity[i, j, kᴺ]

        Tₘ = melting_temperature(liquidus, Sₒ)

        τx = net_fluxes.u
        τy = net_fluxes.v
        Jᵀ = net_fluxes.T
        Jˢ = net_fluxes.S

        sea_ice = Tₒ < Tₘ
        cooling_sea_ice = sea_ice & (Jᵀ[i, j, 1] > 0)

        @show i, j, cooling_sea_ice, sea_ice
        # Don't allow the ocean to cool below the minimum temperature! (make sure it heats up though!)
        Jᵀ[i, j, 1] = ifelse(cooling_sea_ice, zero(grid), Jᵀ[i, j, 1]) 

        # If we are in a "sea ice" region we remove all fluxes
        Jˢ[i, j, 1] = ifelse(sea_ice, zero(grid), Jˢ[i, j, 1])
        τx[i, j, 1] = ifelse(sea_ice, zero(grid), τx[i, j, 1])
        τy[i, j, 1] = ifelse(sea_ice, zero(grid), τy[i, j, 1])
    end
end

# Extend interface methods to work with a `FreezingLimitedOceanTemperature`

atmosphere_sea_ice_interface(sea_ice::FreezingLimitedOceanTemperature, args...) = sea_ice

@inline computed_sea_ice_ocean_fluxes(interface::FreezingLimitedOceanTemperature) = interface

@inline function add_sea_ice_ocean_fluxes!(i, j, grid,
                                           net_ocean_fluxes,
                                           sea_ice_ocean_fluxes::FreezingLimitedOceanTemperature,
                                           sea_ice_concentration,
                                           ocean_salinity,
                                           ocean_surface_temperature)
    
    kᴺ = size(grid, 3)
                                           
    @inbounds begin
        Tₒ = ocean_surface_temperature[i, j, kᴺ]
        Sₒ = ocean_salinity[i, j, kᴺ]

        Tₘ = melting_temperature(sea_ice_ocean_fluxes.liquidus, Sₒ)

        τx = net_ocean_fluxes.u
        τy = net_ocean_fluxes.v
        Jᵀ = net_ocean_fluxes.T
        Jˢ = net_ocean_fluxes.S

        sea_ice = Tₒ < Tₘ
        cooling_sea_ice = sea_ice & (Jᵀ[i, j, 1] > 0)

        # Don't allow the ocean to cool below the minimum temperature! (make sure it heats up though!)
        Jᵀ[i, j, 1] = ifelse(cooling_sea_ice, zero(grid), Jᵀ[i, j, 1]) 

        # If we are in a "sea ice" region we remove all fluxes
        Jˢ[i, j, 1] = ifelse(sea_ice, zero(grid), Jˢ[i, j, 1])
        τx[i, j, 1] = ifelse(sea_ice, zero(grid), τx[i, j, 1])
        τy[i, j, 1] = ifelse(sea_ice, zero(grid), τy[i, j, 1])
    end
end
