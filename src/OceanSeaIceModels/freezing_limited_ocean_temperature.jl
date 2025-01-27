using ClimaSeaIce.SeaIceThermodynamics: LinearLiquidus

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
    net_fluxes = cm.interfaces.net_fluxes.ocean_surface
    sea_ice = cm.sea_ice
    
    launch!(arch, grid, :xy, _adjust_fluxes_over_sea_ice!,
            net_fluxes,
            grid,
            sea_ice.liquidus,
            Tₒ, Sₒ)

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

        # Don't allow the ocean to cool below the minimum temperature! (make sure it heats up though!)
        Jᵀ[i, j, 1] = ifelse(cooling_sea_ice, zero(grid), Jᵀ[i, j, 1]) 

        # If we are in a "sea ice" region we remove all fluxes
        Jˢ[i, j, 1] = ifelse(sea_ice, zero(grid), Jˢ[i, j, 1])
        τx[i, j, 1] = ifelse(sea_ice, zero(grid), τx[i, j, 1])
        τy[i, j, 1] = ifelse(sea_ice, zero(grid), τy[i, j, 1])
    end
end
