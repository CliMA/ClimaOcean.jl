using ClimaSeaIce.SeaIceThermodynamics: LinearLiquidus

#####
##### A simple-minded yet effective "sea ice model"
#####

struct FreezingLimitedOceanTemperature{L}
    liquidus :: L
end

"""
    FreezingLimitedOceanTemperature(FT::DataType) = FreezingLimitedOceanTemperature(LinearLiquidus(FT))

The minimal possible sea ice representation, providing an "Insulating layer" on the surface and clipping the 
temperature below to the freezing point. Not really a ``model'' per se, however, 
it is the most simple way to make sure that temperature does not dip below freezing. 
All fluxes are shut down when the surface is below the `T < Tₘ` except for heating to allow temperature to increase.

the melting temperature is a function of salinity and is controlled by the `liquidus`.
"""
FreezingLimitedOceanTemperature(FT::DataType=Float64) = FreezingLimitedOceanTemperature(LinearLiquidus(FT))

const FreezingLimitedCoupledModel = OceanSeaIceModel{<:FreezingLimitedOceanTemperature}

sea_ice_concentration(::FreezingLimitedOceanTemperature) = nothing

function compute_sea_ice_ocean_fluxes!(cm::FreezingLimitedCoupledModel)
    ocean = cm.ocean
    liquidus = cm.sea_ice.liquidus
    grid = ocean.model.grid
    arch = architecture(grid)
    Sₒ = ocean.model.tracers.S
    Tₒ = ocean.model.tracers.T

    launch!(arch, grid, :xyz,  above_freezing_ocean_temperature!, Tₒ, Sₒ, liquidus)

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

function limit_fluxes_over_sea_ice!(grid, kernel_parameters,
                                    sea_ice::FreezingLimitedOceanTemperature,
                                    centered_velocity_fluxes,
                                    net_tracer_fluxes,
                                    ocean_temperature,
                                    ocean_salinity)

    launch!(architecture(grid), grid, kernel_parameters, _limit_fluxes_over_sea_ice!,
            centered_velocity_fluxes,
            net_tracer_fluxes,
            grid,
            sea_ice.liquidus,
            ocean_temperature,
            ocean_salinity)

    return nothing
end

@kernel function _limit_fluxes_over_sea_ice!(centered_velocity_fluxes,
                                             net_tracer_fluxes,
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

        τx = centered_velocity_fluxes.u
        τy = centered_velocity_fluxes.v
        Jᵀ = net_tracer_fluxes.T
        Jˢ = net_tracer_fluxes.S

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
