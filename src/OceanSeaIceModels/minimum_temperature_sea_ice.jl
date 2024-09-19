using Oceananigans.Architectures: architecture

import ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: limit_fluxes_over_sea_ice!

"""
    struct MinimumTemperatureSeaIce{T}

The minimal possible sea ice representation, providing an "Insulating layer" on the surface.
Not really a ``model'' per se, however, it is the most simple way to make sure that temperature 
does not dip below freezing temperature.
All fluxes are shut down when the surface is below the `minimum_temperature` except for heating.

# Fields
- `minimum_temperature`: The minimum temperature of water.
"""
struct MinimumTemperatureSeaIce{T}
    minimum_temperature :: T
end

MinimumTemperatureSeaIce() = MinimumTemperatureSeaIce(-1.8)

function limit_fluxes_over_sea_ice!(grid, kernel_parameters, sea_ice::MinimumTemperatureSeaIce,
                                    staggered_velocity_fluxes,
                                    net_tracer_fluxes,
                                    ocean_temperature,
                                    ocean_salinity)

    launch!(architecture(grid), grid, kernel_parameters, _cap_fluxes_on_sea_ice!,
            staggered_velocity_fluxes,
            net_tracer_fluxes,
            grid,
            sea_ice.minimum_temperature,
            ocean_temperature)

    return nothing
end

@kernel function _cap_fluxes_on_sea_ice!(centered_velocity_fluxes,
                                         net_tracer_fluxes,
                                         grid,
                                         minimum_temperature,
                                         ocean_temperature)    

    i, j = @index(Global, NTuple)

    @inbounds begin
        Tₒ = ocean_temperature[i, j, 1]

        τx = centered_velocity_fluxes.u
        τy = centered_velocity_fluxes.v
        Jᵀ = net_tracer_fluxes.T
        Jˢ = net_tracer_fluxes.S
    
        sea_ice = Tₒ < minimum_temperature
        cooling_sea_ice = sea_ice & (Jᵀ[i, j, 1] > 0)

        # Don't allow the ocean to cool below the minimum temperature! (make sure it heats up though!)
        Jᵀ[i, j, 1] = ifelse(cooling_sea_ice, zero(grid), Jᵀ[i, j, 1]) 

        # If we are in a "sea ice" region we remove all fluxes
        Jˢ[i, j, 1] = ifelse(sea_ice, zero(grid), Jˢ[i, j, 1])
        τx[i, j, 1] = ifelse(sea_ice, zero(grid), τx[i, j, 1]) 
        τy[i, j, 1] = ifelse(sea_ice, zero(grid), τy[i, j, 1]) 
    end
end
