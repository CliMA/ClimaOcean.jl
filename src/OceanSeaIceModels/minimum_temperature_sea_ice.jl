import ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: limit_fluxes_over_sea_ice!

"""
    struct MinimumTemperatureSeaIce{T}

A struct representing the minimum temperature of sea ice.

# Fields
- `minimum_temperature`: The minimum temperature of sea ice.
"""
struct MinimumTemperatureSeaIce{T}
    minimum_temperature :: T
end

MinimumTemperatureSeaIce() = MinimumTemperatureSeaIce(-1.8)

function limit_fluxes_over_sea_ice!(grid, kernel_parameters, sea_ice::MinimumTemperatureSeaIce,
                                    similarity_theory_fields,
                                    ocean_temperature,
                                    ocean_salinity)

    launch!(grid, kernel_parameters, _cap_fluxes_on_sea_ice!,
            similarity_theory_fields,
            grid, 
            sea_ice.minimum_temperature,
            ocean_temperature)

    return nothing
end

@kernel function _cap_fluxes_on_sea_ice!(similarity_theory_fields,
                                         grid, 
                                         minimum_temperature,
                                         ocean_temperature)    

    i, j = @index(Global, NTuple)
    fields = similarity_theory_fields

    @inbounds begin
        Tₒ = ocean_temperature[i, j, 1]

        Qc = fields.sensible_heat[i, j, 1] # sensible or "conductive" heat flux
        Qv = fields.latent_heat[i, j, 1]   # latent heat flux
        Mv = fields.water_vapor[i, j, 1]   # mass flux of water vapor
        τx = fields.x_momentum[i, j, 1]    # zonal momentum flux
        τy = fields.y_momentum[i, j, 1]    # meridional momentum flux

        sea_ice = Tₒ < minimum_temperature
        cooling_sea_ice = sea_ice & (Qc > 0)
        evaporating_sea_ice = sea_ice & (Qv > 0)

        # Don't allow the ocean to cool below the minimum temperature!
        fields.sensible_heat[i, j, 1] = ifelse(cooling_sea_ice, zero(grid), Qc) # sensible or "conductive" heat flux
        fields.latent_heat[i, j, 1]   = ifelse(evaporating_sea_ice, zero(grid), Qv) # latent heat flux

        # If we are in a "sea ice" region we remove all fluxes
        fields.water_vapor[i, j, 1]   = ifelse(sea_ice, zero(grid), Mv) # mass flux of water vapor
        fields.x_momentum[i, j, 1]    = ifelse(sea_ice, zero(grid), τx) # zonal momentum flux
        fields.y_momentum[i, j, 1]    = ifelse(sea_ice, zero(grid), τy) # meridional momentum flux
    end
end
