using ClimaSeaIce.SeaIceThermodynamics: LinearLiquidus

#####
##### A workaround when you don't have a sea ice model
#####

struct FreezingLimitedOceanTemperature{L, C}
    liquidus :: L
    ice_concentration :: C
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
function FreezingLimitedOceanTemperature(grid; FT::DataType=Float64) 
    ice_concentration = Field{Center, Center, Nothing}(grid)
    return FreezingLimitedOceanTemperature(LinearLiquidus(FT), ice_concentration)
end

const FreezingLimitedCoupledModel = OceanSeaIceModel{<:FreezingLimitedOceanTemperature}

# Extend interface methods to work with a `FreezingLimitedOceanTemperature`
sea_ice_concentration(::FreezingLimitedOceanTemperature) = ZeroField()
sea_ice_thickness(::FreezingLimitedOceanTemperature) = nothing

# does not matter
reference_density(::FreezingLimitedOceanTemperature) = 0
heat_capacity(::FreezingLimitedOceanTemperature) = 0

function compute_sea_ice_ocean_fluxes!(cm::FreezingLimitedCoupledModel)
    ocean = cm.ocean
    ℵ = cm.sea_ice.ice_concentration
    liquidus = cm.sea_ice.liquidus
    grid = ocean.model.grid
    arch = architecture(grid)
    Sₒ = ocean.model.tracers.S
    Tₒ = ocean.model.tracers.T

    launch!(arch, grid, :xyz, _above_freezing_ocean_temperature!, Tₒ, Sₒ, liquidus)
    launch!(arch, grid, :xy, _compute_ice_concentration!, ℵ, grid, ocean.model.tracers.T, ocean.model.tracers.S, liquidus)
   
    return nothing
end

@kernel function _compute_ice_concentration!(ℵ, grid, 
                                             ocean_surface_temperature, 
                                             ocean_salinity, 
                                             liquidus)

    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3)
                                           
    @inbounds begin
        Tₒ = ocean_surface_temperature[i, j, kᴺ]
        Sₒ = ocean_salinity[i, j, kᴺ]

        Tₘ = melting_temperature(liquidus, Sₒ)

        sea_ice = Tₒ < Tₘ
        ℵ[i, j, 1] = ifelse(sea_ice, one(grid), zero(grid))
    end
end

@kernel function _above_freezing_ocean_temperature!(Tₒ, Sₒ, liquidus)

    i, j, k = @index(Global, NTuple)

    @inbounds begin
        Sᵢ = Sₒ[i, j, k]
        Tᵢ = Tₒ[i, j, k]
    end

    Tₘ = melting_temperature(liquidus, Sᵢ)
    @inbounds Tₒ[i, j, k] = ifelse(Tᵢ < Tₘ, Tₘ, Tᵢ)
end