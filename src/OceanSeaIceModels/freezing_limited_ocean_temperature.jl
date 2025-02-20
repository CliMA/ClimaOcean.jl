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
   
    return nothing
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