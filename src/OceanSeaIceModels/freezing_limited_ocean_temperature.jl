using ClimaSeaIce.SeaIceThermodynamics: LinearLiquidus

#####
##### A workaround when you don't have a sea ice model
#####

struct FreezingLimitedOceanTemperature{L}
    liquidus :: L
end

"""
    FreezingLimitedOceanTemperature(FT=Float64; liquidus=LinearLiquidus(FT))

The minimal possible sea ice representation, clipping the temperature below to the freezing point.
Not really a "model" per se, however, it is the most simple way to make sure that temperature
does not dip below freezing.

The melting temperature is a function of salinity and is controlled by the `liquidus`.
"""
FreezingLimitedOceanTemperature(FT::DataType=Oceananigans.defaults.FloatType; liquidus=LinearLiquidus(FT)) =
    FreezingLimitedOceanTemperature(liquidus)

const FreezingLimitedCoupledModel = OceanSeaIceModel{<:FreezingLimitedOceanTemperature}

# Extend interface methods to work with a `FreezingLimitedOceanTemperature`
sea_ice_concentration(::FreezingLimitedOceanTemperature) = ZeroField()
sea_ice_thickness(::FreezingLimitedOceanTemperature) = ZeroField()

# does not matter
reference_density(::FreezingLimitedOceanTemperature) = 0
heat_capacity(::FreezingLimitedOceanTemperature) = 0
time_step!(::FreezingLimitedOceanTemperature, Δt) = nothing

# FreezingLimitedOceanTemperature handles temperature limiting in compute_sea_ice_ocean_fluxes!
above_freezing_ocean_temperature!(ocean, ::FreezingLimitedOceanTemperature) = nothing

# No atmosphere-sea ice or sea ice-ocean interface for FreezingLimitedOceanTemperature
InterfaceComputations.default_ai_temperature(::FreezingLimitedOceanTemperature) = nothing
InterfaceComputations.atmosphere_sea_ice_interface(grid, atmos, ::FreezingLimitedOceanTemperature, args...) = nothing
InterfaceComputations.sea_ice_ocean_interface(grid, ::FreezingLimitedOceanTemperature, ocean; kwargs...) = nothing
InterfaceComputations.net_fluxes(::FreezingLimitedOceanTemperature) = nothing

