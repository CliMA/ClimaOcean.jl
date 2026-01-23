using ClimaSeaIce.SeaIceThermodynamics: melting_temperature
using ClimaSeaIce.SeaIceThermodynamics: LinearLiquidus
using ClimaOcean.OceanSeaIceModels
using ClimaOcean.OceanSeaIceModels: NoSeaIceInterface
using ClimaOcean.OceanSeaIceModels.InterfaceComputations

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

const FreezingLimitedCoupledModel = OceanSeaIceModel{<:FreezingLimitedOceanTemperature, A, O, <:NoSeaIceInterface} where {A, O}

# Extend interface methods to work with a `FreezingLimitedOceanTemperature`
sea_ice_concentration(::FreezingLimitedOceanTemperature) = ZeroField()
sea_ice_thickness(::FreezingLimitedOceanTemperature) = ZeroField()

# does not matter
reference_density(::FreezingLimitedOceanTemperature) = 0
heat_capacity(::FreezingLimitedOceanTemperature) = 0
time_step!(::FreezingLimitedOceanTemperature, Δt) = nothing

# FreezingLimitedOceanTemperature handles temperature limiting in compute_sea_ice_ocean_fluxes!
OceanSeaIceModels.above_freezing_ocean_temperature!(ocean, grid, ::FreezingLimitedOceanTemperature) = nothing

# No atmosphere-sea ice or sea ice-ocean interface for FreezingLimitedOceanTemperature
InterfaceComputations.default_ai_temperature(::FreezingLimitedOceanTemperature) = nothing
InterfaceComputations.ThreeEquationHeatFlux(::FreezingLimitedOceanTemperature) = nothing
InterfaceComputations.atmosphere_sea_ice_interface(grid, atmos, ::FreezingLimitedOceanTemperature, args...) = nothing
InterfaceComputations.atmosphere_sea_ice_interface(grid, ::Nothing, ::FreezingLimitedOceanTemperature, args...) = nothing
InterfaceComputations.sea_ice_ocean_interface(grid, ::FreezingLimitedOceanTemperature, ocean, flux_formulation; kwargs...) = nothing
InterfaceComputations.sea_ice_ocean_interface(grid, ::FreezingLimitedOceanTemperature, ::Nothing, flux_formulation; kwargs...) = nothing
InterfaceComputations.sea_ice_ocean_interface(grid, ::FreezingLimitedOceanTemperature, ocean, ::ThreeEquationHeatFlux; kwargs...) = nothing
InterfaceComputations.sea_ice_ocean_interface(grid, ::FreezingLimitedOceanTemperature, ::Nothing, ::ThreeEquationHeatFlux; kwargs...) = nothing

InterfaceComputations.net_fluxes(::FreezingLimitedOceanTemperature) = nothing

const OnlyOceanwithFreezingLimited      = OceanSeaIceModel{<:FreezingLimitedOceanTemperature, <:Nothing, <:Any}
const OnlyAtmospherewithFreezingLimited = OceanSeaIceModel{<:FreezingLimitedOceanTemperature, <:Any,     <:Nothing}
const SingleComponentPlusFreezingLimited = Union{OnlyAtmospherewithFreezingLimited, OnlyOceanwithFreezingLimited}

# Also for the ocean nothing really happens here
OceanSeaIceModels.update_net_fluxes!(::SingleComponentPlusFreezingLimited, ocean::Simulation{<:HydrostaticFreeSurfaceModel}) = nothing

# No need to compute fluxes for this "sea ice model"
InterfaceComputations.compute_atmosphere_sea_ice_fluxes!(cm::FreezingLimitedCoupledModel) = nothing

# Same for the sea_ice ocean fluxes
function InterfaceComputations.compute_sea_ice_ocean_fluxes!(cm::FreezingLimitedCoupledModel)
    ocean = cm.ocean
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

#####
##### Chekpointing (not needed for FreezingLimitedOceanTemperature)
#####

import Oceananigans: prognostic_state, restore_prognostic_state!

prognostic_state(::FreezingLimitedOceanTemperature) = nothing
restore_prognostic_state!(flt::FreezingLimitedOceanTemperature, state) = flt
restore_prognostic_state!(flt::FreezingLimitedOceanTemperature, ::Nothing) = flt
