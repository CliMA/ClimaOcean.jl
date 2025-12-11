using ClimaOcean.OceanSeaIceModels: FreezingLimitedOceanTemperature, FreezingLimitedCoupledModel
using ClimaSeaIce.SeaIceThermodynamics: melting_temperature

# No need to compute fluxes for this "sea ice model"
function compute_sea_ice_ocean_fluxes!(cm::FreezingLimitedCoupledModel)
    ocean = cm.ocean
    liquidus = cm.sea_ice.liquidus
    grid = ocean.model.grid
    arch = architecture(grid)
    Sₒ = ocean.model.tracers.S
    Tₒ = ocean.model.tracers.T

    launch!(arch, grid, :xyz, _above_freezing_ocean_temperature!, Tₒ, Sₒ, liquidus)

    return nothing
end

compute_atmosphere_sea_ice_fluxes!(cm::FreezingLimitedCoupledModel) = nothing

@kernel function _above_freezing_ocean_temperature!(Tₒ, Sₒ, liquidus)

    i, j, k = @index(Global, NTuple)

    @inbounds begin
        Sᵢ = Sₒ[i, j, k]
        Tᵢ = Tₒ[i, j, k]
    end

    Tₘ = melting_temperature(liquidus, Sᵢ)
    @inbounds Tₒ[i, j, k] = ifelse(Tᵢ < Tₘ, Tₘ, Tᵢ)
end
