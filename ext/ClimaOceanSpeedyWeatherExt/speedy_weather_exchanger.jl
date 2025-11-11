using Oceananigans
using Oceananigans.BoundaryConditions
using Oceananigans.Grids: architecture
using Oceananigans.Utils: launch!
using Oceananigans.Operators: intrinsic_vector
using XESMF

using ClimaOcean.OceanSeaIceModels: sea_ice_concentration

# TODO: Implement conservative regridding when ready
# using ConservativeRegridding 
# using GeoInterface: Polygon, LinearRing
import ClimaOcean.OceanSeaIceModels:
    compute_net_atmosphere_fluxes!

import ClimaOcean.OceanSeaIceModels.InterfaceComputations:
    atmosphere_exchanger,
    initialize!,
    StateExchanger,
    interpolate_atmosphere_state!

# For the moment the workflow is:
# 1. Perform the regridding on the CPU
# 2. Eventually copy the regridded fields to the GPU
# If this work we can
# 1. Copy speedyweather gridarrays to the GPU
# 2. Perform the regridding on the GPU
function atmosphere_exchanger(atmosphere::SpeedySimulation, exchange_grid, exchange_atmosphere_state)

    # Figure this out:
    spectral_grid = atmosphere.model.spectral_grid
    FT = eltype(exchange_atmosphere_state.u)

    # TODO: Implement a conservative regridder when ready
    ocean_atmosphere_regridder = XESMF.Regridder(spectral_grid, exchange_grid)
    atmosphere_ocean_regridder = XESMF.Regridder(exchange_grid, spectral_grid)
    exchanger = (; ocean_atmosphere_regridder, atmosphere_ocean_regridder)

    return exchanger
end

# Regrid the atmospheric state on the exchange grid
function interpolate_atmosphere_state!(interfaces, atmos::SpeedySimulation, coupled_model)
    atmosphere_exchanger = interfaces.exchanger.atmosphere_exchanger
    regrid! = atmosphere_exchanger.atmosphere_ocean_regridder
    exchange_grid  = interfaces.exchanger.exchange_grid
    exchange_state = interfaces.exchanger.exchange_atmosphere_state
    surface_layer = atmos.model.spectral_grid.nlayers

    ua  = RingGrids.field_view(atmos.diagnostic_variables.grid.u_grid,     :, surface_layer).data
    va  = RingGrids.field_view(atmos.diagnostic_variables.grid.v_grid,     :, surface_layer).data
    Ta  = RingGrids.field_view(atmos.diagnostic_variables.grid.temp_grid,  :, surface_layer).data
    qa  = RingGrids.field_view(atmos.diagnostic_variables.grid.humid_grid, :, surface_layer).data
    pa  = exp.(atmos.diagnostic_variables.grid.pres_grid.data)
    Qsa = atmos.diagnostic_variables.physics.surface_shortwave_down.data
    Qla = atmos.diagnostic_variables.physics.surface_longwave_down.data
    Mpa = atmos.diagnostic_variables.physics.total_precipitation_rate.data

    regrid!(vec(interior(exchange_state.u)),  ua)
    regrid!(vec(interior(exchange_state.v)),  va)
    regrid!(vec(interior(exchange_state.T)),  Ta)
    regrid!(vec(interior(exchange_state.q)),  qa)
    regrid!(vec(interior(exchange_state.p)),  pa)
    regrid!(vec(interior(exchange_state.Qs)), Qsa)
    regrid!(vec(interior(exchange_state.Qℓ)), Qla)
    regrid!(vec(interior(exchange_state.Mp)), Mpa)

    arch = architecture(exchange_grid)

    u = exchange_state.u
    v = exchange_state.v

    launch!(arch, exchange_grid, :xy, _rotate_winds!, u, exchange_grid, v)

    fill_halo_regions!((u, v))
    fill_halo_regions!(exchange_state.T)
    fill_halo_regions!(exchange_state.q)
    fill_halo_regions!(exchange_state.p)
    fill_halo_regions!(exchange_state.Qs)
    fill_halo_regions!(exchange_state.Qℓ)
    fill_halo_regions!(exchange_state.Mp)

    return nothing
end

@kernel function _rotate_winds!(u, grid, v)
    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3) 
    uₑ, vₑ = intrinsic_vector(i, j, kᴺ, grid, u, v)
    @inbounds u[i, j, kᴺ] = uₑ
    @inbounds v[i, j, kᴺ] = vₑ
end

# TODO: Fix the coupling with the sea ice model and make sure that 
# the this function works also for sea_ice=nothing and on GPUs without
# needing to allocate memory.
function compute_net_atmosphere_fluxes!(coupled_model::SpeedyCoupledModel)
    atmos = coupled_model.atmosphere
    grid  = coupled_model.interfaces.exchanger.exchange_grid
    regrid! = coupled_model.interfaces.exchanger.atmosphere_exchanger.ocean_atmosphere_regridder
    ao_fluxes = coupled_model.interfaces.atmosphere_ocean_interface.fluxes
    ai_fluxes = coupled_model.interfaces.atmosphere_sea_ice_interface.fluxes

    Qco = ao_fluxes.sensible_heat
    Qci = ai_fluxes.sensible_heat
    Mvo = ao_fluxes.water_vapor
    Mvi = ai_fluxes.water_vapor
    ℵ   = interior(sea_ice_concentration(coupled_model.sea_ice))

    # All the location of these fluxes will change
    Qca = atmos.prognostic_variables.ocean.sensible_heat_flux.data
    Mva = atmos.prognostic_variables.ocean.surface_humidity_flux.data
    sst = atmos.prognostic_variables.ocean.sea_surface_temperature.data
    To  = coupled_model.interfaces.atmosphere_ocean_interface.temperature
    Ti  = coupled_model.interfaces.atmosphere_sea_ice_interface.temperature

    # TODO: Figure out how we are going to deal with upwelling radiation
    # TODO: regrid longwave rather than a mixed surface temperature 
    regrid!(Qca, interior(Qco) .* (1 - ℵ) .+ ℵ .* interior(Qci))
    regrid!(Mva, interior(Mvo) .* (1 - ℵ) .+ ℵ .* interior(Mvi))
    regrid!(sst, interior(To) .* (1 - ℵ) .+ ℵ .* interior(Ti) .+ 273.15)

    return nothing
end

# Simple case -> there is no sea ice!
function compute_net_atmosphere_fluxes!(coupled_model::SpeedyNoSeaIceCoupledModel)
    atmos = coupled_model.atmosphere
    grid  = coupled_model.interfaces.exchanger.exchange_grid
    regridder = coupled_model.interfaces.exchanger.atmosphere_exchanger.ocean_atmosphere_regridder
    ao_fluxes = coupled_model.interfaces.atmosphere_ocean_interface.fluxes
    Qco = ao_fluxes.sensible_heat
    Mvo = ao_fluxes.water_vapor

    # All the location of these fluxes will change
    Qca = atmos.prognostic_variables.ocean.sensible_heat_flux.data
    Mva = atmos.prognostic_variables.ocean.surface_humidity_flux.data
    sst = atmos.prognostic_variables.ocean.sea_surface_temperature.data
    To  = coupled_model.interfaces.atmosphere_ocean_interface.temperature

    # TODO: Figure out how we are going to deal with upwelling radiation
    # TODO: regrid longwave rather than a mixed surface temperature 
    regrid!(Qca, regridder, interior(Qco))
    regrid!(Mva, regridder, interior(Mvo))
    regrid!(sst, regridder, interior(To) .+ 273.15)

    return nothing
end
