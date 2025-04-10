# This document lays out the functions that must be extended to
# use an atmospheric simulation in ClimaOcean.

import ClimaOcean.OceanSeaIceModels:
    compute_net_atmosphere_fluxes!

import ClimaOcean.OceanSeaIceModels.InterfaceComputations:
    atmosphere_exchanger,
    initialize!,
    StateExchanger,
    interpolate_atmosphere_state!

using ClimaOcean.OceanSeaIceModels: OceanSeaIceModel
using Oceananigans
using Thermodynamics

"""
    interpolate_atmospheric_state!(surface_atmosphere_state, 
                                        interpolated_prescribed_freshwater_flux, 
                                        atmos::ClimaAtmosSimulation, 
                                        grid, clock)

Interpolate the atmospheric state in `atmos` to `surface_atmospheric_state`, a
the collection of `Field`s needed to compute turbulent fluxes.
"""
function interpolate_atmosphere_state!(interfaces, atmosphere::SpeedySimulation, coupled_model)

    # Plan:
    # 1. transfer atmos data to GPU (req ability to represent atmos on GPU)
    # 2. interpolate from atmos grid to ocean grid

    # RingGrids.interpolate!(vec(view(ua, :, :, 1)), atmos.diagnostic_variables.grid.u_grid[:, end],            interpolator)
    # RingGrids.interpolate!(vec(view(va, :, :, 1)), atmos.diagnostic_variables.grid.v_grid[:, end],            interpolator)
    # RingGrids.interpolate!(vec(view(Ta, :, :, 1)), atmos.diagnostic_variables.grid.temp_grid[:, end],         interpolator)
    # RingGrids.interpolate!(vec(view(qa, :, :, 1)), atmos.diagnostic_variables.grid.humid_grid[:, end],        interpolator)
    # RingGrids.interpolate!(vec(view(pa, :, :, 1)), exp.(atmos.diagnostic_variables.grid.pres_grid[:, end]),   interpolator)
    # RingGrids.interpolate!(vec(view(Qs, :, :, 1)), atmos.diagnostic_variables.physics.surface_shortwave_down, interpolator)
    # RingGrids.interpolate!(vec(view(Qℓ, :, :, 1)), atmos.diagnostic_variables.physics.surface_longwave_down,  interpolator)

    # Get the atmospheric state on the ocean grid
    # ua = on_architecture(Oceananigans.CPU(), surface_atmosphere_state.u)
    # va = on_architecture(Oceananigans.CPU(), surface_atmosphere_state.v)
    # Ta = on_architecture(Oceananigans.CPU(), surface_atmosphere_state.T)
    # qa = on_architecture(Oceananigans.CPU(), surface_atmosphere_state.q)
    # pa = on_architecture(Oceananigans.CPU(), surface_atmosphere_state.p)
    # Qs = on_architecture(Oceananigans.CPU(), surface_atmosphere_state.Qs)
    # Qℓ = on_architecture(Oceananigans.CPU(), surface_atmosphere_state.Qℓ)
    # Mp = on_architecture(Oceananigans.CPU(), interpolated_prescribed_freshwater_flux)
    # ρf = fluxes.freshwater_density

    # interpolator = interfaces.exchanger.atmosphere_exchanger.to_exchange_interp
    # exchange_atmosphere_state = interfaces.exchanger.exchange_atmosphere_state

    # ue  = parent(exchange_atmosphere_state.u)
    # ve  = parent(exchange_atmosphere_state.v)
    # Te  = parent(exchange_atmosphere_state.T)
    # qe  = parent(exchange_atmosphere_state.q)
    # pe  = parent(exchange_atmosphere_state.p)

    # ue = dropdims(ue, dims=3)
    # ve = dropdims(ve, dims=3)
    # Te = dropdims(Te, dims=3)
    # qe = dropdims(qe, dims=3)
    # pe = dropdims(pe, dims=3)

    return nothing
end

function atmosphere_exchanger(atmosphere::SpeedySimulation, exchange_grid)
    cpu_grid = on_architecture(CPU(), exchange_grid)
    cpu_surface_state = (
        u  = Field{Center, Center, Nothing}(cpu_grid),
        v  = Field{Center, Center, Nothing}(cpu_grid),
        T  = Field{Center, Center, Nothing}(cpu_grid),
        q  = Field{Center, Center, Nothing}(cpu_grid),
        p  = Field{Center, Center, Nothing}(cpu_grid),
        Qs = Field{Center, Center, Nothing}(cpu_grid),
        Qℓ = Field{Center, Center, Nothing}(cpu_grid),
    )

    # Figure this out:
    spectral_grid = atmosphere.model.spectral_grid
    # 1/4 degree?
    interpolator = SpeedyWeather.RingGrids.AnvilInterpolator(Float32,
        SpeedyWeather.FullClenshawGrid, 90, spectral_grid.npoints)

    arch = exchange_grid.architecture
    tmp_grid   = LatitudeLongitudeGrid(arch; size=(360, 179, 1), latitude=(-89.5, 89.5), longitude=(0, 360), z = (0, 1))
    londs, latds = SpeedyWeather.RingGrids.get_londlatds(spectral_grid.Grid, spectral_grid.nlat_half)
    SpeedyWeather.RingGrids.update_locator!(interpolator, londs, latds)

    exchanger = (; cpu_surface_state)

    return something
end

initialize!(::StateExchanger, ::SpeedySimulation) = nothing

function compute_net_atmosphere_fluxes!(coupled_model::SpeedyCoupledModel)
    return nothing
end
