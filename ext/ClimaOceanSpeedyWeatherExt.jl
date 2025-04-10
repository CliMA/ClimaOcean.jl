module ClimaOceanSpeedyWeatherExt

using OffsetArrays
using KernelAbstractions
using Statistics

import SpeedyWeather as SW
import ClimaOcean as CO
import Oceananigans as OC
import CUDA

# This document lays out the functions that must be extended to
# use an atmospheric simulation in ClimaOcean.

import Oceananigans.TimeSteppers: time_step!
import Oceananigans.Models: update_model_field_time_series!

# Make sure the atmospheric parameters from SpeedyWeather can be used in the compute fluxes function
import ClimaOcean.OceanSeaIceModels.PrescribedAtmospheres: 
    thermodynamics_parameters, 
    boundary_layer_height, 
    surface_layer_height

import ClimaOcean.OceanSeaIceModels:
    compute_net_atmosphere_fluxes!

import ClimaOcean.OceanSeaIceModels.InterfaceComputations:
    atmosphere_exchanger,
    initialize!,
    StateExchanger,
    interpolate_atmosphere_state!

using ClimaOcean.OceanSeaIceModels: OceanSeaIceModel

const SpeedySimulation = SpeedyWeather.Simulation
const ClimaCoupledModel = OceanSeaIceModel{<:Any, <:SpeedySimulation}
Base.summary(::SpeedySimulation) = "SpeedyWeather.Simulation"

# This can be left blank:
update_model_field_time_series!(::SpeedySimulation, time) = nothing

# Take one time-step
function time_step!(atmos::SpeedySimulation, Δt)
    # TODO: check if the time-step can be changed.
    @assert Δt == atmos.integrator.dt
    CA.SciMLBase.step!(atmos.integrator)
    return nothing
end

# The height of near-surface variables used in the turbulent flux solver
function surface_layer_height(s::SpeedySimulation)
    T = s.model.atmosphere.temp_ref
    g = s.model.planet.gravity
    Φ = s.model.geopotential.Δp_geopot_full
    return Φ[end] * T / g
end

# This is a parameter that is used in the computation of the fluxes,
# It probably should not be here but in the similarity theory type.
boundary_layer_height(atmos::SpeedySimulation) = 600

Base.eltype(::EarthAtmosphere{FT}) where FT = FT

# This is a _hack_!! The parameters should be consistent with what is specified in SpeedyWeather
thermodynamics_parameters(atmos::SpeedyWeather.Simulation) = 
    PrescribedAtmosphereThermodynamicsParameters(eltype(atmos.model.atmosphere))

using Oceananigans.Grids: λnodes, φnodes, LatitudeLongitudeGrid
using Oceananigans.Fields: Center
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

    # RingGrids.interpolate!(vec(view(ua, :, :, 1)), atmos.diagnostic_variables.grid.u_grid[:, end],            interpolator)
    # RingGrids.interpolate!(vec(view(va, :, :, 1)), atmos.diagnostic_variables.grid.v_grid[:, end],            interpolator)
    # RingGrids.interpolate!(vec(view(Ta, :, :, 1)), atmos.diagnostic_variables.grid.temp_grid[:, end],         interpolator)
    # RingGrids.interpolate!(vec(view(qa, :, :, 1)), atmos.diagnostic_variables.grid.humid_grid[:, end],        interpolator)
    # RingGrids.interpolate!(vec(view(pa, :, :, 1)), exp.(atmos.diagnostic_variables.grid.pres_grid[:, end]),   interpolator)
    # RingGrids.interpolate!(vec(view(Qs, :, :, 1)), atmos.diagnostic_variables.physics.surface_shortwave_down, interpolator)
    # RingGrids.interpolate!(vec(view(Qℓ, :, :, 1)), atmos.diagnostic_variables.physics.surface_longwave_down,  interpolator)
    # RingGrids.interpolate!(vec(view(Mp, :, :, 1)), atmosphere_precipitation,                                  interpolator)

    interpolator = interfaces.exchanger.atmosphere_exchanger.to_exchange_interp
    exchange_atmosphere_state = interfaces.exchanger.exchange_atmosphere_state

    ue  = parent(exchange_atmosphere_state.u)
    ve  = parent(exchange_atmosphere_state.v)
    Te  = parent(exchange_atmosphere_state.T)
    qe  = parent(exchange_atmosphere_state.q)
    pe  = parent(exchange_atmosphere_state.p)

    ue = dropdims(ue, dims=3)
    ve = dropdims(ve, dims=3)
    Te = dropdims(Te, dims=3)
    qe = dropdims(qe, dims=3)
    pe = dropdims(pe, dims=3)

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

    exchanger = (; cpu_surface_state)

    return something
end

initialize!(::StateExchanger, ::SpeedySimulation) = nothing

function compute_net_atmosphere_fluxes!(coupled_model::SpeedyCoupledModel)
    return nothing
end

end # module ClimaOceanSpeedyWeatherExt