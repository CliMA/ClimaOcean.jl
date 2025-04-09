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

# const ClimaCoupledModel = OceanSeaIceModel{<:Any, <:SpeedySimulation}

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
surface_layer_height(s::SpeedySimulation) = 10 # meters, for example

# This is a parameter that is used in the computation of the fluxes,
# It probably should not be here but in the similarity theory type.
boundary_layer_height(atmos::SpeedySimulation) = 600

# Note: possibly, can use the atmos thermodynamic parameters directly here.
thermodynamics_parameters(atmos::SpeedySimulation) = 
    atmos.integrator.p.params.thermodynamics_params

Base.summary(::SpeedySimulation) = "ClimaAtmos.AtmosSimulation"

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
    return something
end

initialize!(::StateExchanger, ::SpeedySimulation) = nothing

function compute_net_atmosphere_fluxes!(coupled_model::SpeedyCoupledModel)
    return nothing
end

end # module ClimaOceanSpeedyWeatherExt