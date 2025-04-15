using ConservativeRegridding
using GeometryOps: CutAtAntimeridianAndPoles, ClosedRing
using GeoInterface: Polygon, LinearRing
using Oceananigans
using Oceananigans.Grids: architecture

import ClimaOcean.OceanSeaIceModels:
    compute_net_atmosphere_fluxes!

import ClimaOcean.OceanSeaIceModels.InterfaceComputations:
    atmosphere_exchanger,
    initialize!,
    StateExchanger,
    interpolate_atmosphere_state!

const OCRExt = Base.get_extension(Oceananigans, :OceananigansConservativeRegriddingExt)
const SWGExt = Base.get_extension(SpeedyWeather, :SpeedyWeatherGeoMakieExt)

get_cell_matrix(grid::SpeedyWeather.SpectralGrid) = SWGExt.get_faces(grid; add_nans=false)
get_cell_matrix(grid::Oceananigans.AbstractGrid)  = OCRExt.compute_cell_matrix(grid, Center(), Center(), nothing)

# For the moment the workflow is:
# 1. Perform the regridding on the CPU
# 2. Eventually copy the regridded fields to the GPU
# If this work we can
# 1. Copy speedyweather gridarrays to the GPU
# 2. Perform the regridding on the GPU
function atmosphere_exchanger(atmosphere::SpeedySimulation, exchange_grid, exchange_atmosphere_state)

    # Figure this out:
    spectral_grid = atmosphere.model.spectral_grid
    atmosphere_cell_matrix = get_cell_matrix(spectral_grid)
    exchange_cell_matrix = get_cell_matrix(exchange_grid)
    arch = architecture(exchange_grid)

    if arch isa Oceananigans.CPU # In case of a CPU grid, we reuse the already allocated fields
        cpu_surface_state = (
            u  = vec(interior(exchange_atmosphere_state.u)),
            v  = vec(interior(exchange_atmosphere_state.v)),
            T  = vec(interior(exchange_atmosphere_state.T)),
            q  = vec(interior(exchange_atmosphere_state.q)),
            p  = vec(interior(exchange_atmosphere_state.p)),
            Qs = vec(interior(exchange_atmosphere_state.Qs)),
            Qℓ = vec(interior(exchange_atmosphere_state.Qℓ))
        )
    else # Otherwise we allocate new CPU fields
        cpu_surface_state = (
            u  = zeros(Float32, length(exchange_cell_matrix)),
            v  = zeros(Float32, length(exchange_cell_matrix)),
            T  = zeros(Float32, length(exchange_cell_matrix)),
            q  = zeros(Float32, length(exchange_cell_matrix)),
            p  = zeros(Float32, length(exchange_cell_matrix)),
            Qs = zeros(Float32, length(exchange_cell_matrix)),
            Qℓ = zeros(Float32, length(exchange_cell_matrix)),
        )
    end
    
    regridder = ConservativeRegridding.Regridder(exchange_cell_matrix, atmosphere_cell_matrix)

    exchanger = (; cpu_surface_state, regridder)

    return exchanger
end

fill_exchange_fields!(::Oceananigans.CPU, args...) = nothing

# TODO: add GPU support
function fill_exchange_fields!(::Oceananigans.GPU, exchange_state, cpu_surface_state)
    # Can I just copyto! here?
end

# Regrid the atmospheric state on the exchange grid
function interpolate_atmosphere_state!(interfaces, atmos::SpeedySimulation, coupled_model)
    exchanger = interfaces.exchanger
    regridder = exchanger.regridder
    exchange_grid = interfaces.exchanger.exchange_grid
    exchange_state = exchanger.exchange_atmosphere_state

    ua  = atmos.diagnostic_variables.grid.u_grid[:, end]           
    va  = atmos.diagnostic_variables.grid.v_grid[:, end]           
    Ta  = atmos.diagnostic_variables.grid.temp_grid[:, end]        
    qa  = atmos.diagnostic_variables.grid.humid_grid[:, end]       
    pa  = exp.(atmos.diagnostic_variables.grid.pres_grid[:, end])  
    Qsa = atmos.diagnostic_variables.physics.surface_shortwave_down
    Qla = atmos.diagnostic_variables.physics.surface_longwave_down 

    ue  = exchanger.cpu_surface_state.u
    ve  = exchanger.cpu_surface_state.v
    Te  = exchanger.cpu_surface_state.T
    qe  = exchanger.cpu_surface_state.q
    pe  = exchanger.cpu_surface_state.p
    Qse = exchanger.cpu_surface_state.Qs
    Qle = exchanger.cpu_surface_state.Qℓ

    ConservativeRegridding.regrid!(ue,  regridder, ua)
    ConservativeRegridding.regrid!(ve,  regridder, va)
    ConservativeRegridding.regrid!(Te,  regridder, Ta)
    ConservativeRegridding.regrid!(qe,  regridder, qa)
    ConservativeRegridding.regrid!(pe,  regridder, pa)
    ConservativeRegridding.regrid!(Qse, regridder, Qsa)
    ConservativeRegridding.regrid!(Qle, regridder, Qla)
    
    arch = architecture(exchange_grid)
    fill_exchange_fields!(arch, exchange_state, exchanger.cpu_surface_state)

    return nothing
end

function compute_net_atmosphere_fluxes!(coupled_model::SpeedyCoupledModel)
    atmos = coupled_model.atmosphere

    # All the location of these fluxes will change
    Qs = similarity_theory_fields.sensible_heat
    Mv = similarity_theory_fields.water_vapor
    Qu = total_fluxes.ocean.heat.upwelling_radiation 

    regridder = coupled_model.interfaces.exchanger.regridder

    ConservativeRegridding.regrid!(atmos.diagnostic_variables.physics.sensible_heat_flux,  transpose(regridder), vec(interior(Qs)))
    ConservativeRegridding.regrid!(atmos.diagnostic_variables.physics.evaporative_flux,    transpose(regridder), vec(interior(Mv)))
    ConservativeRegridding.regrid!(atmos.diagnostic_variables.physics.surface_longwave_up, transpose(regridder), vec(interior(Qu)))

    return nothing
end
