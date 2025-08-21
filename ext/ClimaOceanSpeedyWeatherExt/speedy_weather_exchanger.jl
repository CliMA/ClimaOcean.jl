using Oceananigans
using Oceananigans.Grids: architecture
using Oceananigans.Utils: launch!
using Oceananigans.Operators: intrinsic_vector

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
    arch = architecture(exchange_grid)
    FT = eltype(exchange_atmosphere_state.u)

    if arch isa Oceananigans.CPU # In case of a CPU grid, we reuse the already allocated fields
        cpu_surface_state = (
            u  = vec(interior(exchange_atmosphere_state.u)),
            v  = vec(interior(exchange_atmosphere_state.v)),
            T  = vec(interior(exchange_atmosphere_state.T)),
            q  = vec(interior(exchange_atmosphere_state.q)),
            p  = vec(interior(exchange_atmosphere_state.p)),
            Qs = vec(interior(exchange_atmosphere_state.Qs)),
            Qℓ = vec(interior(exchange_atmosphere_state.Qℓ)),
            Mp = vec(interior(exchange_atmosphere_state.Mp))
        )
    else # Otherwise we allocate new CPU fields
        cpu_surface_state = (
            u  = zeros(FT, size(exchange_grid, 1) * size(exchange_grid, 2)),
            v  = zeros(FT, size(exchange_grid, 1) * size(exchange_grid, 2)),
            T  = zeros(FT, size(exchange_grid, 1) * size(exchange_grid, 2)),
            q  = zeros(FT, size(exchange_grid, 1) * size(exchange_grid, 2)),
            p  = zeros(FT, size(exchange_grid, 1) * size(exchange_grid, 2)),
            Qs = zeros(FT, size(exchange_grid, 1) * size(exchange_grid, 2)),
            Qℓ = zeros(FT, size(exchange_grid, 1) * size(exchange_grid, 2)),
            Mp = zeros(FT, size(exchange_grid, 1) * size(exchange_grid, 2))
        )
    end
    
    # TODO: Implement a conservative regridder when ready
    regridder = BilinearInterpolator(exchange_grid, spectral_grid)
    exchanger = (; cpu_surface_state, regridder)

    return exchanger
end

fill_exchange_fields!(::Oceananigans.CPU, args...) = nothing

# TODO: improve GPU support
function fill_exchange_fields!(::Oceananigans.GPU, state, cpustate)
    set!(state.u,  reshape(cpustate.u,  size(state.u)))
    set!(state.v,  reshape(cpustate.v,  size(state.v)))
    set!(state.T,  reshape(cpustate.T,  size(state.T)))
    set!(state.q,  reshape(cpustate.q,  size(state.q)))
    set!(state.p,  reshape(cpustate.p,  size(state.p)))
    set!(state.Qs, reshape(cpustate.Qs, size(state.Qs)))
    set!(state.Qℓ, reshape(cpustate.Qℓ, size(state.Qℓ)))
    set!(state.Mp, reshape(cpustate.Mp, size(state.Mp)))
    return nothing
end

# Regrid the atmospheric state on the exchange grid
function interpolate_atmosphere_state!(interfaces, atmos::SpeedySimulation, coupled_model)
    atmosphere_exchanger = interfaces.exchanger.atmosphere_exchanger
    regridder = atmosphere_exchanger.regridder
    exchange_grid = interfaces.exchanger.exchange_grid
    exchange_state = interfaces.exchanger.exchange_atmosphere_state
    surface_layer = atmos.model.spectral_grid.nlayers

    ua  = RingGrids.field_view(atmos.diagnostic_variables.grid.u_grid, :, surface_layer)
    va  = RingGrids.field_view(atmos.diagnostic_variables.grid.v_grid, :, surface_layer)
    Ta  = RingGrids.field_view(atmos.diagnostic_variables.grid.temp_grid, :, surface_layer)
    qa  = RingGrids.field_view(atmos.diagnostic_variables.grid.humid_grid, :, surface_layer)
    pa  = exp.(atmos.diagnostic_variables.grid.pres_grid[:, end])
    Qsa = atmos.diagnostic_variables.physics.surface_shortwave_down
    Qla = atmos.diagnostic_variables.physics.surface_longwave_down
    Mpa = atmos.diagnostic_variables.physics.total_precipitation_rate

    ue  = atmosphere_exchanger.cpu_surface_state.u
    ve  = atmosphere_exchanger.cpu_surface_state.v
    Te  = atmosphere_exchanger.cpu_surface_state.T
    qe  = atmosphere_exchanger.cpu_surface_state.q
    pe  = atmosphere_exchanger.cpu_surface_state.p
    Qse = atmosphere_exchanger.cpu_surface_state.Qs
    Qle = atmosphere_exchanger.cpu_surface_state.Qℓ
    Mpe = atmosphere_exchanger.cpu_surface_state.Mp

    regrid!(ue,  regridder.set1, ua)
    regrid!(ve,  regridder.set1, va)
    regrid!(Te,  regridder.set1, Ta)
    regrid!(qe,  regridder.set1, qa)
    regrid!(pe,  regridder.set1, pa)
    regrid!(Qse, regridder.set1, Qsa)
    regrid!(Qle, regridder.set1, Qla)
    regrid!(Mpe, regridder.set1, Mpa)
    
    arch = architecture(exchange_grid)
    fill_exchange_fields!(arch, exchange_state, atmosphere_exchanger.cpu_surface_state)

    u = exchange_state.u
    v = exchange_state.v

    launch!(arch, exchange_grid, :xy, _rotate_winds!, u, exchange_grid, v)

    return nothing
end

@kernel function _rotate_winds!(u, grid, v)
    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3) 
    uₑ, vₑ = intrinsic_vector(i, j, kᴺ, grid, u, v)
    @inbounds u[i, j, kᴺ] = uₑ
    @inbounds v[i, j, kᴺ] = vₑ
end

# TODO: For the moment this is just ciupling between ocean and atmosphere.
# we will also need to add the coupling with the sea-ice model
function compute_net_atmosphere_fluxes!(coupled_model::SpeedyCoupledModel)
    atmos = coupled_model.atmosphere
    regridder = coupled_model.interfaces.exchanger.atmosphere_exchanger.regridder

    # All the location of these fluxes will change
    Qc = coupled_model.interfaces.atmosphere_ocean_interface.fluxes.sensible_heat
    Mv = coupled_model.interfaces.atmosphere_ocean_interface.fluxes.water_vapor

    Qca = atmos.prognostic_variables.ocean.sensible_heat_flux
    Mva = atmos.prognostic_variables.ocean.surface_humidity_flux

    # TODO: Figure out how we are going to deal with upwelling radiation
    regrid!(Qca, regridder.set2, interior(Qc))
    regrid!(Mva, regridder.set2, interior(Mv))

    return nothing
end
