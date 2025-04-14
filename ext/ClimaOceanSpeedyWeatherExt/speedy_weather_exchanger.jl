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

# For the moment the workflow is:
# 1. Perform the regridding on the CPU
# 2. Eventually copy the regridded fields to the GPU
# If this work we can
# 1. Copy speedyweather gridarrays to the GPU
# 2. Perform the regridding on the GPU
function atmosphere_exchanger(atmosphere::SpeedySimulation, exchange_grid, exchange_atmosphere_state)

    # Figure this out:
    spectral_grid = atmosphere.model.spectral_grid
    grid_faces = get_faces(spectral_grid)
    exchange_faces = get_faces(exchange_grid)
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
            u  = zeros(Float32, length(exchange_faces)),
            v  = zeros(Float32, length(exchange_faces)),
            T  = zeros(Float32, length(exchange_faces)),
            q  = zeros(Float32, length(exchange_faces)),
            p  = zeros(Float32, length(exchange_faces)),
            Qs = zeros(Float32, length(exchange_faces)),
            Qℓ = zeros(Float32, length(exchange_faces)),
        )
    end
    
    # Magical incantation from ConservativeRegridding
    polys1 = map(ClosedRing() ∘ CutAtAntimeridianAndPoles(), (Polygon([LinearRing(f)]) for f in eachcol(grid_faces)))
    polys2 = map(ClosedRing() ∘ CutAtAntimeridianAndPoles(), (Polygon([LinearRing(f)]) for f in eachcol(exchange_faces)))

    regridder = ConservativeRegridding.intersection_areas(polys1, polys2)

    exchange_areas   = vec(sum(A, dims=2))
    atmosphere_areas = vec(sum(A, dims=1))

    exchanger = (; cpu_surface_state, regridder, exchange_areas, atmosphere_areas)

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

    exchange_areas = exchanger.exchange_areas

    ue  = exchanger.cpu_surface_state.u
    ve  = exchanger.cpu_surface_state.v
    Te  = exchanger.cpu_surface_state.T
    qe  = exchanger.cpu_surface_state.q
    pe  = exchanger.cpu_surface_state.p
    Qse = exchanger.cpu_surface_state.Qs
    Qle = exchanger.cpu_surface_state.Qℓ

    ConservativeRegridding.regrid!(ua,  ue,  regridder, exchange_areas)
    ConservativeRegridding.regrid!(va,  ve,  regridder, exchange_areas)
    ConservativeRegridding.regrid!(Ta,  Te,  regridder, exchange_areas)
    ConservativeRegridding.regrid!(qa,  qe,  regridder, exchange_areas)
    ConservativeRegridding.regrid!(pa,  pe,  regridder, exchange_areas)
    ConservativeRegridding.regrid!(Qsa, Qse, regridder, exchange_areas)
    ConservativeRegridding.regrid!(Qla, Qle, regridder, exchange_areas)
    
    arch = architecture(exchange_grid)
    fill_exchange_fields!(arch, exchange_state, exchanger.cpu_surface_state)

    return nothing
end