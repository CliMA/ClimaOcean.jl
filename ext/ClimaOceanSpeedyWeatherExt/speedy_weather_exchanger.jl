using ConservativeRegridding
using GeometryOps: CutAtAntimeridianAndPoles, ClosedRing
using GeoInterface: Polygon, LinearRing
using Oceananigans

function atmosphere_exchanger(atmosphere::SpeedySimulation, exchange_grid, exchange_atmosphere_state)

    # Figure this out:
    spectral_grid = atmosphere.model.spectral_grid
    grid_faces = ConservativeRegridding.get_faces(spectral_grid)
    exchange_faces = list_cell_vertices(exchange_grid)
    arch = architecture(exchange_grid)

    if arch isa CPU # In case of a CPU grid, we reuse the already allocated fields
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
    
    regularizing_function = ClosedRing() ∘ CutAtAntimeridianAndPoles()

    # Magical incantation from ConservativeRegridding
    polys1 = map(regularizing_function, (Polygon([LinearRing(f)]) for f in eachcol(grid_faces)))
    polys2 = map(regularizing_function, (Polygon([LinearRing(f)]) for f in eachcol(exchange_faces)))

    regridder = ConservativeRegridding.intersection_areas(polys1, polys2)

    exchange_areas = vec(sum(A, dims=2))
    atmosphere_areas = vec(sum(A, dims=1))

    exchanger = (; cpu_surface_state, regridder, exchange_areas, atmosphere_areas)

    return exchanger
end

fill_exchange_fields!(::Oceananigans.CPU, args...) = nothing

function fill_exchange_fields!(::Oceananigans.GPU, exchange_state, cpu_surface_state)
    # Can I just copyto! here?
end

function interpolate_atmosphere_state!(interfaces, atmos::SpeedySimulation, coupled_model)
    # Get the atmospheric state on the ocean grid
    exchanger = interfaces.exchanger
    regridder = exchanger.regridder
    exchange_grid = interfaces.exchanger.exchange_grid
    exchange_state = exchanger.exchange_atmosphere_state
    arch = architecture(exchange_grid)

    ua = atmos.diagnostic_variables.grid.u_grid[:, end]           
    va = atmos.diagnostic_variables.grid.v_grid[:, end]           
    Ta = atmos.diagnostic_variables.grid.temp_grid[:, end]        
    qa = atmos.diagnostic_variables.grid.humid_grid[:, end]       
    pa = exp.(atmos.diagnostic_variables.grid.pres_grid[:, end])  
    Qsa = atmos.diagnostic_variables.physics.surface_shortwave_down
    Qla = atmos.diagnostic_variables.physics.surface_longwave_down 

    exchange_areas = exchanger.exchange_areas

    ue = exchanger.cpu_surface_state.u
    ve = exchanger.cpu_surface_state.v
    Te = exchanger.cpu_surface_state.T
    qe = exchanger.cpu_surface_state.q
    pe = exchanger.cpu_surface_state.p
    Qse = exchanger.cpu_surface_state.Qs
    Qle = exchanger.cpu_surface_state.Qℓ

    regrid!(ua, ue, regridder, exchange_areas)
    regrid!(va, ve, regridder, exchange_areas)
    regrid!(Ta, Te, regridder, exchange_areas)
    regrid!(qa, qe, regridder, exchange_areas)
    regrid!(pa, pe, regridder, exchange_areas)
    regrid!(Qsa, Qse, regridder, exchange_areas)
    regrid!(Qla, Qle, regridder, exchange_areas)
    
    fill_exchange_fields!(arch, exchange_state, exchanger.cpu_surface_state)

    return nothing
end