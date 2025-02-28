using Oceananigans.Operators: intrinsic_vector
using Oceananigans.Grids: _node
using Oceananigans.OutputReaders: InterpolatingTimeIndices

using ...OceanSimulations: forcing_barotropic_potential

"""Interpolate the atmospheric state onto the ocean / sea-ice grid."""
function interpolate_atmospheric_state!(coupled_model)
    ocean = coupled_model.ocean
    atmosphere = coupled_model.atmosphere
    atmosphere_grid = atmosphere.grid

    # Basic model properties
    grid = ocean.model.grid
    arch = architecture(grid)
    clock = coupled_model.clock

    #####
    ##### First interpolate atmosphere time series
    ##### in time and to the ocean grid.
    #####
    
    # We use .data here to save parameter space (unlike Field, adapt_structure for
    # fts = FieldTimeSeries does not return fts.data)
    atmosphere_velocities = (u = atmosphere.velocities.u.data,
                             v = atmosphere.velocities.v.data)

    atmosphere_tracers = (T = atmosphere.tracers.T.data,
                          q = atmosphere.tracers.q.data)

    Qs = atmosphere.downwelling_radiation.shortwave
    Qℓ = atmosphere.downwelling_radiation.longwave
    downwelling_radiation = (shortwave=Qs.data, longwave=Qℓ.data)
    freshwater_flux = map(ϕ -> ϕ.data, atmosphere.freshwater_flux)
    atmosphere_pressure = atmosphere.pressure.data

    # Extract info for time-interpolation
    u = atmosphere.velocities.u # for example
    atmosphere_times = u.times
    atmosphere_backend = u.backend
    atmosphere_time_indexing = u.time_indexing

    atmosphere_fields = coupled_model.interfaces.exchanger.exchange_atmosphere_state
    space_interp_indices = coupled_model.interfaces.exchanger.atmosphere_exchanger
    exchange_grid = coupled_model.interfaces.exchanger.exchange_grid

    # Simplify NamedTuple to reduce parameter space consumption.
    # See https://github.com/CliMA/ClimaOcean.jl/issues/116.
    atmosphere_data = (u = atmosphere_fields.u.data,
                       v = atmosphere_fields.v.data,
                       T = atmosphere_fields.T.data,
                       p = atmosphere_fields.p.data,
                       q = atmosphere_fields.q.data,
                       Qs = atmosphere_fields.Qs.data,
                       Qℓ = atmosphere_fields.Qℓ.data,
                       Mp = atmosphere_fields.Mp.data)

    kernel_parameters = interface_kernel_parameters(grid)

    ua = atmosphere.velocities.u
    time_interp_indices = InterpolatingTimeIndices(ua, clock.time)
    
    launch!(arch, grid, kernel_parameters,
            _interpolate_primary_atmospheric_state!,
            atmosphere_data,
            space_interp_indices,
            time_interp_indices,
            exchange_grid,
            atmosphere_velocities,
            atmosphere_tracers,
            atmosphere_pressure,
            downwelling_radiation,
            freshwater_flux,
            atmosphere_backend,
            atmosphere_time_indexing)

    # Separately interpolate the auxiliary freshwater fluxes, which may
    # live on a different grid than the primary fluxes and atmospheric state.
    auxiliary_freshwater_flux = atmosphere.auxiliary_freshwater_flux
    interpolated_prescribed_freshwater_flux = atmosphere_data.Mp

    if !isnothing(auxiliary_freshwater_flux)
        # TODO: do not assume that `auxiliary_freshater_flux` is a tuple
        auxiliary_data = map(ϕ -> ϕ.data, auxiliary_freshwater_flux)

        first_auxiliary_flux    = first(auxiliary_freshwater_flux)
        auxiliary_grid          = first_auxiliary_flux.grid
        auxiliary_times         = first_auxiliary_flux.times
        auxiliary_backend       = first_auxiliary_flux.backend
        auxiliary_time_indexing = first_auxiliary_flux.time_indexing

        launch!(arch, grid, kernel_parameters,
                _interpolate_auxiliary_freshwater_flux!,
                interpolated_prescribed_freshwater_flux,
                grid,
                clock,
                auxiliary_data,
                auxiliary_grid,
                auxiliary_times,
                auxiliary_backend,
                auxiliary_time_indexing)
    end

    # Set ocean barotropic pressure forcing
    barotropic_potential = forcing_barotropic_potential(ocean)
    ρₒ = coupled_model.interfaces.ocean_properties.reference_density
    if !isnothing(barotropic_potential)
        parent(barotropic_potential) .= parent(atmosphere_data.p) ./ ρₒ
    end
end
    
@kernel function _interpolate_primary_atmospheric_state!(surface_atmos_state,
                                                         space_interp_indices,
                                                         time_interp_indices,
                                                         exchange_grid,
                                                         atmos_velocities,
                                                         atmos_tracers,
                                                         atmos_pressure,
                                                         downwelling_radiation,
                                                         prescribed_freshwater_flux,
                                                         atmos_backend,
                                                         atmos_time_indexing)

    i, j = @index(Global, NTuple)

    @inbounds begin
        x_itp = space_interp_indices[i, j, 1]
        t_itp = time_interp_indices
        atmos_args = (x_itp, t_itp, atmos_backend, atmos_time_indexing)

        uₐ = interp_atmos_time_series(atmos_velocities.u, atmos_args...)
        vₐ = interp_atmos_time_series(atmos_velocities.v, atmos_args...)
        Tₐ = interp_atmos_time_series(atmos_tracers.T,    atmos_args...)
        qₐ = interp_atmos_time_series(atmos_tracers.q,    atmos_args...)
        pₐ = interp_atmos_time_series(atmos_pressure,     atmos_args...)

        Qs = interp_atmos_time_series(downwelling_radiation.shortwave, atmos_args...)
        Qℓ = interp_atmos_time_series(downwelling_radiation.longwave,  atmos_args...)

        # Usually precipitation
        Mh = interp_atmos_time_series(prescribed_freshwater_flux, atmos_args...)

        # Convert atmosphere velocities (usually defined on a latitude-longitude grid) to 
        # the frame of reference of the native grid
        kᴺ = size(exchange_grid, 3) # index of the top ocean cell
        uₐ, vₐ = intrinsic_vector(i, j, kᴺ, exchange_grid, uₐ, vₐ)

        surface_atmos_state.u[i, j, 1] = uₐ
        surface_atmos_state.v[i, j, 1] = vₐ
        surface_atmos_state.T[i, j, 1] = Tₐ
        surface_atmos_state.p[i, j, 1] = pₐ
        surface_atmos_state.q[i, j, 1] = qₐ
        surface_atmos_state.Qs[i, j, 1] = Qs
        surface_atmos_state.Qℓ[i, j, 1] = Qℓ
        surface_atmos_state.Mp[i, j, 1] = Mh
    end
end

@kernel function _interpolate_auxiliary_freshwater_flux!(freshwater_flux,
                                                         interface_grid,
                                                         clock,
                                                         auxiliary_freshwater_flux,
                                                         auxiliary_grid,
                                                         auxiliary_times,
                                                         auxiliary_backend,
                                                         auxiliary_time_indexing)

    i, j = @index(Global, NTuple)
    kᴺ = size(interface_grid, 3) # index of the top ocean cell

    @inbounds begin
        X = _node(i, j, kᴺ + 1, interface_grid, c, c, f)
        time = Time(clock.time)
        Mr = interp_atmos_time_series(auxiliary_freshwater_flux, X, time,
                                      auxiliary_grid,
                                      auxiliary_times,
                                      auxiliary_backend,
                                      auxiliary_time_indexing)

        freshwater_flux[i, j, 1] += Mr
    end
end

#####
##### Utility for interpolating tuples of fields
#####

# Note: assumes loc = (c, c, nothing) (and the third location should not matter.)
@inline interp_atmos_time_series(J::AbstractArray, x_itp::FractionalIndices, t_itp, args...) =
    interpolate(x_itp, t_itp, J, args...)

@inline interp_atmos_time_series(J::AbstractArray, X, time, grid, args...) =
    interpolate(X, time, J, (c, c, nothing), grid, args...)

@inline interp_atmos_time_series(ΣJ::NamedTuple, args...) =
    interp_atmos_time_series(values(ΣJ), args...)

@inline interp_atmos_time_series(ΣJ::Tuple{<:Any}, args...) =
    interp_atmos_time_series(ΣJ[1], args...) +
    interp_atmos_time_series(ΣJ[2], args...)

@inline interp_atmos_time_series(ΣJ::Tuple{<:Any, <:Any}, args...) =
    interp_atmos_time_series(ΣJ[1], args...) +
    interp_atmos_time_series(ΣJ[2], args...)

@inline interp_atmos_time_series(ΣJ::Tuple{<:Any, <:Any, <:Any}, args...) =
    interp_atmos_time_series(ΣJ[1], args...) +
    interp_atmos_time_series(ΣJ[2], args...) +
    interp_atmos_time_series(ΣJ[3], args...)

@inline interp_atmos_time_series(ΣJ::Tuple{<:Any, <:Any, <:Any, <:Any}, args...) =
    interp_atmos_time_series(ΣJ[1], args...) +
    interp_atmos_time_series(ΣJ[2], args...) +
    interp_atmos_time_series(ΣJ[3], args...) +
    interp_atmos_time_series(ΣJ[4], args...)

