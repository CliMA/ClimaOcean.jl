using Oceananigans.Operators: extrinsic_vector, intrinsic_vector
using Oceananigans.Grids: _node

const c = Center()
const f = Face()

#####
##### Surface flux computation
#####

function compute_atmosphere_ocean_fluxes!(coupled_model)
    ocean = coupled_model.ocean
    atmosphere = coupled_model.atmosphere
    atmosphere_grid = atmosphere.grid
    sea_ice = coupled_model.sea_ice

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

    atmosphere_pressure = atmosphere.pressure.data

    Qs = atmosphere.downwelling_radiation.shortwave
    Qℓ = atmosphere.downwelling_radiation.longwave
    downwelling_radiation = (shortwave=Qs.data, longwave=Qℓ.data)

    freshwater_flux = map(ϕ -> ϕ.data, atmosphere.freshwater_flux)

    # Extract info for time-interpolation
    u = atmosphere.velocities.u # for example
    atmosphere_times = u.times
    atmosphere_backend = u.backend
    atmosphere_time_indexing = u.time_indexing

    # kernel parameters that compute fluxes in 0:Nx+1 and 0:Ny+1
    Nx, Ny, Nz = size(grid)
    single_column_grid = Nx == 1 && Ny == 1

    if single_column_grid
        kernel_parameters = KernelParameters(1:1, 1:1)
    else
        kernel_parameters = KernelParameters(0:Nx+1, 0:Ny+1)
    end

    surface_atmosphere_state = coupled_model.fluxes.surface_atmosphere_state

    # Simplify NamedTuple to reduce parameter space consumption.
    # See https://github.com/CliMA/ClimaOcean.jl/issues/116.
    surface_atmosphere_state = (u = surface_atmosphere_state.u.data,
                                v = surface_atmosphere_state.v.data,
                                T = surface_atmosphere_state.T.data,
                                p = surface_atmosphere_state.p.data,
                                q = surface_atmosphere_state.q.data,
                                Qs = surface_atmosphere_state.Qs.data,
                                Qℓ = surface_atmosphere_state.Qℓ.data,
                                Mp = surface_atmosphere_state.Mp.data)

    launch!(arch, grid, kernel_parameters,
            _interpolate_primary_atmospheric_state!,
            surface_atmosphere_state,
            grid,
            clock,
            atmosphere_velocities,
            atmosphere_tracers,
            atmosphere_pressure,
            downwelling_radiation,
            freshwater_flux,
            atmosphere_grid,
            atmosphere_times,
            atmosphere_backend,
            atmosphere_time_indexing)

    # Separately interpolate the auxiliary freshwater fluxes, which may
    # live on a different grid than the primary fluxes and atmospheric state.
    
    auxiliary_freshwater_flux = atmosphere.auxiliary_freshwater_flux
    interpolated_prescribed_freshwater_flux = surface_atmosphere_state.Mp

    if !isnothing(auxiliary_freshwater_flux)
        # TODO: do not assume that `auxiliary_freshater_flux` is a tuple
        auxiliary_data = map(ϕ -> ϕ.data, auxiliary_freshwater_flux)

        first_auxiliary_flux = first(auxiliary_freshwater_flux)
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

    #####
    ##### Next compute turbulent fluxes.
    #####

     # Fluxes, and flux contributors
    centered_velocity_fluxes = (u = coupled_model.fluxes.total.ocean.momentum.uᶜᶜᶜ,
                                v = coupled_model.fluxes.total.ocean.momentum.vᶜᶜᶜ)

    staggered_velocity_fluxes = (u = coupled_model.fluxes.total.ocean.momentum.u,
                                 v = coupled_model.fluxes.total.ocean.momentum.v)
    net_tracer_fluxes    = coupled_model.fluxes.total.ocean.tracers
    similarity_theory    = coupled_model.fluxes.turbulent.coefficients.ocean
    radiation_properties = coupled_model.fluxes.radiation
    turbulent_fluxes     = coupled_model.fluxes.turbulent.fields.ocean

    ocean_state = (u = ocean.model.velocities.u,
                   v = ocean.model.velocities.v,
                   T = ocean.model.tracers.T,
                   S = ocean.model.tracers.S)

    launch!(arch, grid, kernel_parameters,
            _compute_atmosphere_surface_similarity_theory_fluxes!,
            turbulent_fluxes,
            similarity_theory,
            grid,
            clock,
            ocean_state,
            AtmosphericThermodynamics.Liquid(), # surface phase
            coupled_model.fluxes.ocean_reference_density,
            coupled_model.fluxes.ocean_heat_capacity,
            coupled_model.fluxes.ocean_temperature_units,
            surface_atmosphere_state,
            radiation_properties,
            atmosphere.reference_height, # height at which the state is known
            atmosphere.boundary_layer_height,
            atmosphere.thermodynamics_parameters)   

    #####
    ##### Finally cobble together and properly interpolate fluxes
    ##### to be used by the ocean model.
    #####

    interpolated_downwelling_radiation = (shortwave = surface_atmosphere_state.Qs,
                                          longwave = surface_atmosphere_state.Qℓ)
    
    launch!(arch, grid, kernel_parameters,
            _assemble_atmosphere_ocean_fluxes!,
            centered_velocity_fluxes,
            net_tracer_fluxes,
            grid,
            clock,
            ocean_state.T,
            ocean_state.S,
            coupled_model.fluxes.ocean_temperature_units,
            fluxes_fields, 
            interpolated_downwelling_radiation,
            interpolated_prescribed_freshwater_flux,
            radiation_properties,
            coupled_model.fluxes.ocean_reference_density,
            coupled_model.fluxes.ocean_heat_capacity,
            coupled_model.fluxes.freshwater_density)
                
    limit_fluxes_over_sea_ice!(grid, kernel_parameters, sea_ice,
                               centered_velocity_fluxes,
                               net_tracer_fluxes,
                               ocean_state.T,
                               ocean_state.S)

    launch!(arch, grid, :xy, reconstruct_momentum_fluxes!,
            grid, staggered_velocity_fluxes, centered_velocity_fluxes)

    return nothing
end

@kernel function _interpolate_primary_atmospheric_state!(surface_atmos_state,
                                                         grid,
                                                         clock,
                                                         atmos_velocities,
                                                         atmos_tracers,
                                                         atmos_pressure,
                                                         downwelling_radiation,
                                                         prescribed_freshwater_flux,
                                                         atmos_grid,
                                                         atmos_times,
                                                         atmos_backend,
                                                         atmos_time_indexing)

    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3) # index of the top ocean cell

    @inbounds begin
        # Atmos state, which is _assumed_ to exist at location = (c, c, nothing)
        # The third index "k" should not matter but we put the correct index to get
        # a surface node anyways.
        atmos_args = (atmos_grid, atmos_times, atmos_backend, atmos_time_indexing)
        X = _node(i, j, kᴺ + 1, grid, c, c, f)
        time = Time(clock.time)

        uₐ = interp_atmos_time_series(atmos_velocities.u, X, time, atmos_args...)
        vₐ = interp_atmos_time_series(atmos_velocities.v, X, time, atmos_args...)
        Tₐ = interp_atmos_time_series(atmos_tracers.T,    X, time, atmos_args...)
        qₐ = interp_atmos_time_series(atmos_tracers.q,    X, time, atmos_args...)
        pₐ = interp_atmos_time_series(atmos_pressure,     X, time, atmos_args...)

        Qs = interp_atmos_time_series(downwelling_radiation.shortwave, X, time, atmos_args...)
        Qℓ = interp_atmos_time_series(downwelling_radiation.longwave,  X, time, atmos_args...)

        # Usually precipitation
        Mh = interp_atmos_time_series(prescribed_freshwater_flux, X, time, atmos_args...)

        # Convert atmosphere velocities (defined on a latitude-longitude grid) to 
        # the frame of reference of the native grid
        uₐ, vₐ = intrinsic_vector(i, j, kᴺ, grid, uₐ, vₐ)
    
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
                                                         grid,
                                                         clock,
                                                         auxiliary_freshwater_flux,
                                                         auxiliary_grid,
                                                         auxiliary_times,
                                                         auxiliary_backend,
                                                         auxiliary_time_indexing)

    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3) # index of the top ocean cell

    @inbounds begin
        X = _node(i, j, kᴺ + 1, grid, c, c, f)
        time = Time(clock.time)
        Mr = interp_atmos_time_series(auxiliary_freshwater_flux, X, time,
                                      auxiliary_grid,
                                      auxiliary_times,
                                      auxiliary_backend,
                                      auxiliary_time_indexing)

        freshwater_flux[i, j, 1] += Mr
    end
end

@kernel function _assemble_atmosphere_ocean_fluxes!(centered_velocity_fluxes,
                                                    net_tracer_fluxes,
                                                    grid,
                                                    clock,
                                                    ocean_temperature,
                                                    ocean_salinity,
                                                    ocean_temperature_units,
                                                    turbulent_fluxes,
                                                    downwelling_radiation,
                                                    prescribed_freshwater_flux,
                                                    radiation_properties,
                                                    ocean_reference_density,
                                                    ocean_heat_capacity,
                                                    freshwater_density)

    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3)
    time = Time(clock.time)

    @inbounds begin
        Tₒ = ocean_temperature[i, j, kᴺ]
        Tₒ = convert_to_kelvin(ocean_temperature_units, Tₒ)
        Sₒ = ocean_salinity[i, j, kᴺ]

        Qs = downwelling_radiation.shortwave[i, j, 1]
        Qℓ = downwelling_radiation.longwave[i, j, 1]

        Mp = prescribed_freshwater_flux[i, j, 1]

        Qc  = turbulent_fluxes.sensible_heat[i, j, 1] # sensible or "conductive" heat flux
        Qv  = turbulent_fluxes.latent_heat[i, j, 1]   # latent heat flux
        Mv  = turbulent_fluxes.water_vapor[i, j, 1]   # mass flux of water vapor
        ρτx = turbulent_fluxes.x_momentum[i, j, 1]    # zonal momentum flux
        ρτy = turbulent_fluxes.y_momentum[i, j, 1]    # meridional momentum flux
    end

    # Compute heat fluxes, bulk flux first
    Qd = net_downwelling_radiation(i, j, grid, time, radiation_properties, Qs, Qℓ)
    Qu = net_upwelling_radiation(i, j, grid, time, radiation_properties, Tₒ)

    ΣQ = Qd + Qu + Qc + Qv

    # Convert from a mass flux to a volume flux (aka velocity)
    # by dividing with the density of freshwater.
    # Also switch the sign, for some reason we are given freshwater flux as positive down.
    ρf⁻¹ = 1 / freshwater_density
    ΣF   = - Mp * ρf⁻¹

    # Add the contribution from the turbulent water vapor flux, which has
    # a different sign convention as the prescribed water mass fluxes (positive upwards)
    Fv = Mv * ρf⁻¹
    ΣF += Fv

    # Compute fluxes for u, v, T, S from momentum, heat, and freshwater fluxes
    τx = centered_velocity_fluxes.u
    τy = centered_velocity_fluxes.v
    Jᵀ = net_tracer_fluxes.T
    Jˢ = net_tracer_fluxes.S

    ρₒ⁻¹ = 1 / ocean_reference_density
    cₒ   = ocean_heat_capacity

    atmos_ocean_τx = ρτx * ρₒ⁻¹
    atmos_ocean_τy = ρτy * ρₒ⁻¹
    atmos_ocean_Jᵀ = ΣQ  * ρₒ⁻¹ / cₒ
    atmos_ocean_Jˢ = - Sₒ * ΣF

    # Mask fluxes over land for convenience
    inactive = inactive_node(i, j, kᴺ, grid, c, c, c)

    @inbounds begin
        τx[i, j, 1] = ifelse(inactive, 0, atmos_ocean_τx)
        τy[i, j, 1] = ifelse(inactive, 0, atmos_ocean_τy)
        Jᵀ[i, j, 1] = ifelse(inactive, 0, atmos_ocean_Jᵀ)
        Jˢ[i, j, 1] = ifelse(inactive, 0, atmos_ocean_Jˢ)
    end
end

@kernel function reconstruct_momentum_fluxes!(grid, J, Jᶜᶜᶜ)
    i, j = @index(Global, NTuple)

    @inbounds begin
        J.u[i, j, 1] = ℑxᶠᵃᵃ(i, j, 1, grid, Jᶜᶜᶜ.u) 
        J.v[i, j, 1] = ℑyᵃᶠᵃ(i, j, 1, grid, Jᶜᶜᶜ.v) 
    end
end

# Fallback for a `Nothing` radiation scheme
@inline   net_upwelling_radiation(i, j, grid, time, ::Nothing, Tₒ)     = zero(grid)
@inline net_downwelling_radiation(i, j, grid, time, ::Nothing, Qs, Qℓ) = zero(grid)

@inline function net_downwelling_radiation(i, j, grid, time, radiation, Qs, Qℓ)
    α = stateindex(radiation.reflection.ocean, i, j, 1, grid, time)
    ϵ = stateindex(radiation.emission.ocean, i, j, 1, grid, time)
    
    return @inbounds - (1 - α) * Qs - ϵ * Qℓ
end

@inline function net_upwelling_radiation(i, j, grid, time, radiation, Tₒ)
    σ = radiation.stefan_boltzmann_constant
    ϵ = stateindex(radiation.emission.ocean, i, j, 1, grid, time)

    # Note: positive implies _upward_ heat flux, and therefore cooling.
    return ϵ * σ * Tₒ^4
end
