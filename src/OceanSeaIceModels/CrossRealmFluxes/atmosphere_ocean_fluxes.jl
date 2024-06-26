
#####
##### Surface flux computation
#####

const c = Center()
const f = Face()

function compute_atmosphere_ocean_fluxes!(coupled_model)
    ocean = coupled_model.ocean
    atmosphere = coupled_model.atmosphere
    atmosphere_grid = atmosphere.grid
    sea_ice = coupled_model.sea_ice

    # Basic model properties
    grid = ocean.model.grid
    arch = architecture(grid)
    clock = coupled_model.clock

    # Ocean, atmosphere, and sea ice state
    ocean_velocities  = surface_velocities(ocean)
    ocean_tracers     = surface_tracers(ocean)

    # Fluxes, and flux contributors
    centered_velocity_fluxes = (u = coupled_model.fluxes.total.ocean.momentum.uᶜᶜᶜ,
                                v = coupled_model.fluxes.total.ocean.momentum.vᶜᶜᶜ)

    staggered_velocity_fluxes = (u = coupled_model.fluxes.total.ocean.momentum.u,
                                 v = coupled_model.fluxes.total.ocean.momentum.v)

    net_tracer_fluxes    = coupled_model.fluxes.total.ocean.tracers
    similarity_theory    = coupled_model.fluxes.turbulent
    radiation_properties = coupled_model.fluxes.radiation

    ocean_state = merge(ocean_velocities, ocean_tracers)

    atmosphere_velocities = map(u -> u.data, atmosphere.velocities)
    atmosphere_tracers    = map(c -> c.data, atmosphere.tracers)
    atmosphere_pressure   = atmosphere.pressure.data

    atmosphere_state = merge(atmosphere_velocities, atmosphere_tracers, (; p=atmosphere_pressure))
    freshwater_flux  = map(ϕ -> ϕ.data, atmosphere.freshwater_flux)

    u = atmosphere.velocities.u # for example 
    atmosphere_times = u.times
    atmosphere_backend = u.backend
    atmosphere_time_indexing = u.time_indexing

    runoff_args = get_runoff_args(atmosphere.runoff_flux)

    Qs = atmosphere.downwelling_radiation.shortwave
    Ql = atmosphere.downwelling_radiation.longwave

    downwelling_radiation = (shortwave=Qs.data, longwave=Ql.data)

    kernel_size = (size(grid, 1) + 2, size(grid, 2) + 2)

    # kernel parameters that compute fluxes in 0:Nx+1 and 0:Ny+1
    kernel_parameters = KernelParameters(kernel_size, (-1, -1))

    launch!(arch, grid, kernel_parameters, _compute_atmosphere_ocean_similarity_theory_fluxes!,
            similarity_theory,
            grid,
            clock,
            ocean_state,
            coupled_model.fluxes.ocean_temperature_units,
            atmosphere_state,
            atmosphere_grid,
            atmosphere_times,
            atmosphere_backend,
            atmosphere_time_indexing,
            atmosphere.reference_height, # height at which the state is known
            atmosphere.boundary_layer_height,
            atmosphere.thermodynamics_parameters)   
    
    launch!(arch, grid, kernel_parameters, _assemble_atmosphere_ocean_fluxes!,
            centered_velocity_fluxes,
            net_tracer_fluxes,
            grid,
            clock,
            ocean_state.T,
            ocean_state.S,
            coupled_model.fluxes.ocean_temperature_units,
            similarity_theory.fields,
            downwelling_radiation,
            freshwater_flux,
            atmosphere_grid,
            atmosphere_times,
            atmosphere_backend,
            atmosphere_time_indexing,
            runoff_args,
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

# Fallback
@inline convert_to_intrinsic_reference_frame(i, j, k, grid, uₒ, vₒ) = uₒ, vₒ
@inline convert_to_extrinsic_reference_frame(i, j, k, grid, uₒ, vₒ) = uₒ, vₒ

# Fallback!
limit_fluxes_over_sea_ice!(args...) = nothing

@kernel function _compute_atmosphere_ocean_similarity_theory_fluxes!(similarity_theory,
                                                                     grid,
                                                                     clock,
                                                                     ocean_state,
                                                                     ocean_temperature_units,
                                                                     atmos_state,
                                                                     atmos_grid,
                                                                     atmos_times,
                                                                     atmos_backend,
                                                                     atmos_time_indexing,
                                                                     atmosphere_reference_height,
                                                                     atmosphere_boundary_layer_height,
                                                                     atmos_thermodynamics_parameters)

    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3)

    time = Time(clock.time)

    # Extract state variables at cell centers
    @inbounds begin
        # Ocean state
        uₒ = ℑxᶜᵃᵃ(i, j, 1, grid, ocean_state.u)
        vₒ = ℑyᵃᶜᵃ(i, j, 1, grid, ocean_state.v)
        Tₒ = ocean_state.T[i, j, 1]
        Tₒ = convert_to_kelvin(ocean_temperature_units, Tₒ)
        Sₒ = ocean_state.S[i, j, 1]
    end

    kᴺ = size(grid, 3) # index of the top ocean cell

    # Convert the native grid velocities to a zonal - meridional 
    # frame of reference (assuming the frame of reference is 
    # latitude - longitude here, we might want to change it)
    uₒ, vₒ = convert_to_intrinsic_reference_frame(i, j, kᴺ, grid, uₒ, vₒ)
        
    @inbounds begin
        # Atmos state, which is _assumed_ to exist at location = (c, c, nothing)
        # The third index "k" should not matter but we put the correct index to get
        # a surface node anyways.
        X = node(i, j, kᴺ + 1, grid, c, c, f)
        atmos_args = (atmos_grid, atmos_times, atmos_backend, atmos_time_indexing)

        uₐ = interp_atmos_time_series(atmos_state.u, X, time, atmos_args...)
        vₐ = interp_atmos_time_series(atmos_state.v, X, time, atmos_args...)

        Tₐ = interp_atmos_time_series(atmos_state.T, X, time, atmos_args...)
        pₐ = interp_atmos_time_series(atmos_state.p, X, time, atmos_args...)
        qₐ = interp_atmos_time_series(atmos_state.q, X, time, atmos_args...)
    end

    # Build thermodynamic and dynamic states in the atmosphere and surface.
    # Notation:
    #   ⋅ 𝒬 ≡ thermodynamic state vector
    #   ⋅ 𝒰 ≡ "dynamic" state vector (thermodynamics + reference height + velocity)
    ℂₐ = atmos_thermodynamics_parameters
    𝒬ₐ = thermodynamic_atmospheric_state = AtmosphericThermodynamics.PhaseEquil_pTq(ℂₐ, pₐ, Tₐ, qₐ)

    hₐ = atmosphere_reference_height # elevation of atmos variables relative to surface
    Uₐ = SVector(uₐ, vₐ)
    𝒰ₐ = dynamic_atmos_state = SurfaceFluxes.StateValues(hₐ, Uₐ, 𝒬ₐ)

    # Build surface state with saturated specific humidity
    surface_type = AtmosphericThermodynamics.Liquid()
    qₒ = seawater_saturation_specific_humidity(ℂₐ, Tₒ, Sₒ, 𝒬ₐ,
                                               similarity_theory.water_mole_fraction,
                                               similarity_theory.water_vapor_saturation,
                                               surface_type)
    
    # Thermodynamic and dynamic surface state
    𝒬₀ = thermodynamic_surface_state = AtmosphericThermodynamics.PhaseEquil_pTq(ℂₐ, pₐ, Tₒ, qₒ)

    h₀ = zero(grid) # surface height
    Uₒ = SVector(uₒ, vₒ)
    𝒰₀ = dynamic_ocean_state = SurfaceFluxes.StateValues(h₀, Uₒ, 𝒬₀)

    Qv = similarity_theory.fields.latent_heat
    Qc = similarity_theory.fields.sensible_heat
    Fv = similarity_theory.fields.water_vapor
    τx = similarity_theory.fields.x_momentum
    τy = similarity_theory.fields.y_momentum

    g = default_gravitational_acceleration
    ϰ = similarity_theory.von_karman_constant
    
    inactive = inactive_node(i, j, kᴺ, grid, c, c, c)
    maxiter  = ifelse(inactive, 1, similarity_theory.maxiter)

    turbulent_fluxes = compute_similarity_theory_fluxes(similarity_theory,
                                                        dynamic_ocean_state,
                                                        dynamic_atmos_state,
                                                        atmosphere_boundary_layer_height,
                                                        ℂₐ, g, ϰ, maxiter)

    # Convert back from a zonal - meridional flux to the frame of 
    # reference of the native ocean grid
    τˣ, τʸ = convert_to_extrinsic_reference_frame(i, j, kᴺ, grid, turbulent_fluxes.x_momentum, 
                                                                  turbulent_fluxes.y_momentum)

    @inbounds begin
        # +0: cooling, -0: heating
        Qv[i, j, 1] = ifelse(inactive, 0, turbulent_fluxes.latent_heat)
        Qc[i, j, 1] = ifelse(inactive, 0, turbulent_fluxes.sensible_heat)
        Fv[i, j, 1] = ifelse(inactive, 0, turbulent_fluxes.water_vapor)
        τx[i, j, 1] = ifelse(inactive, 0, τˣ)
        τy[i, j, 1] = ifelse(inactive, 0, τʸ)
    end
end

@kernel function _assemble_atmosphere_ocean_fluxes!(centered_velocity_fluxes,
                                                    net_tracer_fluxes,
                                                    grid,
                                                    clock,
                                                    ocean_temperature,
                                                    ocean_salinity,
                                                    ocean_temperature_units,
                                                    similarity_theory_fields,
                                                    downwelling_radiation,
                                                    prescribed_freshwater_flux,
                                                    atmos_grid,
                                                    atmos_times,
                                                    atmos_backend,
                                                    atmos_time_indexing,
                                                    runoff_args,
                                                    radiation_properties,
                                                    ocean_reference_density,
                                                    ocean_heat_capacity,
                                                    freshwater_density)

    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3)
    time = Time(clock.time)

    @inbounds begin
        Tₒ = ocean_temperature[i, j, 1]
        Tₒ = convert_to_kelvin(ocean_temperature_units, Tₒ)
        Sₒ = ocean_salinity[i, j, 1]

        X = node(i, j, kᴺ + 1, grid, c, c, f)
        atmos_args = (atmos_grid, atmos_times, atmos_backend, atmos_time_indexing)

        Qs = interp_atmos_time_series(downwelling_radiation.shortwave, X, time, atmos_args...)
        Qℓ = interp_atmos_time_series(downwelling_radiation.longwave,  X, time, atmos_args...)

        # Accumulate mass fluxes of freshwater due to rain, snow, rivers, icebergs, and whatever else.
        # Rememeber runoff fluxes could be `nothing` if rivers and icebergs are not included in the forcing
        Mp = interp_atmos_time_series(prescribed_freshwater_flux, X, time, atmos_args...)
        Mr = get_runoff_flux(X, time, runoff_args) 

        Qc = similarity_theory_fields.sensible_heat[i, j, 1] # sensible or "conductive" heat flux
        Qv = similarity_theory_fields.latent_heat[i, j, 1]   # latent heat flux
        Mv = similarity_theory_fields.water_vapor[i, j, 1]   # mass flux of water vapor
        τx = similarity_theory_fields.x_momentum[i, j, 1]    # zonal momentum flux
        τy = similarity_theory_fields.y_momentum[i, j, 1]    # meridional momentum flux
    end

    # Compute heat fluxes, bulk flux first
    Qd = net_downwelling_radiation(i, j, grid, time, radiation_properties, Qs, Qℓ)
    Qu = net_upwelling_radiation(i, j, grid, time, radiation_properties, Tₒ)

    ΣQ = Qd + Qu + Qc + Qv

    # Convert from a mass flux to a volume flux (aka velocity)
    # by dividing by the density of freshwater.
    # Also switch the sign, for some reason we are given freshwater flux as positive down.
    ρᶠ = freshwater_density
    ΣF = - (Mp + Mr) / ρᶠ

    # Add the contribution from the turbulent water vapor flux
    Fv = Mv / ρᶠ
    ΣF += Fv

    # Compute fluxes for u, v, T, S from momentum, heat, and freshwater fluxes
    Jᵘ = centered_velocity_fluxes.u
    Jᵛ = centered_velocity_fluxes.v
    Jᵀ = net_tracer_fluxes.T
    Jˢ = net_tracer_fluxes.S

    ρₒ = ocean_reference_density
    cₒ = ocean_heat_capacity

    atmos_ocean_Jᵘ = τx / ρₒ
    atmos_ocean_Jᵛ = τy / ρₒ
    atmos_ocean_Jᵀ = ΣQ / (ρₒ * cₒ)
    atmos_ocean_Jˢ = - Sₒ * ΣF

    # Mask fluxes over land for convenience
    inactive = inactive_node(i, j, kᴺ, grid, c, c, c)

    @inbounds begin
        Jᵘ[i, j, 1] = ifelse(inactive, 0, atmos_ocean_Jᵘ)
        Jᵛ[i, j, 1] = ifelse(inactive, 0, atmos_ocean_Jᵛ)
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

# Retrieve the details of runoff fluxes (rivers and icebergs, if present in the simulation).
# Note that these forcing fields are different in terms of frequency (daily instead of three-hourly)
# and gridsize (1/4 degree instead of 1/2 degree) when compared to the other prescribed fluxes
# So they need to be interpolated using their own grid / times / backend / time_indexing
@inline get_runoff_args(::Nothing) = nothing

@inline function get_runoff_args(runoff_flux)

    data    = map(ϕ -> ϕ.data, runoff_flux)
    grid    = runoff_flux.rivers.grid
    times   = runoff_flux.rivers.times
    backend = runoff_flux.rivers.backend
    time_indexing = runoff_flux.rivers.time_indexing

    return (data, grid, times, backend, time_indexing)
end

@inline get_runoff_flux(X, time, ::Nothing) = zero(eltype(X))

@inline function get_runoff_flux(X, time, runoff_args)
    
    @inbounds runoff_flux = runoff_args[1] # The data is located at position 1 of the tuple
    @inbounds other_args  = runoff_args[2:end] # Other args contain grid, times, backend and time_indexing
    
    Mr = interp_atmos_time_series(runoff_flux, X, time, other_args...)

    return Mr
end
