using ClimaSeaIce: SeaIceModel

sea_ice_thickness(sea_ice::Simulation{<:SeaIceModel}) = sea_ice.model.ice_thickness
sea_ice_thickness(::Nothing) = nothing

# Nothing yet...
function compute_atmosphere_sea_ice_fluxes!(coupled_model)
    ocean = coupled_model.ocean
    atmosphere = coupled_model.atmosphere
    atmosphere_grid = atmosphere.grid
    sea_ice = coupled_model.sea_ice.model

    # Basic model properties
    grid = ocean.model.grid
    arch = architecture(grid)
    clock = coupled_model.clock

    atmosphere_velocities = map(u -> u.data, atmosphere.velocities)
    atmosphere_pressure   = atmosphere.pressure.data

    atmosphere_state = merge(atmosphere_velocities, (; p=atmosphere_pressure))

    u = atmosphere.velocities.u # for example 
    atmosphere_times = u.times
    atmosphere_backend = u.backend
    atmosphere_time_indexing = u.time_indexing

    launch!(arch, grid, :xy, _compute_atmosphere_sea_ice_stress!,
            sea_ice.external_momentum_stresses,
            grid,
            clock,
            sea_ice.ice_dynamics.ocean_ice_drag_coefficient,
            atmosphere_state,
            atmosphere_grid,
            atmosphere_times,
            atmosphere_backend,
            atmosphere_time_indexing)   
    
    return nothing
end

@kernel function _compute_atmosphere_sea_ice_stress!(τ, grid, clock, Cᴰ,
                                                     atmos_state, 
                                                     atmos_grid,
                                                     atmos_times,
                                                     atmos_backend,
                                                     atmos_time_indexing)

    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3)
    time = Time(clock.time)
            
    @inbounds begin
        # Atmos state, which is _assumed_ to exist at location = (c, c, nothing)
        # The third index "k" should not matter but we put the correct index to get
        # a surface node anyways.
        Xᵢ⁻ = node(i,   j,   kᴺ+1, grid, f, c, f)
        Xⱼ⁻ = node(i,   j,   kᴺ+1, grid, c, f, f)
        Xᵢ⁺ = node(i+1, j,   kᴺ+1, grid, f, c, f)
        Xⱼ⁺ = node(i,   j+1, kᴺ+1, grid, c, f, f)

        atmos_args = (atmos_grid, atmos_times, atmos_backend, atmos_time_indexing)

        uₐ⁻ = interp_atmos_time_series(atmos_state.u, Xᵢ⁻, time, atmos_args...)
        vₐ⁻ = interp_atmos_time_series(atmos_state.v, Xⱼ⁻, time, atmos_args...)
        uₐ⁺ = interp_atmos_time_series(atmos_state.u, Xᵢ⁺, time, atmos_args...)
        vₐ⁺ = interp_atmos_time_series(atmos_state.v, Xⱼ⁺, time, atmos_args...)

        uₐ = (uₐ⁻ + uₐ⁺) / 2
        vₐ = (vₐ⁻ + vₐ⁺) / 2

        𝒰ₐ = sqrt(uₐ^2 + vₐ^2) 

        τu = 1.3 * Cᴰ * 𝒰ₐ * uₐ
        τv = 1.3 * Cᴰ * 𝒰ₐ * vₐ

        τˣ, τʸ = intrinsic_vector(i, j, kᴺ, grid, τu, τv)

        τ.u[i, j, 1] = τˣ
        τ.v[i, j, 1] = τʸ
    end
end