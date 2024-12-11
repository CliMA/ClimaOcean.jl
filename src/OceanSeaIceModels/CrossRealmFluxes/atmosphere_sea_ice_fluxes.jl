using ClimaSeaIce: SeaIceModel

sea_ice_thickness(sea_ice::Simulation{<:SeaIceModel}) = sea_ice.model.ice_thickness
sea_ice_thickness(::Nothing) = nothing

sea_ice_concentration(sea_ice::Simulation{<:SeaIceModel}) = sea_ice.model.ice_concentration
sea_ice_concentration(::Nothing) = nothing

function limit_fluxes_over_sea_ice!(grid, kernel_parameters, sea_ice,
                                    centered_velocity_fluxes,
                                    net_tracer_fluxes, args...)
    
    ice_concentration = sea_ice_concentration(sea_ice)
    launch!(architecture(grid), grid, kernel_parameters, _limit_fluxes_over_sea_ice!,
            centered_velocity_fluxes, net_tracer_fluxes, ice_concentration)
end

@kernel function _limit_fluxes_over_sea_ice!(centered_velocity_fluxes, net_tracer_fluxes, ice_concentration)

    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3)
    
    τx = centered_velocity_fluxes.u
    τy = centered_velocity_fluxes.v
    Jᵀ = net_tracer_fluxes.T
    Jˢ = net_tracer_fluxes.S

    @inbounds begin
        ℵ = ice_concentration[i, j, kᴺ]
        τx[i, j, kᴺ] = τx[i, j, kᴺ] * (1 - ℵ)
        τy[i, j, kᴺ] = τy[i, j, kᴺ] * (1 - ℵ)
        Jᵀ[i, j, kᴺ] = Jᵀ[i, j, kᴺ] * (1 - ℵ)
        Jˢ[i, j, kᴺ] = Jˢ[i, j, kᴺ] * (1 - ℵ)
    end
end

function compute_atmosphere_sea_ice_fluxes!(coupled_model)
    surface_atmosphere_state = coupled_model.fluxes.surface_atmosphere_state
    sea_ice = coupled_model.sea_ice

    # Basic model properties
    grid = sea_ice.model.grid
    arch = architecture(grid)
    clock = coupled_model.clock

    #####
    ##### First interpolate atmosphere time series
    ##### in time and to the ocean grid.
    #####
    
    # Assumption: The sea ice and the ocean are on the same horizontal 
    # grid, so we do not need to reinterpolate the atmospheric state
    surface_atmosphere_state = (u = surface_atmosphere_state.u.data,
                                v = surface_atmosphere_state.v.data,
                                T = surface_atmosphere_state.T.data,
                                p = surface_atmosphere_state.p.data,
                                q = surface_atmosphere_state.q.data,
                                Qs = surface_atmosphere_state.Qs.data,
                                Qℓ = surface_atmosphere_state.Qℓ.data,
                                Mp = surface_atmosphere_state.Mp.data)

    #####
    ##### Next compute turbulent fluxes.
    #####

     # Fluxes, and flux contributors
    centered_velocity_fluxes = (u = coupled_model.fluxes.total.sea_ice.momentum.uᶜᶜᶜ,
                                v = coupled_model.fluxes.total.sea_ice.momentum.vᶜᶜᶜ)

    staggered_velocity_fluxes = (u = coupled_model.fluxes.total.sea_ice.momentum.u,
                                 v = coupled_model.fluxes.total.sea_ice.momentum.v)

    net_tracer_fluxes    = coupled_model.fluxes.total.sea_ice.tracers
    similarity_theory    = coupled_model.fluxes.sea_ice.turbulent
    radiation_properties = coupled_model.fluxes.radiation

    sea_ice_state = (u = sea_ice.model.velocities.u,
                     v = sea_ice.model.velocities.v,
                     T = sea_ice.model.tracers.T,
                     S = sea_ice.model.tracers.S)

    launch!(arch, grid, kernel_parameters,
            _compute_atmosphere_surface_similarity_theory_fluxes!,
            similarity_theory,
            grid,
            clock,
            sea_ice_state,
            AtmosphericThermodynamics.Solid(), # surface phase
            coupled_model.fluxes.sea_ice_reference_density,
            coupled_model.fluxes.sea_ice_heat_capacity,
            coupled_model.fluxes.sea_ice_temperature_units,
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
            _assemble_atmosphere_sea_ice_fluxes!,
            centered_velocity_fluxes,
            net_tracer_fluxes,
            grid,
            clock,
            # sea_state.T,
            # sea_state.S,
            coupled_model.fluxes.ocean_temperature_units,
            similarity_theory.fields,
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

# Nothing yet...