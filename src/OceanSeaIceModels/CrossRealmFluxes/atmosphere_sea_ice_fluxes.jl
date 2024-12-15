limit_fluxes_over_sea_ice!(grid, kernel_parameters, ::Nothing, args...) = nothing

function limit_fluxes_over_sea_ice!(grid, kernel_parameters, sea_ice,
                                    centered_velocity_fluxes,
                                    net_tracer_fluxes, args...)
    
    ice_concentration = sea_ice_concentration(sea_ice)
    launch!(architecture(grid), grid, kernel_parameters, _limit_fluxes_over_sea_ice!,
            centered_velocity_fluxes, net_tracer_fluxes, grid, ice_concentration)
end

@kernel function _limit_fluxes_over_sea_ice!(centered_velocity_fluxes, net_tracer_fluxes, grid, ice_concentration)

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

# Nothing yet...
compute_atmosphere_sea_ice_fluxes!(coupled_model) = nothing
