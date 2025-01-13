using Oceananigans.Operators: intrinsic_vector
using Oceananigans.Grids: _node

adjust_fluxes_over_sea_ice!(velocity_fluxes, tracer_fluxes, ::Nothing, grid, kernel_parameters) = nothing

function adjust_fluxes_over_sea_ice!(centered_velocity_fluxes, net_tracer_fluxes, sea_ice,
                                     grid, kernel_parameters)
    
    ice_concentration = sea_ice_concentration(sea_ice)

    launch!(architecture(grid), grid, kernel_parameters, _adjust_fluxes_over_sea_ice!,
            centered_velocity_fluxes, net_tracer_fluxes, grid, ice_concentration)
end

@kernel function _adjust_fluxes_over_sea_ice!(centered_velocity_fluxes, net_tracer_fluxes, grid, ice_concentration)

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

#####
##### Surface flux computation
#####

function compute_atmosphere_ocean_fluxes!(coupled_model)
    ocean = coupled_model.ocean
    atmosphere = coupled_model.atmosphere
    sea_ice = coupled_model.sea_ice
    grid = ocean.model.grid
    arch = architecture(grid)
    clock = coupled_model.clock

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

    surface_phase = AtmosphericThermodynamics.Liquid()

    surface_atmosphere_fields = coupled_model.fluxes.surface_atmosphere_state

    # Simplify NamedTuple to reduce parameter space consumption.
    # See https://github.com/CliMA/ClimaOcean.jl/issues/116.
    surface_atmosphere_data = (u = surface_atmosphere_fields.u.data,
                               v = surface_atmosphere_fields.v.data,
                               T = surface_atmosphere_fields.T.data,
                               p = surface_atmosphere_fields.p.data,
                               q = surface_atmosphere_fields.q.data,
                               Qs = surface_atmosphere_fields.Qs.data,
                               Qℓ = surface_atmosphere_fields.Qℓ.data,
                               Mp = surface_atmosphere_fields.Mp.data)

    kernel_parameters = surface_computations_kernel_parameters(grid)

    launch!(arch, grid, kernel_parameters,
            _compute_atmosphere_surface_similarity_theory_fluxes!,
            turbulent_fluxes,
            similarity_theory,
            grid,
            clock,
            ocean_state,
            surface_phase, 
            coupled_model.fluxes.ocean_reference_density,
            coupled_model.fluxes.ocean_heat_capacity,
            coupled_model.fluxes.ocean_temperature_units,
            surface_atmosphere_data,
            radiation_properties.stefan_boltzmann_constant,
            radiation_properties.reflection.ocean,
            radiation_properties.emission.ocean,
            coupled_model.fluxes.turbulent.water_mole_fraction,
            coupled_model.fluxes.turbulent.water_vapor_saturation,
            atmosphere.reference_height, # height at which the state is known
            atmosphere.boundary_layer_height,
            atmosphere.thermodynamics_parameters)

    #####
    ##### Finally cobble together and properly interpolate fluxes
    ##### to be used by the ocean model.
    #####

    interpolated_prescribed_freshwater_flux = surface_atmosphere_data.Mp
    interpolated_downwelling_radiation = (shortwave = surface_atmosphere_data.Qs,
                                          longwave = surface_atmosphere_data.Qℓ)
    
    launch!(arch, grid, kernel_parameters,
            _assemble_atmosphere_ocean_fluxes!,
            centered_velocity_fluxes,
            net_tracer_fluxes,
            grid,
            clock,
            ocean_state.T,
            ocean_state.S,
            coupled_model.fluxes.ocean_temperature_units,
            turbulent_fluxes, 
            interpolated_downwelling_radiation,
            interpolated_prescribed_freshwater_flux,
            radiation_properties.stefan_boltzmann_constant,
            radiation_properties.reflection.ocean,
            radiation_properties.emission.ocean,
            coupled_model.fluxes.ocean_reference_density,
            coupled_model.fluxes.ocean_heat_capacity,
            coupled_model.fluxes.freshwater_density)
                
    adjust_fluxes_over_sea_ice!(centered_velocity_fluxes, net_tracer_fluxes,
                                sea_ice, grid, kernel_parameters)

    launch!(arch, grid, :xy, reconstruct_momentum_fluxes!,
            grid, staggered_velocity_fluxes, centered_velocity_fluxes)

    return nothing
end

