#####
##### Surface flux computation
#####

function compute_atmosphere_sea_ice_fluxes!(coupled_model)
    ocean = coupled_model.ocean
    atmosphere = coupled_model.atmosphere
    sea_ice = coupled_model.sea_ice
    grid = ocean.model.grid
    arch = architecture(grid)
    clock = coupled_model.clock

    heat_flux            = coupled_model.fluxes.total.sea_ice.top_heat
    similarity_theory    = coupled_model.fluxes.turbulent.coefficients.sea_ice
    radiation_properties = coupled_model.fluxes.radiation
    turbulent_fluxes     = coupled_model.fluxes.turbulent.fields.sea_ice

    sea_ice_state = (u = ocean.model.velocities.u,
                     v = ocean.model.velocities.v,
                     T = ocean.model.tracers.T,
                     S = ocean.model.tracers.S)

    surface_phase = AtmosphericThermodynamics.Ice()

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
            _compute_atmosphere_sea_ice_turbulent_fluxes!,
            turbulent_fluxes,
            similarity_theory,
            grid,
            clock,
            sea_ice_state,
            surface_phase, 
            coupled_model.fluxes.ocean_reference_density,
            coupled_model.fluxes.ocean_heat_capacity,
            coupled_model.fluxes.ocean_temperature_units,
            surface_atmosphere_data,
            radiation_properties.stefan_boltzmann_constant,
            radiation_properties.reflection.sea_ice,
            radiation_properties.emission.sea_ice,
            coupled_model.fluxes.turbulent.water_mole_fraction,
            coupled_model.fluxes.turbulent.water_vapor_saturation,
            atmosphere.reference_height, # height at which the state is known
            atmosphere.boundary_layer_height,
            atmosphere.thermodynamics_parameters)

    #####
    ##### Finally cobble together and properly interpolate fluxes
    ##### to be used by the ocean model.
    #####

    interpolated_downwelling_radiation = (shortwave = surface_atmosphere_data.Qs,
                                          longwave = surface_atmosphere_data.Qℓ)
    
    launch!(arch, grid, kernel_parameters,
            _assemble_atmosphere_sea_ice_fluxes!,
            heat_flux,
            grid,
            clock,
            sea_ice_state.T,
            coupled_model.fluxes.ocean_temperature_units,
            turbulent_fluxes, 
            interpolated_downwelling_radiation,
            radiation_properties.stefan_boltzmann_constant,
            radiation_properties.reflection.sea_ice,
            radiation_properties.emission.sea_ice)
                
    return nothing
end

using Oceananigans.Grids: inactive_node

""" Compute turbulent fluxes between an atmosphere and a surface state using similarity theory """
@kernel function _compute_atmosphere_sea_ice_turbulent_fluxes!(fluxes_fields,
                                                               coefficients,
                                                               grid,
                                                               clock,
                                                               surface_state,
                                                               surface_phase,
                                                               surface_density,
                                                               surface_heat_capacity,
                                                               surface_temperature_units,
                                                               surface_atmos_state,
                                                               stefan_boltzmann_constant,
                                                               albedo,
                                                               emissivity,
                                                               mole_fraction,
                                                               vapor_saturation,
                                                               atmosphere_reference_height,
                                                               atmosphere_boundary_layer_height,
                                                               atmos_thermodynamics_parameters)

    i, j = @index(Global, NTuple)
    kᴺ   = size(grid, 3) # index of the top ocean cell
    time = Time(clock.time)

    @inbounds begin
        uₐ = surface_atmos_state.u[i, j, 1]
        vₐ = surface_atmos_state.v[i, j, 1]
        Tₐ = surface_atmos_state.T[i, j, 1]
        pₐ = surface_atmos_state.p[i, j, 1]
        qₐ = surface_atmos_state.q[i, j, 1]
        Qs = surface_atmos_state.Qs[i, j, 1]
        Qℓ = surface_atmos_state.Qℓ[i, j, 1]

        # Extract state variables at cell centers
        # Ocean state
        uₛ = ℑxᶜᵃᵃ(i, j, kᴺ, grid, surface_state.u)
        vₛ = ℑyᵃᶜᵃ(i, j, kᴺ, grid, surface_state.v)
        Tₛ = surface_state.T[i, j, kᴺ]
        Tₛ = convert_to_kelvin(surface_temperature_units, Tₛ)
        Sₛ = surface_state.S[i, j, kᴺ]
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
    qₛ = seawater_saturation_specific_humidity(ℂₐ, Tₛ, Sₛ, 𝒬ₐ,
                                               mole_fraction,
                                               vapor_saturation,
                                               surface_phase)
    
    # Thermodynamic and dynamic surface state:
    Uₛ = SVector(uₛ, vₛ)
     
    𝒬₀ = thermodynamic_surface_state = AtmosphericThermodynamics.PhaseEquil_pTq(ℂₐ, pₐ, Tₛ, qₛ)
    h₀ = zero(grid) # surface height
    𝒰₀ = dynamic_surface_state = SurfaceFluxes.StateValues(h₀, Uₛ, 𝒬₀)

    # Some parameters
    g = default_gravitational_acceleration
    
    surface_salinity = Sₛ
    σ = stefan_boltzmann_constant
    α = stateindex(albedo, i, j, 1, grid, time)
    ϵ = stateindex(emissivity, i, j, 1, grid, time)
    prescribed_heat_fluxes = net_downwelling_radiation(i, j, grid, time, α, ϵ, Qs, Qℓ) 
    inactive_cell = inactive_node(i, j, kᴺ, grid, Center(), Center(), Center())

    turbulent_fluxes, surface_temperature = compute_turbulent_fluxes(coefficients,
                                                                     dynamic_surface_state, 
                                                                     dynamic_atmos_state, 
                                                                     prescribed_heat_fluxes,
                                                                     σ, α, ϵ,
                                                                     surface_phase,
                                                                     surface_salinity,
                                                                     surface_density,
                                                                     surface_heat_capacity,
                                                                     mole_fraction,
                                                                     vapor_saturation,
                                                                     atmosphere_boundary_layer_height,
                                                                     ℂₐ, g, inactive_cell)

    # Store fluxes
    Qv = fluxes_fields.latent_heat
    Qc = fluxes_fields.sensible_heat
    Fv = fluxes_fields.water_vapor
    ρτx = fluxes_fields.x_momentum
    ρτy = fluxes_fields.y_momentum
    Ts  = fluxes_fields.T_surface

    @inbounds begin
        # +0: cooling, -0: heating
        Qv[i, j, 1]  = ifelse(inactive_cell, 0, turbulent_fluxes.latent_heat)
        Qc[i, j, 1]  = ifelse(inactive_cell, 0, turbulent_fluxes.sensible_heat)
        Fv[i, j, 1]  = ifelse(inactive_cell, 0, turbulent_fluxes.water_vapor)
        ρτx[i, j, 1] = ifelse(inactive_cell, 0, turbulent_fluxes.x_momentum)
        ρτy[i, j, 1] = ifelse(inactive_cell, 0, turbulent_fluxes.y_momentum)
        Ts[i, j, 1]  = surface_temperature
    end
end

