using Oceananigans.Operators: intrinsic_vector
using Oceananigans.Grids: inactive_node
using ClimaOcean.OceanSeaIceModels.PrescribedAtmospheres: thermodynamics_parameters,
                                                          surface_layer_height,
                                                          boundary_layer_height

function compute_atmosphere_ocean_fluxes!(coupled_model)
    ocean = coupled_model.ocean
    atmosphere = coupled_model.atmosphere
    grid = ocean.model.grid
    arch = architecture(grid)
    clock = coupled_model.clock

    ocean_state = (u = ocean.model.velocities.u,
                   v = ocean.model.velocities.v,
                   T = ocean.model.tracers.T,
                   S = ocean.model.tracers.S)

    atmosphere_fields = coupled_model.interfaces.exchanger.exchange_atmosphere_state

    # Simplify NamedTuple to reduce parameter space consumption.
    # See https://github.com/CliMA/ClimaOcean.jl/issues/116.
    atmosphere_data = (u = atmosphere_fields.u.data,
                       v = atmosphere_fields.v.data,
                       T = atmosphere_fields.T.data,
                       p = atmosphere_fields.p.data,
                       q = atmosphere_fields.q.data,
                       Qs = atmosphere_fields.Qs.data,
                       Qℓ = atmosphere_fields.Qℓ.data,
                       Mp = atmosphere_fields.Mp.data,
                       h_bℓ = boundary_layer_height(atmosphere))

    flux_formulation = coupled_model.interfaces.atmosphere_ocean_interface.flux_formulation
    interface_fluxes = coupled_model.interfaces.atmosphere_ocean_interface.fluxes
    interface_temperature = coupled_model.interfaces.atmosphere_ocean_interface.temperature
    interface_properties = coupled_model.interfaces.atmosphere_ocean_interface.properties
    ocean_properties = coupled_model.interfaces.ocean_properties
    atmosphere_properties = (thermodynamics_parameters = thermodynamics_parameters(atmosphere),
                             surface_layer_height = surface_layer_height(atmosphere))

    kernel_parameters = interface_kernel_parameters(grid)

    launch!(arch, grid, kernel_parameters,
            _compute_atmosphere_ocean_interface_state!,
            interface_fluxes,
            interface_temperature,
            grid,
            clock,
            flux_formulation,
            ocean_state,
            atmosphere_data,
            interface_properties,
            atmosphere_properties,
            ocean_properties)

    return nothing
end

""" Compute turbulent fluxes between an atmosphere and a interface state using similarity theory """
@kernel function _compute_atmosphere_ocean_interface_state!(interface_fluxes,
                                                            interface_temperature,
                                                            grid,
                                                            clock,
                                                            turbulent_flux_formulation,
                                                            interior_state,
                                                            atmosphere_state,
                                                            interface_properties,
                                                            atmosphere_properties,
                                                            ocean_properties)

    i, j = @index(Global, NTuple)
    kᴺ   = size(grid, 3) # index of the top ocean cell
    time = Time(clock.time)

    @inbounds begin
        uₐ = atmosphere_state.u[i, j, 1]
        vₐ = atmosphere_state.v[i, j, 1]
        Tₐ = atmosphere_state.T[i, j, 1]
        pₐ = atmosphere_state.p[i, j, 1]
        qₐ = atmosphere_state.q[i, j, 1]
        Qs = atmosphere_state.Qs[i, j, 1]
        Qℓ = atmosphere_state.Qℓ[i, j, 1]

        # Extract state variables at cell centers
        # Ocean state
        uᵢ = ℑxᶜᵃᵃ(i, j, kᴺ, grid, interior_state.u)
        vᵢ = ℑyᵃᶜᵃ(i, j, kᴺ, grid, interior_state.v)
        Tᵢ = interior_state.T[i, j, kᴺ]
        Tᵢ = convert_to_kelvin(ocean_properties.temperature_units, Tᵢ)
        Sᵢ = interior_state.S[i, j, kᴺ]
    end

    # Build thermodynamic and dynamic states in the atmosphere and interface.
    # Notation:
    #   ⋅ 𝒬 ≡ thermodynamic state vector
    #   ⋅ 𝒰 ≡ "dynamic" state vector (thermodynamics + reference height + velocity)
    ℂₐ = atmosphere_properties.thermodynamics_parameters
    𝒬ₐ = thermodynamic_atmospheric_state = AtmosphericThermodynamics.PhaseEquil_pTq(ℂₐ, pₐ, Tₐ, qₐ)
    zₐ = atmosphere_properties.surface_layer_height # elevation of atmos variables relative to interface

    local_atmosphere_state = (z = zₐ,
                              u = uₐ,
                              v = vₐ,
                              𝒬 = 𝒬ₐ,
                              h_bℓ = atmosphere_state.h_bℓ)

    local_interior_state = (u=uᵢ, v=vᵢ, T=Tᵢ, S=Sᵢ)
    downwelling_radiation = (; Qs, Qℓ)

    # Estimate initial interface state
    FT = eltype(grid)
    u★ = convert(FT, 1e-4)

    # Estimate interface specific humidity using interior temperature
    q_formulation = interface_properties.specific_humidity_formulation
    qₛ = saturation_specific_humidity(q_formulation, ℂₐ, 𝒬ₐ.ρ, Tᵢ, Sᵢ)
    initial_interface_state = InterfaceState(u★, u★, u★, uᵢ, vᵢ, Tᵢ, Sᵢ, qₛ)

    # Don't use convergence criteria in an inactive cell
    is_water = !(inactive_node(i, j, kᴺ, grid, Center(), Center(), Center()))
    
    interface_state = compute_interface_state(turbulent_flux_formulation,
                                              is_water,
                                              initial_interface_state,
                                              local_atmosphere_state,
                                              local_interior_state,
                                              downwelling_radiation,
                                              interface_properties,
                                              atmosphere_properties,
                                              ocean_properties)

    # In the case of FixedIterations, make sure interface state is zero'd
    interface_state = ifelse(not_water, zero_interface_state(FT), interface_state)

    u★ = interface_state.u★
    θ★ = interface_state.θ★
    q★ = interface_state.q★

    Ψₛ = interface_state
    Ψₐ = local_atmosphere_state
    Δu, Δv = velocity_difference(interface_properties.velocity_formulation, Ψₐ, Ψₛ)
    ΔU = sqrt(Δu^2 + Δv^2)

    τx = ifelse(ΔU == 0, zero(grid), - u★^2 * Δu / ΔU)
    τy = ifelse(ΔU == 0, zero(grid), - u★^2 * Δv / ΔU)

    ρₐ = AtmosphericThermodynamics.air_density(ℂₐ, 𝒬ₐ)
    cₚ = AtmosphericThermodynamics.cp_m(ℂₐ, 𝒬ₐ) # moist heat capacity
    ℰv = AtmosphericThermodynamics.latent_heat_vapor(ℂₐ, 𝒬ₐ)

    # Store fluxes
    Qv  = interface_fluxes.latent_heat
    Qc  = interface_fluxes.sensible_heat
    Fv  = interface_fluxes.water_vapor
    ρτx = interface_fluxes.x_momentum
    ρτy = interface_fluxes.y_momentum
    Ts  = interface_temperature

    @inbounds begin
        # +0: cooling, -0: heating
        Qv[i, j, 1]  = - ρₐ * u★ * q★ * ℰv
        Qc[i, j, 1]  = - ρₐ * cₚ * u★ * θ★
        Fv[i, j, 1]  = - ρₐ * u★ * q★
        ρτx[i, j, 1] = + ρₐ * τx
        ρτy[i, j, 1] = + ρₐ * τy
        Ts[i, j, 1]  = convert_from_kelvin(ocean_properties.temperature_units, Ψₛ.T)
    end
end
