using Oceananigans.Operators: intrinsic_vector
using Oceananigans.Grids: inactive_node

function compute_atmosphere_ocean_fluxes!(coupled_model)
    exchanger = coupled_model.interfaces.exchanger
    grid = exchanger.grid
    arch = architecture(grid)
    clock = coupled_model.clock
    ocean_state = exchanger.ocean.state
    atmosphere_fields = exchanger.atmosphere.state

    # Simplify NamedTuple to reduce parameter space consumption.
    # See https://github.com/CliMA/ClimaOcean.jl/issues/116.
    atmosphere_data = merge(atmosphere_fields, 
                            (; h_bℓ = boundary_layer_height(coupled_model.atmosphere)))

    flux_formulation = coupled_model.interfaces.atmosphere_ocean_interface.flux_formulation
    interface_fluxes = coupled_model.interfaces.atmosphere_ocean_interface.fluxes
    interface_temperature = coupled_model.interfaces.atmosphere_ocean_interface.temperature
    interface_properties = coupled_model.interfaces.atmosphere_ocean_interface.properties
    ocean_properties = coupled_model.interfaces.ocean_properties
    atmosphere_properties = (thermodynamics_parameters = thermodynamics_parameters(coupled_model.atmosphere),
                             surface_layer_height = surface_layer_height(coupled_model.atmosphere),
                             gravitational_acceleration = coupled_model.interfaces.properties.gravitational_acceleration)

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
    #   ⋅ 𝒰 ≡ "dynamic" state vector (thermodynamics + reference height + velocity)
    ℂₐ = atmosphere_properties.thermodynamics_parameters
    zₐ = atmosphere_properties.surface_layer_height # elevation of atmos variables relative to interface

    local_atmosphere_state = (z = zₐ,
                              u = uₐ,
                              v = vₐ,
                              T = Tₐ,
                              p = pₐ,
                              q = qₐ,
                              h_bℓ = atmosphere_state.h_bℓ)

    local_interior_state = (u=uᵢ, v=vᵢ, T=Tᵢ, S=Sᵢ)
    downwelling_radiation = (; Qs, Qℓ)

    # Estimate initial interface state
    FT = typeof(Tᵢ)
    u★ = convert(FT, 1e-4)

    # Estimate interface specific humidity using interior temperature
    q_formulation = interface_properties.specific_humidity_formulation
    qₛ = surface_specific_humidity(q_formulation, ℂₐ, Tₐ, pₐ, qₐ, Tᵢ, Sᵢ)
    initial_interface_state = InterfaceState(u★, u★, u★, uᵢ, vᵢ, Tᵢ, Sᵢ, qₛ)

    # Don't use convergence criteria in an inactive cell
    stop_criteria = turbulent_flux_formulation.solver_stop_criteria
    needs_to_converge = stop_criteria isa ConvergenceStopCriteria
    not_water = inactive_node(i, j, kᴺ, grid, Center(), Center(), Center())

    # Compute local radiative properties and rebuild the interface properties
    α = stateindex(interface_properties.radiation.α, i, j, kᴺ, grid, time, (Center, Center, Center), Qs)
    ϵ = stateindex(interface_properties.radiation.ϵ, i, j, kᴺ, grid, time, (Center, Center, Center))
    σ = interface_properties.radiation.σ

    interface_properties = InterfaceProperties((; α, ϵ, σ),
                                               interface_properties.specific_humidity_formulation,
                                               interface_properties.temperature_formulation,
                                               interface_properties.velocity_formulation)

    if needs_to_converge && not_water
        interface_state = zero_interface_state(FT)
    else
        interface_state = compute_interface_state(turbulent_flux_formulation,
                                                  initial_interface_state,
                                                  local_atmosphere_state,
                                                  local_interior_state,
                                                  downwelling_radiation,
                                                  interface_properties,
                                                  atmosphere_properties,
                                                  ocean_properties)
    end

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

    ρₐ = AtmosphericThermodynamics.air_density(ℂₐ, Tₐ, pₐ, qₐ)
    cₚ = AtmosphericThermodynamics.cp_m(ℂₐ, qₐ) # moist heat capacity
    ℒv = AtmosphericThermodynamics.latent_heat_vapor(ℂₐ, Tₐ)
    

    # Store fluxes
    Qv  = interface_fluxes.latent_heat
    Qc  = interface_fluxes.sensible_heat
    Fv  = interface_fluxes.water_vapor
    ρτx = interface_fluxes.x_momentum
    ρτy = interface_fluxes.y_momentum
    Ts  = interface_temperature

    @inbounds begin
        # +0: cooling, -0: heating
        Qv[i, j, 1]  = - ρₐ * ℒv * u★ * q★ 
        Qc[i, j, 1]  = - ρₐ * cₚ * u★ * θ★
        Fv[i, j, 1]  = - ρₐ * u★ * q★
        ρτx[i, j, 1] = + ρₐ * τx
        ρτy[i, j, 1] = + ρₐ * τy
        Ts[i, j, 1]  = convert_from_kelvin(ocean_properties.temperature_units, Ψₛ.T)

        interface_fluxes.friction_velocity[i, j, 1] = u★
        interface_fluxes.temperature_scale[i, j, 1] = θ★
        interface_fluxes.water_vapor_scale[i, j, 1] = q★
    end
end
