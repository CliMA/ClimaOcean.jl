using Oceananigans.Operators: intrinsic_vector
using Oceananigans.Grids: inactive_node
using Oceananigans.Fields: ZeroField

function compute_atmosphere_sea_ice_fluxes!(coupled_model)
    ocean = coupled_model.ocean
    atmosphere = coupled_model.atmosphere
    sea_ice = coupled_model.sea_ice
    grid = sea_ice.model.grid
    arch = architecture(grid)
    clock = coupled_model.clock

    interior_state = (u = ZeroField(),
                      v = ZeroField(),
                      h = sea_ice.model.ice_thickness,
                      Tₒ = ocean.model.tracers.T,
                      Sₒ = ocean.model.tracers.S)

    atmosphere_fields = coupled_model.fluxes.near_surface_atmosphere_state

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
                       zb = atmosphere.boundary_layer_height)

    flux_formulation = coupled_model.fluxes.atmosphere_sea_ice_interface.flux_formulation
    interface_fluxes = coupled_model.fluxes.atmosphere_sea_ice_interface.fluxes
    interface_temperature = coupled_model.fluxes.atmosphere_sea_ice_interface.temperature
    interface_properties = coupled_model.fluxes.atmosphere_sea_ice_interface.properties
    sea_ice_properties = coupled_model.fluxes.sea_ice_properties

    atmosphere_properties = (thermodynamics_parameters = atmosphere.thermodynamics_parameters,
                             reference_height = atmosphere.reference_height)

    kernel_parameters = surface_computations_kernel_parameters(grid)

    launch!(arch, grid, kernel_parameters,
            _compute_atmosphere_sea_ice_interface_state!,
            interface_fluxes,
            interface_temperature,
            grid,
            clock,
            flux_formulation,
            interior_state,
            atmosphere_data,
            interface_properties,
            atmosphere_properties,
            sea_ice_properties)

    return nothing
end

""" Compute turbulent fluxes between an atmosphere and a interface state using similarity theory """
@kernel function _compute_atmosphere_sea_ice_interface_state!(interface_fluxes,
                                                              interface_temperature,
                                                              grid,
                                                              clock,
                                                              turbulent_flux_formulation,
                                                              interior_state,
                                                              atmosphere_state,
                                                              interface_properties,
                                                              atmosphere_properties,
                                                              sea_ice_properties)

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
        # Ocean properties below sea ice
        Tᵢ = interior_state.Tₒ[i, j, kᴺ]
        Tᵢ = convert_to_kelvin(sea_ice_properties.temperature_units, Tᵢ)
        Sᵢ = interior_state.Sₒ[i, j, kᴺ]

        # Sea ice properties
        uᵢ = ℑxᶜᵃᵃ(i, j, 1, grid, interior_state.u)
        vᵢ = ℑyᵃᶜᵃ(i, j, 1, grid, interior_state.v)
        hᵢ = interior_state.h[i, j, 1]
        Tₛ = interface_temperature[i, j, 1]
    end

    # Build thermodynamic and dynamic states in the atmosphere and interface.
    # Notation:
    #   ⋅ 𝒬 ≡ thermodynamic state vector
    #   ⋅ 𝒰 ≡ "dynamic" state vector (thermodynamics + reference height + velocity)
    ℂₐ = atmosphere_properties.thermodynamics_parameters
    𝒬ₐ = thermodynamic_atmospheric_state = AtmosphericThermodynamics.PhaseEquil_pTq(ℂₐ, pₐ, Tₐ, qₐ)
    hₐ = atmosphere_reference_height # elevation of atmos variables relative to interface
    Uₐ = SVector(uₐ, vₐ)
    local_atmosphere_state = SurfaceFluxes.StateValues(hₐ, Uₐ, 𝒬ₐ)
    downwelling_radiation = (; Qs, Qℓ)
    local_interior_state = (u=uᵢ, v=vᵢ, T=Tᵢ, S=Sᵢ, h=hᵢ)
    atmosphere_properties = ℂₐ

    # Estimate initial interface state
    FT = eltype(grid)
    u★ = convert(FT, 1e-4)

    # Estimate interface specific humidity using interior temperature
    q_formulation = interface_properties.specific_humidity_formulation
    qₛ = saturation_specific_humidity(q_formulation, ℂₐ, 𝒬ₐ.ρ, Tₛ, Sᵢ) 
    initial_interface_state = InterfaceState(u★, u★, u★, uᵢ, vᵢ, Tₛ, Sᵢ, qₛ)

    if inactive_node(i, j, kᴺ, grid, Center(), Center(), Center())
        interface_state = zero_interface_state(FT)
    else
        interface_state = compute_interface_state(turbulent_flux_formulation,
                                                  initial_interface_state,
                                                  local_atmosphere_state,
                                                  local_interior_state,
                                                  downwelling_radiation,
                                                  interface_properties,
                                                  atmosphere_properties,
                                                  sea_ice_properties)
    end

    u★ = interface_state.u★
    θ★ = interface_state.θ★
    q★ = interface_state.q★

    #=
    Pr = similarity_theory.turbulent_prandtl_number
    θ★ = θ★ / Pr
    q★ = q★ / Pr
    =#

    Ψₛ = interface_state
    Ψₐ = local_atmosphere_state
    Δu, Δv = velocity_difference(turbulent_flux_formulation.bulk_velocity, Ψₐ, Ψₛ)
    ΔU = sqrt(Δu^2 + Δv^2)
    τx = - u★^2 * Δu / ΔU
    τy = - u★^2 * Δv / ΔU

    𝒬ₐ = local_atmosphere_state.ts
    ℂₐ = atmosphere_properties
    ρₐ = AtmosphericThermodynamics.air_density(ℂₐ, 𝒬ₐ)
    cₚ = AtmosphericThermodynamics.cp_m(ℂₐ, 𝒬ₐ) # moist heat capacity
    ℰv = AtmosphericThermodynamics.latent_heat_vapor(ℂₐ, 𝒬ₐ)

    # Store fluxes
    Qv = interface_fluxes.latent_heat
    Qc = interface_fluxes.sensible_heat
    Fv = interface_fluxes.water_vapor
    ρτx = interface_fluxes.x_momentum
    ρτy = interface_fluxes.y_momentum
    Ts = interface_temperature

    @inbounds begin
        # +0: cooling, -0: heating
        Qv[i, j, 1]  = - ρₐ * u★ * q★ * ℰv
        Qc[i, j, 1]  = - ρₐ * cₚ * u★ * θ★
        Fv[i, j, 1]  = - ρₐ * u★ * q★
        ρτx[i, j, 1] = + ρₐ * τx
        ρτy[i, j, 1] = + ρₐ * τy
        Ts[i, j, 1]  = interface_state.T
    end
end

