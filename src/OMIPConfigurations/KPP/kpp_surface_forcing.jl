# Surface-forcing primitives: u★, non-solar buoyancy flux Bo, two-band SW
# penetration, and the depth-integrated buoyancy forcing Bf(d). KPP sign
# convention: positive Bo = stabilizing (negation of Oceananigans' top flux).

# Carrier for surface BCs passed kernel-side; built CPU-side per compute call.
struct KPPTopBoundaryConditions{V, T}
    velocities :: V
    tracers    :: T
end

Adapt.adapt_structure(to, b::KPPTopBoundaryConditions) =
    KPPTopBoundaryConditions(adapt(to, b.velocities), adapt(to, b.tracers))

#####
##### Friction velocity
#####

@inline function friction_velocity(i, j, grid, clock, fields, top_velocity_bcs, params)
    τx = getbc(top_velocity_bcs.u, i, j, grid, clock, fields)
    τy = getbc(top_velocity_bcs.v, i, j, grid, clock, fields)
    return max(sqrt(sqrt(τx^2 + τy^2)), params.minimum_friction_velocity)
end

#####
##### Non-solar surface buoyancy flux Bo (KPP sign: stabilizing positive)
#####

@inline non_solar_buoyancy(i, j, grid, clock, fields, buoyancy, top_tracer_bcs) =
    - top_buoyancy_flux(i, j, grid, buoyancy.formulation, top_tracer_bcs, clock, fields)

#####
##### Two-band SW penetration: fraction remaining at d, integrated buoyancy gain.
#####

@inline shortwave_fraction(d, ::Nothing) = zero(d)

@inline function shortwave_fraction(d, radiation)
    FT = typeof(d)
    ϵ₁ = radiation.first_color_fraction
    κ₁ = radiation.first_absorption_coefficient
    κ₂ = radiation.second_absorption_coefficient
    return ϵ₁ * exp(- κ₁ * d) + (one(FT) - ϵ₁) * exp(- κ₂ * d)
end

@inline solar_buoyancy_above(i, j, d, ::Nothing, α, g) = zero(d)

@inline function solar_buoyancy_above(i, j, d, radiation, α, g)
    FT = typeof(d)
    J₀ = @inbounds radiation.surface_flux[i, j, 1]
    return - g * α * J₀ * (one(FT) - shortwave_fraction(d, radiation))
end

@inline buoyancy_forcing_above(i, j, d, Bo, radiation, α, g) =
    Bo + solar_buoyancy_above(i, j, d, radiation, α, g)

#####
##### Surface-cell EOS coefficients
#####

@inline function αᶜᶜᶜ(i, j, grid, buoyancy, tracers)
    FT  = eltype(grid)
    eos = buoyancy.formulation.equation_of_state
    α   = thermal_expansionᶜᶜᶜ(i, j, grid.Nz, grid, eos, tracers.T, tracers.S)
    return ifelse(inactive_node(i, j, grid.Nz, grid, Center(), Center(), Center()), zero(FT), α)
end

@inline function βᶜᶜᶜ(i, j, grid, buoyancy, tracers)
    FT  = eltype(grid)
    eos = buoyancy.formulation.equation_of_state
    β   = haline_contractionᶜᶜᶜ(i, j, grid.Nz, grid, eos, tracers.T, tracers.S)
    return ifelse(inactive_node(i, j, grid.Nz, grid, Center(), Center(), Center()), zero(FT), β)
end
