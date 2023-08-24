using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities:
    tracer_mixing_lengthᶜᶜᶠ,
    stability_functionᶜᶜᶠ,
    stable_length_scaleᶜᶜᶠ,
    convective_length_scaleᶜᶜᶠ,
    buoyancy_flux,
    shear_production,
    dissipation

function tracer_convective_length_scale_operation(model, closure=model.closure)
    velocities = model.velocities
    tracers = model.tracers
    buoyancy = model.buoyancy
    clock = model.clock
    Qᵇ = model.diffusivity_fields.Qᵇ
    # h = model.diffusivity_fields.h
    top_tracer_bcs = NamedTuple(c => tracers[c].boundary_conditions.top for c in propertynames(tracers))

    Cᶜ = closure.mixing_length.Cᶜc
    Cᵉ = closure.mixing_length.Cᵉc
    Cˢᵖ = closure.mixing_length.Cˢᵖ

    convective_length_scale_args = (closure, Cᶜ, Cᵉ, Cˢᵖ,
                                    velocities,
                                    tracers,
                                    buoyancy,
                                    Qᵇ)
                                    #Qᵇ, h)
                                    #clock, top_tracer_bcs)

    ℓᶜconv = KernelFunctionOperation{Center, Center, Face}(convective_length_scaleᶜᶜᶠ, grid, convective_length_scale_args...)

    return ℓᶜconv
end

function tracer_stable_length_scale_operation(model, closure=model.closure)
    velocities = model.velocities
    tracers = model.tracers
    buoyancy = model.buoyancy
    clock = model.clock

    Cˡᵒ = closure.mixing_length.Cˡᵒc
    Cʰⁱ = closure.mixing_length.Cʰⁱc

    σ = KernelFunctionOperation{Center, Center, Face}(stability_functionᶜᶜᶠ, grid, closure, Cˡᵒ, Cʰⁱ,
                                                      velocities, tracers, buoyancy)

    ℓ_stable = KernelFunctionOperation{Center, Center, Face}(stable_length_scaleᶜᶜᶠ, grid, closure, tracers.e,
                                                             velocities, tracers, buoyancy)

    ℓᶜshear = σ * ℓ_stable

    return ℓᶜshear
end

function tracer_mixing_length_operation(model, closure=model.closure)
    velocities = model.velocities
    tracers = model.tracers
    buoyancy = model.buoyancy
    clock = model.clock
    Qᵇ = model.diffusivity_fields.Qᵇ
    #h = model.diffusivity_fields.h
    top_tracer_bcs = NamedTuple(c => tracers[c].boundary_conditions.top for c in propertynames(tracers))

    ℓᶜ = KernelFunctionOperation{Center, Center, Face}(tracer_mixing_lengthᶜᶜᶠ, grid,
                                                       closure, velocities, tracers, buoyancy, Qᵇ)

    return ℓᶜ
end

function buoyancy_flux_operation(model, closure=model.closure)
    velocities = model.velocities
    tracers = model.tracers
    buoyancy = model.buoyancy
    diffusivities = model.diffusivity_fields
    wb = KernelFunctionOperation{Center, Center, Center}(buoyancy_flux, grid, closure,
                                                         velocities, tracers, buoyancy, diffusivities)
    return wb
end

function dissipation_operation(model, closure=model.closure)
    velocities = model.velocities
    tracers = model.tracers
    buoyancy = model.buoyancy
    diffusivities = model.diffusivity_fields
    ϵ = KernelFunctionOperation{Center, Center, Center}(dissipation, grid, closure,
                                                        velocities, tracers, buoyancy, diffusivities)
    return ϵ
end

function shear_production_operation(model, closure=model.closure)
    velocities = model.velocities
    tracers = model.tracers
    buoyancy = model.buoyancy
    diffusivities = model.diffusivity_fields
    P = KernelFunctionOperation{Center, Center, Center}(shear_production, grid, closure,
                                                        velocities, tracers, buoyancy, diffusivities)
    return P
end

