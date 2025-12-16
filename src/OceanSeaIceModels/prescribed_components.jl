using Oceananigans.TimeSteppers: tick!
using Oceananigans.OutputReaders: extract_field_time_series, update_field_time_series!
import Oceananigans.TimeSteppers: time_step!, update_state!

abstract type AbstractPrescribedComponent end

@inline function update_state!(component::AbstractPrescribedComponent)
    time = Time(component.clock.time)
    ftses = extract_field_time_series(component)

    for fts in ftses
        update_field_time_series!(fts, time)
    end
    return nothing
end

@inline function time_step!(component::AbstractPrescribedComponent, Δt)
    tick!(component.clock, Δt)

    update_state!(component)

    return nothing
end

# No need to compute anything here...
net_fluxes(::AbstractPrescribedComponent) = nothing
update_net_fluxes!(coupled_model, ::AbstractPrescribedComponent) = nothing
