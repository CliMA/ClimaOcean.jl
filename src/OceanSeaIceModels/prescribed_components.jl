using Oceananigans.TimeSteppers: tick!
using Oceananigans.OutputReaders: extract_field_time_series, update_field_time_series!
import Oceananigans.TimeSteppers: time_step!, update_state!

abstract type PrescribedComponent end

@inline function update_state!(component::PrescribedOcean)
    time = Time(component.clock.time)
    ftses = extract_field_time_series(component)

    for fts in ftses
        update_field_time_series!(fts, time)
    end
    return nothing
end

@inline function time_step!(component::PrescribedOcean, Δt)
    tick!(component.clock, Δt)

    update_state!(component)

    return nothing
end

# No need to compute anything here...
update_net_fluxes!(coupled_model, ::PrescribedAtmosphere) = nothing
