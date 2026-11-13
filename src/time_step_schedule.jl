"""
    TimeStepSchedule(initial_Δt, changes...)

A time step that starts at `initial_Δt` and changes at the times of the `time => Δt` pairs `changes`, in
increasing order of time. As a callback, it sets the time step of the simulation from the model time.

```jldoctest
using ClimaOcean, Oceananigans.Units

TimeStepSchedule(4minutes, 60days => 30minutes, 365days => 45minutes)

# output
TimeStepSchedule: 4 minutes, 30 minutes from 60 days, 45 minutes from 365 days
```
"""
struct TimeStepSchedule{T}
    initial_Δt :: T
    change_times :: Vector{T}
    Δts :: Vector{T}
end

function TimeStepSchedule(initial_Δt, changes::Pair...)
    T = promote_type(typeof(float(initial_Δt)), map(c -> typeof(float(last(c))), changes)...)
    return TimeStepSchedule(convert(T, initial_Δt), T[first(c) for c in changes], T[last(c) for c in changes])
end

function time_step(schedule::TimeStepSchedule, time)
    n = searchsortedlast(schedule.change_times, time)
    return n == 0 ? schedule.initial_Δt : schedule.Δts[n]
end

(schedule::TimeStepSchedule)(simulation) = simulation.Δt = time_step(schedule, time(simulation))

Oceananigans.prognostic_state(::TimeStepSchedule) = nothing

function Base.show(io::IO, schedule::TimeStepSchedule)
    changes = [string(prettytime(Δt), " from ", prettytime(t)) for (t, Δt) in zip(schedule.change_times, schedule.Δts)]
    print(io, "TimeStepSchedule: ", join([prettytime(schedule.initial_Δt); changes], ", "))
end
