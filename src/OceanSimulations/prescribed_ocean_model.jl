using Oceananigans.OutputReaders: extract_field_time_series, update_field_time_series!

import Oceananigans.TimeSteppers: time_step!, update_state!, reset!, tick!
import Oceananigans.Models: timestepper, update_model_field_time_series!

import ClimaOcean: reference_density, heat_capacity
import Oceananigans.Architectures: on_architecture

#####
##### A prescribed ocean...
#####

struct PrescribedOcean{G, C, U, T, F, Arch} <: AbstractModel{Nothing, Arch}
    architecture :: Arch      
    grid :: G        
    clock :: Clock{C}
    velocities :: U
    tracers :: T
    timeseries :: F
end

"""
    PrescribedOcean(timeseries=NamedTuple(); grid, clock=Clock{Float64}(time = 0))

Create a prescribed ocean model to be used in combination with ClimaOcean's `OceanSeaIceModel` 
on a `grid` with a `clock`.

Arguments
=========
- `timeseries::NamedTuple`: A named tuple containing time series data for various fields. The named tuple can be empty or 
                            include any combination of the following fields: `u`, `v`, `T`, `S`.
                            All elements provided must be of type `FieldTimeSeries` and reside on the provided `grid`.
"""
function PrescribedOcean(timeseries=NamedTuple(); 
                         grid, 
                         clock=Clock{Float64}(time = 0)) 

    # Make sure all elements of the timeseries are on the same grid
    # If we decide to pass a timeseries
    if !isempty(timeseries)
        for (k, v) in timeseries
            isa(v, FieldTimeSeries) ||
                throw(ArgumentError("All variables in the `timeseries` argument must be `FieldTimeSeries`"))
            v.grid == grid ||
                throw(ArgumentError("All variables in the timeseries reside on the provided grid"))
        end
    end

    τx = Field{Face, Center, Nothing}(grid)
    τy = Field{Center, Face, Nothing}(grid)
    Jᵀ = Field{Center, Center, Nothing}(grid)
    Jˢ = Field{Center, Center, Nothing}(grid)

    u = XFaceField(grid,  boundary_conditions=FieldBoundaryConditions(grid, (Face,   Center, Center), top = FluxBoundaryCondition(τx)))
    v = YFaceField(grid,  boundary_conditions=FieldBoundaryConditions(grid, (Center, Face,   Center), top = FluxBoundaryCondition(τy)))
    T = CenterField(grid, boundary_conditions=FieldBoundaryConditions(grid, (Center, Center, Center), top = FluxBoundaryCondition(Jᵀ)))
    S = CenterField(grid, boundary_conditions=FieldBoundaryConditions(grid, (Center, Center, Center), top = FluxBoundaryCondition(Jˢ)))

    return PrescribedOcean(architecture(grid), grid, clock, (; u, v, w=ZeroField()), (; T, S), timeseries)
end

#####
##### Need to extend a couple of methods
#####

function time_step!(model::PrescribedOcean, Δt; callbacks=[], euler=true)
    tick!(model.clock, Δt)
    time = Time(model.clock.time)

    possible_fts = merge(model.velocities, model.tracers)
    time_series_tuple = extract_field_time_series(possible_fts)

    for fts in time_series_tuple
        update_field_time_series!(fts, time)
    end

    update_u_velocity  = haskey(model.timeseries, :u)
    update_v_velocity  = haskey(model.timeseries, :v)
    update_temperature = haskey(model.timeseries, :T)
    update_salinity    = haskey(model.timeseries, :S)

    update_u_velocity  && parent(model.velocities.u) .= parent(model.timeseries.u[time])
    update_v_velocity  && parent(model.velocities.v) .= parent(model.timeseries.v[time])
    update_temperature && parent(model.tracers.T)    .= parent(model.timeseries.T[time])
    update_salinity    && parent(model.tracers.S)    .= parent(model.timeseries.S[time])

    return nothing
end

update_state!(::PrescribedOcean) = nothing
timestepper(::PrescribedOcean) = nothing

reference_density(ocean::Simulation{<:PrescribedOcean}) = 1025.6
heat_capacity(ocean::Simulation{<:PrescribedOcean}) = 3995.6
