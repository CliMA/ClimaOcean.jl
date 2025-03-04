
using Oceananigans
using Oceananigans.BoundaryConditions
using Oceananigans.Grids: architecture

import Oceananigans.TimeSteppers: time_step!, update_state!, tick!
import Oceananigans.Models: timestepper, update_model_field_time_series!

import ClimaOcean.OceanSeaIceModels: reference_density, heat_capacity

#####
##### A prescribed ocean...
#####

struct PrescribedOcean{A, G, C, U, T, F} <: AbstractModel{Nothing}
    architecture :: A       
    grid :: G        
    clock :: Clock{C}
    velocities :: U
    tracers :: T
    timeseries :: F
end

"""
    PrescribedOcean(timeseries; grid, clock=Clock{Float64}(time = 0))

Create a `PrescribedOcean` model on a specified `grid` with a `clock` that evolves
according to the data passed in `timeseries`.

`timeseries` _must_ be a `NamedTuple` containing `FieldTimeseries` for `u`, `v`, `T`, and `S`.
"""
function PrescribedOcean(timeseries; 
                         grid, 
                         clock=Clock{Float64}(time = 0)) 

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

    parent(model.velocities.u) .= parent(model.timeseries.u[time])
    parent(model.velocities.v) .= parent(model.timeseries.v[time])
    parent(model.tracers.T)    .= parent(model.timeseries.T[time])
    parent(model.tracers.S)    .= parent(model.timeseries.S[time])

    return nothing
end

update_state!(::PrescribedOcean) = nothing
timestepper(::PrescribedOcean) = nothing

reference_density(ocean::Simulation{<:PrescribedOcean}) = 1025.6
heat_capacity(ocean::Simulation{<:PrescribedOcean}) = 3995.6
