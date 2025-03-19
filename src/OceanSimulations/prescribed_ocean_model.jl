using Oceananigans.Models: AbstractModel
using Oceananigans.Fields: ZeroField
using Oceananigans.OutputReaders: extract_field_time_series, update_field_time_series!

import Oceananigans.TimeSteppers: time_step!, update_state!, reset!, tick!
import Oceananigans.Models: timestepper, update_model_field_time_series!

import ClimaOcean: reference_density, heat_capacity
import Oceananigans.Architectures: on_architecture

#####
##### A prescribed ocean...
#####

struct PrescribedOceanModel{G, C, U, T, F, Arch} <: AbstractModel{Nothing, Arch}
    architecture :: Arch      
    grid :: G        
    clock :: Clock{C}
    velocities :: U
    tracers :: T
    timeseries :: F
end

"""
    PrescribedOceanModel(timeseries=NamedTuple(); grid, clock=Clock{Float64}(time = 0))

Create a prescribed ocean model to be used in combination with ClimaOcean's `OceanSeaIceModel` 
on a `grid` with a `clock`.

Arguments
=========
- `timeseries`: A `NamedTuple` containing time series data for various fields. The named tuple can be empty 
                (meaning that the ocean does not evolve in time) or include any combination of the 
                following fields: `u`, `v`, `T`, `S`. All elements provided must be of type `FieldTimeSeries` 
                and reside on the provided `grid`.
"""
function PrescribedOceanModel(timeseries=NamedTuple(); 
                         grid, 
                         clock=Clock{Float64}(time = 0)) 

    # Make sure all elements of the timeseries are on the same grid
    # If we decide to pass a timeseries
    if !isempty(timeseries)
        for k in keys(timeseries)
            f = timeseries[k]
            isa(f, FieldTimeSeries) ||
                throw(ArgumentError("All variables in the `timeseries` argument must be `FieldTimeSeries`"))
            f.grid == grid ||
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

    return PrescribedOceanModel(architecture(grid), grid, clock, (; u, v, w=ZeroField()), (; T, S), timeseries)
end

#####
##### Need to extend a couple of methods
#####

function time_step!(model::PrescribedOceanModel, Δt; callbacks=[], euler=true)
    tick!(model.clock, Δt)
    time = Time(model.clock.time)

    possible_fts = merge(model.velocities, model.tracers)
    time_series_tuple = extract_field_time_series(possible_fts)

    for fts in time_series_tuple
        update_field_time_series!(fts, time)
    end

    # Time stepping the model!

    if haskey(model.timeseries, :u)  
        arent(model.velocities.u) .= parent(model.timeseries.u[time])
    end

    if haskey(model.timeseries, :v)  
        parent(model.velocities.v) .= parent(model.timeseries.v[time])
    end
     
    if haskey(model.timeseries, :T)  
        parent(model.tracers.T) .= parent(model.timeseries.T[time])
    end

    if haskey(model.timeseries, :S)  
        parent(model.tracers.S) .= parent(model.timeseries.S[time])
    end

    return nothing
end

update_state!(::PrescribedOceanModel) = nothing
timestepper(::PrescribedOceanModel) = nothing

reference_density(ocean::Simulation{<:PrescribedOceanModel}) = 1025.6
heat_capacity(ocean::Simulation{<:PrescribedOceanModel}) = 3995.6
