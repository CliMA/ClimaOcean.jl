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

struct PrescribedOcean{G, C, U, T, F, FT} 
    grid :: G        
    clock :: Clock{C}
    velocities :: U
    tracers :: T
    timeseries :: F
    reference_density :: FT
    heat_capacity :: FT
end

"""
    PrescribedOcean(timeseries=NamedTuple(); grid, clock=Clock{Float64}(time = 0))

Create a prescribed ocean model to be used in combination with ClimaOcean's `OceanSeaIceModel` 
on a `grid` with a `clock`.

Arguments
=========
- `timeseries`: A `NamedTuple` containing time series data for various fields. The named tuple can be empty 
                (meaning that the ocean does not evolve in time) or include any combination of the 
                following fields: `u`, `v`, `T`, `S`. All elements provided must be of type `FieldTimeSeries` 
                and reside on the provided `grid`.
"""
function PrescribedOcean(grid; 
                         velocities = default_ocean_velocities(grid), 
                         tracers = default_ocean_tracers(grid),
                         clock = Clock{Float64}(time = 0),
                         reference_density = 1029,
                         heat_capacity = 3998) 

    reference_density = convert(eltype(grid), reference_density)
    heat_capacity = convert(eltype(grid), heat_capacity)

    return PrescribedOcean(grid, 
                           clock, 
                           velocities, 
                           tracers,
                           reference_density,
                           heat_capacity)
end

default_ocean_velocities(grid) = (u = XFaceField(grid), v = YFaceField(grid))
default_ocean_tracers(grid) = (T = CenterField(grid), S = CenterField(grid))

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

    return nothing
end

get_ocean_state()

update_state!(::PrescribedOcean) = nothing
timestepper(::PrescribedOcean) = nothing

reference_density(ocean::PrescribedOcean) = ocean.reference_density
heat_capacity(ocean::PrescribedOcean) = ocean.heat_capacity
