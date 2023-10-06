using Oceananigans
using Oceananigans.Utils: Time
using Oceananigans.Grids: architecture
using Oceananigans.Models: AbstractModel
import Oceananigans.Grids: launch!

launch!(model::AbstractModel, args...; kwargs...) = launch!(architecture(model.grid), model.grid, args...; kwargs...) 

@inline getflux(::Nothing, f::Function,                i, j, grid, clock, fields)  = f(i, j, grid, clock, fields)
@inline getflux(::Nothing, f::AbstractArray{<:Any, 2}, i, j, grid, args...)        = @inbounds f[i, j]
@inline getflux(::Nothing, f::AbstractField,           i, j, grid, args...)        = @inbounds f[i, j, 1]
@inline getflux(::Nothing, f::FieldTimeSeries,         i, j, grid, clock, args...) = @inbounds f[i, j, Time(clock.time)]

# If we have ice, do not compute fluxes!
@inline function get_flux(ice_thickness, f, i, j, args...)
    h = @inbounds ice_thickness[i, j, 1]
    return ifelse(h > 0, getflux(nothing, f, i, j, args...), 0)
end