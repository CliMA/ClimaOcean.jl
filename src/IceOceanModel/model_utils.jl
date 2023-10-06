using Oceananigans
using Oceananigans.Utils: Time
using Oceananigans.Grids: architecture
using Oceananigans.Models: AbstractModel
import Oceananigans.Grids: launch!

launch!(model::AbstractModel, args...; kwargs...) = launch!(architecture(model.grid), model.grid, args...; kwargs...) 

@inline getflux(f::Number,                  i::Int, j::Int, grid::AbstractGrid, clock, fields)  = f
@inline getflux(f::Function,                i::Int, j::Int, grid::AbstractGrid, clock, fields)  = f(i, j, grid, clock, fields)
@inline getflux(f::AbstractArray{<:Any, 2}, i::Int, j::Int, grid::AbstractGrid, args...)        = @inbounds f[i, j]
@inline getflux(f::AbstractField,           i::Int, j::Int, grid::AbstractGrid, args...)        = @inbounds f[i, j, 1]
@inline getflux(f::FieldTimeSeries,         i::Int, j::Int, grid::AbstractGrid, clock, args...) = @inbounds f[i, j, Time(clock.time)]

# If we have ice, do not compute fluxes!
@inline function get_flux(ice_thickness, f, i::Int, j::Int, grid::AbstractGrid, args...)
    h = @inbounds ice_thickness[i, j, 1]
    return ifelse(h > 0, getflux(f, i, j, grid,args...), 0)
end