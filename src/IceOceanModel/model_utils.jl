using Oceananigans.Grids: architecture
using Oceananigans.Models: AbstractModel
import Oceananigans.Grids: launch!

launch!(model::AbstractModel, args...; kwargs...) = launch!(architecture(model.grid), model.grid, args...; kwargs...) 


@inline getflux(f::Function, i, j, grid, clock, fields) = f(i, j, grid, clock, fields)
  