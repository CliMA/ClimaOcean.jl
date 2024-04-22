using ClimaOcean
import ClimaOcean.InitialConditions: three_dimensional_interpolate!

using Oceananigans
using Oceananigans.BoundaryConditions
using Oceananigans.Fields: OneField
using Oceananigans.Grids: peripheral_node
using Oceananigans.Utils: launch!
using Oceananigans.Fields: instantiated_location, interior, CenterField
using Oceananigans.Architectures: architecture, device, GPU, child_architecture

# Implementation of 3-dimensional regridding
# TODO: move all the following to Oceananigans! 

using Oceananigans.Fields: regrid!, interpolate!
using Oceananigans.Grids: cpu_face_constructor_x, 
                          cpu_face_constructor_y, 
                          cpu_face_constructor_z,
                          topology,
                          λnode

using OrthogonalSphericalShellGrids: TRG, WRG

const TField = Field{<:Any, <:Any, <:Any, <:Any, <:TRG}
const WField = Field{<:Any, <:Any, <:Any, <:Any, <:WRG}

function three_dimensional_interpolate!(a::Union{<:TField, <:WField}, b)

    interpolate!(a, b)

    return a
end

import ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: convert_to_latlong, convert_to_native_grid

# Here we assume that the tripolar grid is locally orthogonal
@inline function convert_to_latlong(i, j, grid::WRG, uₒ, vₒ)
    λ₁ = λnode(i, j,   1, grid, Center(), Face(), Center())
    λ₂ = λnode(i, j+1, 1, grid, Center(), Face(), Center())
     
    θ = λ₂ - λ₁
    
    return uₒ * cosd(θ) + vₒ * sind(θ), uₒ * sind(θ) + vₒ * cosd(θ)
end

@inline function convert_to_native_grid(i, j, grid::WRG, uₒ, vₒ) 
    λ₁ = λnode(i, j,   1, grid, Center(), Face(), Center())
    λ₂ = λnode(i, j+1, 1, grid, Center(), Face(), Center())
     
    θ = λ₂ - λ₁
    
    return uₒ * cosd(θ) + vₒ * sind(θ), uₒ * sind(θ) + vₒ * cosd(θ)
end

# Here we assume that the tripolar grid is locally orthogonal
@inline function convert_to_latlong(i, j, grid::TRG, uₒ, vₒ)
    λ₁ = λnode(i, j,   1, grid, Center(), Face(), Center())
    λ₂ = λnode(i, j+1, 1, grid, Center(), Face(), Center())
     
    θ = λ₂ - λ₁
    
    return uₒ * cosd(θ) + vₒ * sind(θ), uₒ * sind(θ) + vₒ * cosd(θ)
end

@inline function convert_to_native_grid(i, j, grid::TRG, uₒ, vₒ) 
    λ₁ = λnode(i, j,   1, grid, Center(), Face(), Center())
    λ₂ = λnode(i, j+1, 1, grid, Center(), Face(), Center())
     
    θ = λ₂ - λ₁
    
    return uₒ * cosd(θ) + vₒ * sind(θ), uₒ * sind(θ) + vₒ * cosd(θ)
end
