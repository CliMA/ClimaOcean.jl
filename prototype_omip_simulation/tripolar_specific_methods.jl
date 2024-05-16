using ClimaOcean
import ClimaOcean.InitialConditions: interpolate!

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
                          λnode, φnode

using OrthogonalSphericalShellGrids: TRG

@inline hack_cosd(φ) = cos(π * φ / 180)
@inline hack_sind(φ) = sin(π * φ / 180)

import ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: convert_to_latlong_frame, convert_to_native_frame

# Here we assume that the tripolar grid is locally orthogonal
@inline function convert_to_latlong_frame(i, j, grid::TRG, uₒ, vₒ)
    φ₁ = φnode(i,   j, 1, grid, Face(), Center(), Center())
    φ₂ = φnode(i+1, j, 1, grid, Face(), Center(), Center())
     
    θ = φ₂ - φ₁
    
    return uₒ * hack_cosd(θ) + vₒ * hack_sind(θ), uₒ * hack_sind(θ) + vₒ * hack_cosd(θ)
end

@inline function convert_to_native_frame(i, j, grid::TRG, uₒ, vₒ) 
    φ₁ = φnode(i, j,   1, grid, Face(), Center(), Center())
    φ₂ = φnode(i, j+1, 1, grid, Face(), Center(), Center())
     
    θ = φ₂ - φ₁
    
    return uₒ * hack_cosd(θ) + vₒ * hack_sind(θ), uₒ * hack_sind(θ) + vₒ * hack_cosd(θ)
end
