using ClimaOcean
import ClimaOcean.InitialConditions: interpolate!

using Oceananigans
using Oceananigans.Operators
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
import ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: convert_to_latlon_frame, convert_to_native_frame

@inline hack_cosd(φ) = cos(π * φ / 180)
@inline hack_sind(φ) = sin(π * φ / 180)

# Here we assume that the tripolar grid is locally orthogonal
@inline function convert_to_latlong_frame(i, j, grid::TRG, uₒ, vₒ)

    φᶜᶠᵃ₊ = φnode(i, j+1, 1, grid, Center(), Face(), Center())
    φᶜᶠᵃ₋ = φnode(i,   j, 1, grid, Center(), Face(), Center())
    Δyᶜᶜᵃ = Δyᶜᶜᶜ(i,   j, 1, grid)

    ũ = deg2rad(φᶜᶠᵃ₊ - φᶜᶠᵃ₋) / Δyᶜᶜᵃ

    φᶠᶜᵃ₊ = φnode(i+1, j, 1, grid, Face(), Center(), Center())
    φᶠᶜᵃ₋ = φnode(i,   j, 1, grid, Face(), Center(), Center())
    Δxᶜᶜᵃ = Δxᶜᶜᶜ(i,   j, 1, grid)

    ṽ = - deg2rad(φᶠᶜᵃ₊ - φᶠᶜᵃ₋) / Δxᶜᶜᵃ

    𝒰 = sqrt(ũ^2 + ṽ^2)

    d₁ = ũ / 𝒰
    d₂ = ṽ / 𝒰

    return uₒ * d₁ - vₒ * d₂, uₒ * d₂ + vₒ * d₁
end

@inline function convert_to_native_frame(i, j, grid::TRG, uₒ, vₒ) 

    φᶜᶠᵃ₊ = φnode(i, j+1, 1, grid, Center(), Face(), Center())
    φᶜᶠᵃ₋ = φnode(i,   j, 1, grid, Center(), Face(), Center())
    Δyᶜᶜᵃ = Δyᶜᶜᶜ(i,   j, 1, grid)

    ũ = deg2rad(φᶜᶠᵃ₊ - φᶜᶠᵃ₋) / Δyᶜᶜᵃ

    φᶠᶜᵃ₊ = φnode(i+1, j, 1, grid, Face(), Center(), Center())
    φᶠᶜᵃ₋ = φnode(i,   j, 1, grid, Face(), Center(), Center())
    Δxᶜᶜᵃ = Δxᶜᶜᶜ(i,   j, 1, grid)

    ṽ = - deg2rad(φᶠᶜᵃ₊ - φᶠᶜᵃ₋) / Δxᶜᶜᵃ

    𝒰 = sqrt(ũ^2 + ṽ^2)

    d₁ = ũ / 𝒰
    d₂ = ṽ / 𝒰

    return uₒ * d₁ + vₒ * d₂, uₒ * d₂ - vₒ * d₁
end
