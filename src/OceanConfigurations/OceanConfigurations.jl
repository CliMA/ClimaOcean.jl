module OceanConfigurations

using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: IsopycnalSkewSymmetricDiffusivity,
                                       ConvectiveAdjustmentVerticalDiffusivity

using NumericalEarth.Oceans: ocean_simulation, default_ocean_closure
using NumericalEarth.Bathymetry: regrid_bathymetry, ORCAGrid
using NumericalEarth.DataWrangling.ORCA: ORCA1

export latitude_longitude_ocean,
       half_degree_tripolar_ocean,
       one_degree_tripolar_ocean,
       sixth_degree_tripolar_ocean,
       orca_ocean,
       simplified_ocean_closure

#####
##### Shared utilities
#####

@inline νhb(i, j, k, grid, ℓx, ℓy, ℓz, clock, fields, λ) =
    Oceananigans.Operators.Az(i, j, k, grid, ℓx, ℓy, ℓz)^2 / λ

# Background tracer diffusivity following Henyey et al. (1986), as implemented
# by Harrison and Hallberg (2008) and used in GFDL OM4 (Adcroft et al. 2019).
# Pivot value κ = 1.5e-5 m²/s at ±30° latitude, equatorial value κ = 2e-6 m²/s.
@inline henyey_diffusivity(x, y, z, t) = max(2e-6, 3e-5 * abs(sind(y)))

"""
    simplified_ocean_closure(FT=Oceananigans.defaults.FloatType)

A minimal closure suitable for testing on memory-limited GPUs (e.g. P100).
Uses `ConvectiveAdjustmentVerticalDiffusivity` with a background vertical
diffusivity and viscosity, avoiding the large parameter space of CATKE +
Gent-McWilliams + biharmonic closures.
"""
function simplified_ocean_closure(FT=Oceananigans.defaults.FloatType)
    return ConvectiveAdjustmentVerticalDiffusivity(FT;
               convective_κz = 1.0,
               background_κz = 1e-5,
               convective_νz = 1.0,
               background_νz = 1e-4)
end

# Standard vertical coordinate for all configurations.
# 60 levels, exponential spacing, 6000 m depth.
function vertical_coordinate(; Nz=60, depth=6000, zstar=false)
    return ExponentialDiscretization(Nz, -depth, 0; mutable=zstar)
end

#####
##### Configuration constructors
#####

include("latitude_longitude.jl")
include("half_degree_tripolar.jl")
include("one_degree_tripolar.jl")
include("sixth_degree_tripolar.jl")
include("orca.jl")

end # module
