module OceanConfigurations

using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: IsopycnalSkewSymmetricDiffusivity

using NumericalEarth.Oceans: ocean_simulation, default_ocean_closure
using NumericalEarth.Bathymetry: regrid_bathymetry, ORCAGrid
using NumericalEarth.DataWrangling.ORCA: ORCA1

export latitude_longitude_ocean,
       half_degree_tripolar_ocean,
       one_degree_tripolar_ocean,
       sixth_degree_tripolar_ocean,
       orca_ocean

#####
##### Shared utilities
#####

@inline νhb(i, j, k, grid, ℓx, ℓy, ℓz, clock, fields, λ) =
    Oceananigans.Operators.Az(i, j, k, grid, ℓx, ℓy, ℓz)^2 / λ

# Background tracer diffusivity following Henyey et al. (1986), as implemented
# by Harrison and Hallberg (2008) and used in GFDL OM4 (Adcroft et al. 2019).
# Pivot value κ = 1.5e-5 m²/s at ±30° latitude, equatorial value κ = 2e-6 m²/s.
@inline henyey_diffusivity(x, y, z, t) = max(2e-6, 3e-5 * abs(sind(y)))

# Standard vertical coordinate for all configurations.
# 60 levels, exponential spacing, 6000 m depth.
function vertical_coordinate(; zstar=false)
    return ExponentialDiscretization(60, -6000, 0; mutable=zstar)
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
