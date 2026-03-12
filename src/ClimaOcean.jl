module ClimaOcean

# Use the README as the module docs
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end ClimaOcean

using Reexport
using Oceananigans
using Oceananigans.Operators: ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ

@reexport using NumericalEarth
@reexport using NumericalEarth.DataWrangling
@reexport using NumericalEarth.EarthSystemModels
@reexport using NumericalEarth.EarthSystemModels.InterfaceComputations

#####
##### Source code
#####

include("InitialConditions/InitialConditions.jl")
# include("Bathymetry/Bathymetry.jl")
include("Diagnostics/Diagnostics.jl")

@reexport using NumericalEarth.DataWrangling: ETOPO, ECCO, GLORYS, EN4, JRA55
@reexport using NumericalEarth.Bathymetry
@reexport using NumericalEarth.EarthSystemModels
@reexport using NumericalEarth.Atmospheres
@reexport using NumericalEarth.Oceans
@reexport using NumericalEarth.SeaIces
using .InitialConditions

end # module
