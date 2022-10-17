module ClimaOcean

using Oceananigans
using DataDeps

function __init__()
    branch_url = "https://github.com/CliMA/OceananigansArtifacts.jl/raw/glw/near-global-data"
    dir = "lat_lon_bathymetry_and_fluxes"
    bathymetry_name = "bathymetry_lat_lon_360x150.jld2"
    initial_conditions_name = "initial_conditions_360x150x48.jld2"
    surface_boundary_conditions_name = "surface_boundary_conditions_360x150.jld2"

    bathymetry_url = joinpath(branch_url, dir, bathymetry_name)
    initial_conditions_url = joinpath(branch_url, dir, initial_conditions_name)
    surface_boundary_conditions_url = joinpath(branch_url, dir, surface_boundary_conditions_name)

    dep = DataDep("near_global_one_degree",
                  "Bathymetry, initial conditions, and surface boundary conditions for " *
                  "near-global one degree simulations",
                  [bathymetry_url, initial_conditions_url, surface_boundary_conditions_url])

    DataDeps.register(dep)
end

include("VerticalGrids.jl")
include("NearGlobalSimulations.jl")
include("Diagnostics.jl")

end # module

