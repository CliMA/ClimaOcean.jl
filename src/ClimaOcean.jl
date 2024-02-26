module ClimaOcean

export regrid_bathymetry
export stretched_vertical_faces
export PowerLawStretching, LinearStretching
export jra55_field_time_series
export ecco2_field, ECCO2Metadata
export initialize!
export OceanSeaIceModel, FreezingLimitedOceanTemperature

using Oceananigans
using Oceananigans.Operators: ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ
using DataDeps
using CubicSplines

function __init__(; remove_existing_data=false)

    ## Data for the one_degree_global_simulation
    branch_url = "https://github.com/CliMA/OceananigansArtifacts.jl/raw/glw/near-global-data"
    dir = "lat_lon_bathymetry_and_fluxes"
    bathymetry_name                  = "bathymetry_lat_lon_360_150.jld2"
    initial_conditions_month_1       = "initial_conditions_month_01_360_150_48.jld2"
    initial_conditions_month_2       = "initial_conditions_month_02_360_150_48.jld2"
    surface_boundary_conditions_name = "surface_boundary_conditions_12_months_360_150.jld2"

    bathymetry_url = joinpath(branch_url, dir, bathymetry_name)
    initial_conditions_1_url = joinpath(branch_url, dir, initial_conditions_month_1)
    initial_conditions_2_url = joinpath(branch_url, dir, initial_conditions_month_2)
    surface_boundary_conditions_url = joinpath(branch_url, dir, surface_boundary_conditions_name)

    dep = DataDep("near_global_one_degree",
                  "Bathymetry, initial conditions, and surface boundary conditions for " *
                  "near-global one degree simulations",
                  [bathymetry_url, initial_conditions_1_url, initial_conditions_2_url, surface_boundary_conditions_url])

    DataDeps.register(dep)

    ## Data for the quarter_degree_global_simulation
    path = "https://github.com/CliMA/OceananigansArtifacts.jl/raw/ss/new_hydrostatic_data_after_cleared_bugs/quarter_degree_near_global_input_data/"

    datanames = ["bathymetry-1440x600",
                "temp-1440x600-latitude-75",
                "salt-1440x600-latitude-75",
                "tau_x-1440x600-latitude-75",
                "tau_y-1440x600-latitude-75",
                "initial_conditions"]

    dh_quarter = DataDep("near_global_quarter_degree",
        "Forcing data for global latitude longitude simulation",
        [path * data * ".jld2" for data in datanames]
    )

    DataDeps.register(dh_quarter)

    remove_existing_data && rm(datadep"near_global_one_degree",     recursive=true, force=true)
    remove_existing_data && rm(datadep"near_global_quarter_degree", recursive=true, force=true)
end

@inline ϕ²(i, j, k, grid, ϕ) = @inbounds ϕ[i, j, k]^2
@inline spᶠᶜᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, ϕ², Φ.v))
@inline spᶜᶠᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, ϕ², Φ.u))

@inline u_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds - μ * Φ.u[i, j, 1] * spᶠᶜᶜ(i, j, 1, grid, Φ)
@inline v_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds - μ * Φ.v[i, j, 1] * spᶜᶠᶜ(i, j, 1, grid, Φ)

@inline u_immersed_bottom_drag(i, j, k, grid, c, Φ, μ) = @inbounds - μ * Φ.u[i, j, k] * spᶠᶜᶜ(i, j, k, grid, Φ)
@inline v_immersed_bottom_drag(i, j, k, grid, c, Φ, μ) = @inbounds - μ * Φ.v[i, j, k] * spᶜᶠᶜ(i, j, k, grid, Φ)

include("OceanSeaIceModels/OceanSeaIceModels.jl")
include("VerticalGrids.jl")
include("InitialConditions/InitialConditions.jl")
include("DataWrangling/DataWrangling.jl")
include("Bathymetry.jl")
include("Diagnostics.jl")
include("NearGlobalSimulations/NearGlobalSimulations.jl")

using .VerticalGrids
using .Bathymetry
using .DataWrangling: JRA55
using .DataWrangling: ECCO2
using .InitialConditions
using .OceanSeaIceModels: OceanSeaIceModel

end # module

