module ClimaOcean

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

struct CubicSplineFunction{d, FT, S} <: Function
    coordinates :: Vector{FT}
    nodes :: Vector{FT}
    spline :: S
    
    @doc """
        CubicSplineFunction{d}(coordinates, nodes, FT=Float64)

    Return a function-like object that interpolates `nodes` between the `coordinates`
    in dimension `d` using cubic splines. `d` can be either `x`, `y`, or `z`.
    """
    function CubicSplineFunction{d}(coordinates, nodes, FT=Float64) where d
        # Hack to enforce no-gradient boundary conditions,
        # since CubicSplines doesn't support natively
        ΔL = coordinates[2] - coordinates[1]
        pushfirst!(coordinates, coordinates[1] - ΔL)

        ΔR = coordinates[end] - coordinates[end-1]
        push!(coordinates, coordinates[end] + ΔR)
    
        pushfirst!(nodes, nodes[1])
        push!(nodes, nodes[end])
    
        coordinates = Vector{FT}(coordinates)
        nodes = Vector{FT}(nodes)

        # Now we can build the spline
        spline = CubicSpline(coordinates, nodes)
        S = typeof(spline)

        d == :x || d == :y || d == :z || error("Dimension 'd' must be :x or :y or :z")
    
        return new{d, FT, S}(coordinates, nodes, spline)
    end
end

(csf::CubicSplineFunction{:x})(x, y=nothing, z=nothing) = csf.spline[x]
(csf::CubicSplineFunction{:y})(x, y, z=nothing)         = csf.spline[y]
(csf::CubicSplineFunction{:y})(y)                       = csf.spline[y]
(csf::CubicSplineFunction{:z})(x, y, z)                 = csf.spline[z]
(csf::CubicSplineFunction{:z})(z)                       = csf.spline[z]

@inline ϕ²(i, j, k, grid, ϕ) = @inbounds ϕ[i, j, k]^2
@inline spᶠᶜᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, ϕ², Φ.v))
@inline spᶜᶠᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, ϕ², Φ.u))

@inline u_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds - μ * Φ.u[i, j, 1] * spᶠᶜᶜ(i, j, 1, grid, Φ)
@inline v_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds - μ * Φ.v[i, j, 1] * spᶜᶠᶜ(i, j, 1, grid, Φ)

@inline u_immersed_bottom_drag(i, j, k, grid, c, Φ, μ) = @inbounds - μ * Φ.u[i, j, k] * spᶠᶜᶜ(i, j, k, grid, Φ)
@inline v_immersed_bottom_drag(i, j, k, grid, c, Φ, μ) = @inbounds - μ * Φ.v[i, j, k] * spᶜᶠᶜ(i, j, k, grid, Φ)

include("VerticalGrids.jl")
include("DataWrangling.jl")
include("Bathymetry.jl")
include("InitialConditions.jl")
include("Diagnostics.jl")
include("NearGlobalSimulations/NearGlobalSimulations.jl")
include("IdealizedSimulations/IdealizedSimulations.jl")

end # module
