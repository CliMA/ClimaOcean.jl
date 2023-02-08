module ClimaOcean

using Oceananigans
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities:
    CATKEVerticalDiffusivity,
    MixingLength,
    TurbulentKineticEnergyEquation

using DataDeps

function __init__(; remove_existing_data=false)

    # Data for one_degree_global_simulation
    base_url            = "https://github.com/glwagner/ClimaOceanData/raw/main/near_global_simulation_data"
    bathymetry_url                  = joinpath(base_url, "near_global_bathymetry_360_150.jld2")
    initial_conditions_url          = joinpath(base_url, "near_global_initial_conditions_360_150_48.jld2")
    surface_boundary_conditions_url = joinpath(base_url, "near_global_boundary_conditions_360_150.jld2")

    one_degree_urls = [bathymetry_url, initial_conditions_url, surface_boundary_conditions_url]
    one_degree_desc = "Bathymetry, initial conditions, and surface boundary conditions " *
                      "for near-global one degree simulations"

    dep = DataDep("near_global_one_degree", one_degree_desc, one_degree_urls)
    DataDeps.register(dep)

    # Data for quarter_degree_global_simulation
    bathymetry_url          = joinpath(base_url, "near_global_bathymetry_1440_600.jld2")
    surface_salinity_url    = joinpath(base_url, "near_global_surface_salinity_1440_600.jld2")
    surface_temperature_url = joinpath(base_url, "near_global_surface_temperature_1440_600.jld2")
    east_momentum_flux_url  = joinpath(base_url, "near_global_east_momentum_flux_1440_600.jld2")
    north_momentum_flux_url = joinpath(base_url, "near_global_north_momentum_flux_1440_600.jld2")

    quarter_degree_urls = [bathymetry_url, surface_salinity_url, surface_temperature_url,
                           east_momentum_flux_url, north_momentum_flux_url]
    quarter_degree_desc =  "Bathymetry, initial conditions, and surface boundary conditions " * 
                           "for near-global quarter degree degree simulations"

    dep = DataDep("near_global_quarter_degree", quarter_degree_desc, quarter_degree_urls)
    DataDeps.register(dep)

    remove_existing_data && rm(datadep"near_global_one_degree",     recursive=true, force=true)
    remove_existing_data && rm(datadep"near_global_quarter_degree", recursive=true, force=true)
end

turbulent_kinetic_energy_equation = TurbulentKineticEnergyEquation(
    C⁻D   = 1.0,
    C⁺D   = 1.0,
    CᶜD   = 0.0,
    CᵉD   = 0.0,
    Cᵂu★  = 1.0,
    CᵂwΔ  = 1.0,
)

mixing_length = MixingLength(
    Cᵇ   = Inf,
    Cˢ   = Inf,
    Cᶜc  = 0.0,
    Cᶜe  = 0.0,
    Cᵉc  = 0.0,
    Cᵉe  = 0.0,
    Cˢᶜ  = 0.0,
    C⁻u  = 1.0,
    C⁺u  = 1.0,
    C⁻c  = 1.0,
    C⁺c  = 1.0,
    C⁻e  = 1.0,
    C⁺e  = 1.0,
    CRiʷ = 1.0,
    CRiᶜ = 0.0,
)

neutral_catke = CATKEVerticalDiffusivity(; mixing_length, turbulent_kinetic_energy_equation)

include("VerticalGrids.jl")
include("NearGlobalSimulations.jl")
include("Diagnostics.jl")
#include("DataWrangling.jl")

end # module

