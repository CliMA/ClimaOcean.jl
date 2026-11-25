module ClimaOceanMakieExt

using ClimaOcean
using ClimaOcean: OMIPOutput, output_path, convert_to_teos10!
using ClimaOcean.Observations: Observations, BoyerMontegutMonthly, HadISSTMonthly, NCEPMonthly

using Makie
using GeoMakie: GeoMakie, GeoAxis

using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: λnode, φnode, znodes, λnodes, φnodes
using Oceananigans.ImmersedBoundaries: inactive_node
using Oceananigans.Operators: Azᶜᶜᶜ, Δxᶜᶠᶜ, Vᶜᶜᶜ
using Oceananigans.OutputReaders: OnDisk

using NumericalEarth
using NumericalEarth.Bathymetry: atlantic_ocean_basin
using NumericalEarth.DataWrangling: all_dates, download_with_retries
using NumericalEarth.EarthSystemModels.InterfaceComputations: vector_component_boundary_conditions

using ConservativeRegridding: ConservativeRegridding
using DelimitedFiles: readdlm
using JLD2: jldopen
using NCDatasets: NCDataset
using Dates: Dates, DateTime, Millisecond, month, year
using Statistics: mean, std

include("style.jl")
include("output_reading.jl")
include("regridding.jl")
include("observations.jl")
include("surface.jl")
include("sea_ice.jl")
include("global_series.jl")
include("sections.jl")
include("overturning.jl")
include("straits.jl")
include("figures.jl")

end # module
