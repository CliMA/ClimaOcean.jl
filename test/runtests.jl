# Common test setup file to make stand-alone tests easy
include("runtests_setup.jl")

using CUDA

test_group = get(ENV, "TEST_GROUP", :all)
test_group = Symbol(test_group)

using ClimaOcean.DataWrangling: download_dataset

if test_group == :init || test_group == :all
    using CUDA
    CUDA.set_runtime_version!(v"12.6"; local_toolkit = true)
    CUDA.precompile_runtime()

    ###
    ### Download bathymetry data
    ###

    download_bathymetry()

    ####
    #### Download JRA55 data
    ####

    atmosphere = JRA55PrescribedAtmosphere(backend=JRA55NetCDFBackend(2))

    ####
    #### Download Dataset data
    ####

    # Download few datasets for tests
    for dataset in test_datasets
        time_resolution = dataset isa ECCO2Daily ? Day(1) : Month(1)
        end_date = start_date + 2 * time_resolution
        dates = start_date:time_resolution:end_date

        temperature_metadata = Metadata(:temperature; dataset, dates)
        salinity_metadata    = Metadata(:salinity; dataset, dates)

        download_dataset(temperature_metadata)
        download_dataset(salinity_metadata)
    end
end

# Tests JRA55 utilities, plus some DataWrangling utilities
if test_group == :JRA55 || test_group == :all
    include("test_jra55.jl")
end

if test_group == :ecco2 || test_group == :all
    include("test_ecco2.jl")
end

if test_group == :ecco4_en4 || test_group == :all
    include("test_ecco4_en4.jl")
end

# Tests that we can download JRA55 utilities
if test_group == :downloading || test_group == :all
    include("test_downloading.jl")
end

if test_group == :fluxes || test_group == :all
    include("test_surface_fluxes.jl")
end

if test_group == :bathymetry || test_group == :all
    include("test_bathymetry.jl")
end

if test_group == :ocean_sea_ice_model || test_group == :all
    include("test_ocean_sea_ice_model.jl")
    include("test_diagnostics.jl")
end

if test_group == :distributed || test_group == :all
    include("test_distributed_utils.jl")
end

if test_group == :reactant || test_group == :all
    include("test_reactant.jl")
end
