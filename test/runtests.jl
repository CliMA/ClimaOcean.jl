# Common test setup file to make stand-alone tests easy
include("runtests_setup.jl")

using CUDA

test_group = get(ENV, "TEST_GROUP", :all)
test_group = Symbol(test_group)

using ClimaOcean.ECCO: download_dataset

if test_group == :init || test_group == :all
    using CUDA
    CUDA.set_runtime_version!(v"12.6"; local_toolkit = true)
    CUDA.precompile_runtime()

    ####
    #### Download bathymetry data
    ####
    
    download_bathymetry() 

    ####
    #### Download JRA55 data 
    ####
    
    atmosphere = JRA55PrescribedAtmosphere()

    ####
    #### Download ECCO data 
    ####

    download_dataset(temperature_metadata)
    download_dataset(salinity_metadata)
end

# Tests JRA55 utilities, plus some DataWrangling utilities
if test_group == :jra55 || test_group == :all
    include("test_jra55.jl")
end

if test_group == :ecco || test_group == :all
    include("test_ecco.jl")
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

if test_group == :simulations || test_group == :all
    CUDA.set_runtime_version!(v"12.2", local_toolkit = true) # Seems to help in finding the correct CUDA version
    include("test_simulations.jl")
end

if test_group == :distributed || test_group == :all
    include("test_distributed_utils.jl")
end

