# Common test setup file to make stand-alone tests easy
include("runtests_setup.jl")

using CUDA
using PythonCall
using CondaPkg
using Scratch

test_group = get(ENV, "TEST_GROUP", :all)
test_group = Symbol(test_group)

using ClimaOcean.DataWrangling: download_dataset

function delete_inpainted_files(dir)
    @info "Cleaning inpainted files..."
    for (root, _, files) in walkdir(dir)
        for file in files
            if endswith(file, "_inpainted.jld2")
                filepath = joinpath(root, file)
                rm(filepath; force=true)
                @info "    Deleted: $filepath"
            end
        end
    end
end


if test_group == :init || test_group == :all
    #####
    ##### Delete inpainted files
    #####

    delete_inpainted_files(@get_scratch!("."))

    #####
    ##### Download bathymetry data
    #####

    download_bathymetry()

    #####
    ##### Download JRA55 data
    #####

    atmosphere = JRA55PrescribedAtmosphere(backend=JRA55NetCDFBackend(2))

    #####
    ##### Download Dataset data
    #####

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

if test_group == :ecco2_monthly || test_group == :all
    include("test_ecco2_monthly.jl")
end

if test_group == :ecco2_daily || test_group == :all
    include("test_ecco2_daily.jl")
end

if test_group == :ecco4_en4 || test_group == :all
    include("test_ecco4_en4.jl")
end

# Tests that we can download JRA55 utilities
if test_group == :downloading || test_group == :all
    include("test_downloading.jl")
end

# Tests that we can download JRA55 utilities
if test_group == :copernicus_downloading || test_group == :all
    include("test_copernicus_downloading.jl")
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
