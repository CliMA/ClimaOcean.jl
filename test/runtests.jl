# Common test setup file to make stand-alone tests easy
using ClimaOcean
using Oceananigans
using CUDA
using Test

using Oceananigans.Architectures: architecture, on_architecture
using Oceananigans.OutputReaders: interpolate!
using Dates

using CUDA: @allowscalar

gpu_test = parse(Bool, get(ENV, "GPU_TEST", "false"))
test_architectures = gpu_test ? [GPU()] : [CPU()]
start_date = DateTime(1993, 1, 1)

using CUDA
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

    ETOPOmetadata = Metadatum(:bottom_height, dataset=ClimaOcean.ETOPO.ETOPO2022())
    ClimaOcean.download_dataset(ETOPOmetadata)


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

        ClimaOcean.download_dataset(temperature_metadata)
        ClimaOcean.download_dataset(salinity_metadata)

        if dataset isa Union{ECCO2DarwinMonthly, ECCO4DarwinMonthly}
            PO₄_metadata = Metadata(:phosphate; dataset, dates)
            ClimaOcean.download_dataset(PO₄_metadata)
        end
    end
end

if test_group == :bathymetry || test_group == :all
    include("test_bathymetry.jl")
end

if test_group == :ocean_sea_ice_model || test_group == :all
    include("test_ocean_sea_ice_model.jl")
    include("test_diagnostics.jl")
end
