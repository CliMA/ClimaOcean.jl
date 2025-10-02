include("runtests_setup.jl")

using ClimaOcean
using ClimaOcean.ECCO
using ClimaOcean.DataWrangling: NearestNeighborInpainting, metadata_path, native_times, download_dataset

using Dates
using Oceananigans.Grids: topology
using Oceananigans.OutputReaders: time_indices
using Oceananigans.TimeSteppers: update_state!
using Oceananigans.Units

using CUDA: @allowscalar

# Inpaint only the first two cells inside the missing mask
inpainting = NearestNeighborInpainting(2)

test_ecco_datasets = tuple((ds for ds in test_datasets if occursin(r"^ECCO2.*Monthly",string(typeof(ds)),))...)

start_date = DateTime(1993, 1, 1)

for arch in test_architectures, dataset in test_ecco_datasets
    A = typeof(arch)
    D = typeof(dataset)

    if dataset isa ECCO2DarwinMonthly
        @info "Skipping tests because of failure (see https://github.com/CliMA/ClimaOcean.jl/issues/636)"
    else
        @testset "$A metadata tests for $D" begin
            @info "Running Metadata tests for $D on $A..."

            time_resolution = dataset isa ECCO2Daily ? Day(1) : Month(1)
            end_date = start_date + 4 * time_resolution
            dates = start_date : time_resolution : end_date

            @testset "Fields utilities" begin
                for name in test_names[dataset]
                    metadata = Metadata(name; dates, dataset)

                    @show "I am downloading the dataset"
                    download_dataset(metadata) # just in case is not downloaded
                    @show "dataset downloaded"
                    for datum in metadata
                        @test isfile(metadata_path(datum))
                    end

                    datum = first(metadata)
                    ψ = Field(datum, arch, inpainting=NearestNeighborInpainting(2))
                    @test ψ isa Field
                    datapath = ClimaOcean.DataWrangling.inpainted_metadata_path(datum)
                    @test isfile(datapath)
                end
            end

            @testset "Setting a field from a dataset" begin
                test_setting_from_metadata(arch, dataset, start_date, inpainting, 
                                        varnames=test_names[dataset])
            end

            @testset "Field utilities" begin
                test_ocean_metadata_utilities(arch, dataset, dates, inpainting,
                                            varnames=test_names[dataset])
            end

            @testset "DatasetRestoring with LinearlyTaperedPolarMask" begin
                test_dataset_restoring(arch, dataset, dates, inpainting, 
                                    varnames=test_names[dataset],
                                    fldnames=test_fields[dataset])
            end

            @testset "Timestepping with DatasetRestoring" begin
                test_timestepping_with_dataset_restoring(arch, dataset, dates, inpainting, 
                                                        varnames=test_names[dataset], 
                                                        fldnames=test_fields[dataset])
            end

            # @testset "Dataset cycling boundaries" begin
            #     test_cycling_dataset_restoring(arch, dataset, dates, inpainting)
            # end

            # Expensive due to the high resolution of ECCO2
            # @testset "Inpainting algorithm" begin
            #     test_inpainting_algorithm(arch, dataset, start_date, inpainting)
            # end
        end
    end
end
