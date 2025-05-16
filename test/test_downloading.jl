include("runtests_setup.jl")

using ClimaOcean.DataWrangling: metadata_path
using ClimaOcean.JRA55: JRA55NetCDFBackend

@testset "Availability of JRA55 data" begin
    @info "Testing that we can download all the JRA55 data..."
    for name in ClimaOcean.DataWrangling.JRA55.JRA55_variable_names
        fts = ClimaOcean.JRA55.JRA55FieldTimeSeries(name; backend=JRA55NetCDFBackend(2))
    end
end

@testset "Availability of ECCO/EN4 data" begin
    for dataset in test_datasets

        @info "Testing that we can download $(typeof(dataset)) data..."

        variables = dataset isa ECCO4Monthly ? keys(ClimaOcean.ECCO.ECCO4_dataset_variable_names) :
                    dataset isa ECCO2Monthly ? keys(ClimaOcean.ECCO.ECCO2_dataset_variable_names) :
                    dataset isa ECCO2Daily ? keys(ClimaOcean.ECCO.ECCO2_dataset_variable_names) :
                    dataset isa EN4Monthly ? keys(ClimaOcean.EN4.EN4_dataset_variable_names) :
                    error("what am I supposed to download?")

        for variable in variables
            metadata = Metadata(variable; dates=DateTimeProlepticGregorian(1993, 1, 1), dataset)
            filepath = metadata_path(metadata)
            isfile(filepath) && rm(filepath; force=true)
            ClimaOcean.DataWrangling.download_dataset(metadata)
        end
    end
end

@testset "Availability of the ETOPO2022 Bathymetry" begin
    @info "Testing that we can download the bathymetry..."
    ETOPOmetadata = Metadatum(:bottom_height, dataset=ETOPO2022())
    filepath = metadata_path(metadata)
    isfile(filepath) && rm(filepath; force=true)
    ClimaOcean.DataWrangling.download_dataset(metadata)
end
