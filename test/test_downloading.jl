include("runtests_setup.jl")

using ClimaOcean.ECCO: metadata_path

@testset "Availability of JRA55 data" begin
    @info "Testing that we can download all the JRA55 data..."
    for name in ClimaOcean.DataWrangling.JRA55.JRA55_variable_names
        
        fts = ClimaOcean.JRA55.JRA55_field_time_series(name; time_indices=2:3)
    end
end

@testset "Availability of ECCO data" begin
    @info "Testing that we can download ECCO data..."
    for variable in keys(ClimaOcean.ECCO.ECCO4_short_names)
        metadata = ECCOMetadata(variable)
        filepath = metadata_path(metadata)
        isfile(filepath) && rm(filepath; force=true)
        ClimaOcean.ECCO.download_dataset(metadata)
    end
end

@testset "Availability of the Bathymetry" begin
    @info "Testing that we can download the bathymetry..."
    dir="./"
    filename="ETOPO_2022_v1_60s_N90W180_surface.nc"
    filepath=joinpath(dir, filename)
    isfile(filepath) && rm(filepath; force=true)
    download_bathymetry(; dir, filename)
end
