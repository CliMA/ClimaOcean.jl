include("runtests_setup.jl")

@testset "Availability of JRA55 data" begin
    @info "Testing that we can download all the JRA55 data..."
    for name in ClimaOcean.DataWrangling.JRA55.JRA55_variable_names
        fts = ClimaOcean.JRA55.JRA55_field_time_series(name; time_indices=2:3)
    end
end
