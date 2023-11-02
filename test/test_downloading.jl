include("runtests_setup.jl")

@testset "Availability of JRA55 data" begin
    @info "Testing that we can download all the JRA55 data..."
    for name in ClimaOcean.JRA55.jra55_variable_names 
        fts = ClimaOcean.JRA55.jra55_field_time_series(name, CPU(); time_indices=1:1)
    end
end

