# Common test setup file to make stand-alone tests easy
include("runtests_setup.jl")

test_group = get(ENV, "TEST_GROUP", :all)
test_group = Symbol(test_group)

if test_group == :set_runtime_version
    using CUDA
    CUDA.set_runtime_version!(v"12.2")
end

if test_group == :precompile_runtime
    using CUDA
    CUDA.precompile_runtime()
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

if test_group == :turbulent_fluxes || test_group == :all
    include("test_surface_fluxes.jl")
end

if test_group == :simulations || test_group == :all
    include("test_simulations.jl")
end

