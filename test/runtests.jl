# Common test setup file to make stand-alone tests easy
include("runtests_setup.jl")

test_group = get(ENV, "TEST_GROUP", :all)
test_group = Symbol(test_group)

# Tests JRA55 utilities, plus some DataWrangling utilities
if test_group == :jra55 || test_group == :all
    include("test_jra55.jl")
end

if test_group == :ecco4 || test_group == :all
    include("test_ecco4.jl")
end

# Tests that we can download JRA55 utilities
if test_group == :downloading || test_group == :all
    include("test_downloading.jl")
end
