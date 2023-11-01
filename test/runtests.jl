# Common test setup file to make stand-alone tests easy
include("runtests_setup.jl")

@test 1==1

# Tests JRA55 utilities, plus some DataWrangling utilities
include("test_jra55.jl")

