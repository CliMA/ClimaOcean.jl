# Common test setup file to make stand-alone tests easy
include("runtests_setup.jl")

test_group = get(ENV, "TEST_GROUP", :all)
test_group = Symbol(test_group)

# Fictitious grid that triggers bathymetry download
function download_bathymetry()
    grid = LatitudeLongitudeGrid(size = (10, 10, 1), 
                                 longitude = (0, 100), 
                                 latitude = (0, 50),
                                 z = (-6000, 0))

    bottom = regrid_bathymetry(grid)

    return nothing
end

if test_group == :init || test_group == :all
    using CUDA
    CUDA.precompile_runtime()

    ####
    #### Download bathymetry data
    ####
    
    download_bathymetry() 

    # Download JRA55 data
    atmosphere = JRA55_prescribed_atmosphere()

    # Download ECCO data
    start_date = DateTimeProlepticGregorian(1993, 1, 1)
    end_date = DateTimeProlepticGregorian(1993, 4, 1)
    dates = start_date : Month(1) : end_date

    temperature = ECCOMetadata(:temperature, dates)
    salinity = ECCOMetadata(:salinity, dates)

    download_dataset!(temperature)
    download_dataset!(salinity)
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

if test_group == :bathymetry || test_group == :all
    include("test_bathymetry.jl")
end

if test_group == :simulations || test_group == :all
    include("test_simulations.jl")
end

