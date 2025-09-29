using PythonCall
using CondaPkg
using Oceananigans.Grids: AbstractGrid
using SparseArrays
using LinearAlgebra

const Grids = Union{SpeedyWeather.SpectralGrid, AbstractGrid}

struct BilinearInterpolator{W1, W2}
    set1 :: W1
    set2 :: W2
end

function BilinearInterpolator(grid1::Grids, grid2::Grids) 
    W1 = regridder_weights(grid1, grid2; method="bilinear")
    W2 = regridder_weights(grid2, grid1; method="bilinear")
    return BilinearInterpolator(W1, W2)
end

regrid!(dst, weights, scr) = LinearAlgebra.mul!(vec(dst), weights, vec(scr))

function install_package(package)
    @info "Installing the $package package..."
    CondaPkg.add(package; channel = "conda-forge")
    return nothing
end

# Import and store as constants for submodules
function get_package(package)
    try
        return pyimport(package)
    catch
        install_package(package)
        return pyimport(package)
    end
end

two_dimensionalize(lat::Matrix, lon::Matrix) = lat, lon

function two_dimensionalize(lat::AbstractVector, lon::AbstractVector) 
    Nx = length(lon)
    Ny = length(lat)
    lat = repeat(lat', Nx)
    lon = repeat(lon, 1, Ny)
    return lat, lon
end

coordinate_dataset(grid::SpeedyWeather.SpectralGrid) = coordinate_dataset(grid.grid)
    
function coordinate_dataset(grid::SpeedyWeather.RingGrids.AbstractGrid)
    numpy = get_package("numpy")
    lon, lat = RingGrids.get_londlatds(grid)
    return numpy.array(Dict("lat" => lat, "lon" => lon))
end

function coordinate_dataset(grid::SpeedyWeather.RingGrids.AbstractFullGrid)
    lon  = RingGrids.get_lond(grid)
    lat  = RingGrids.get_latd(grid)
    dlon = lon[2] - lon[1]

    lat_b = [90, 0.5 .* (lat[1:end-1] .+ lat[2:end])..., -90]
    lon_b = [lon[1] - dlon / 2, lon .+ dlon / 2...]

    lat,   lon   = two_dimensionalize(lat,   lon)
    lat_b, lon_b = two_dimensionalize(lat_b, lon_b)

    return structured_coordinate_dataset(lat, lon, lat_b, lon_b)
end

function coordinate_dataset(grid::AbstractGrid)
    lat = Array(φnodes(grid, Center(), Center(), Center()))
    lon = Array(λnodes(grid, Center(), Center(), Center()))

    lat_b = Array(φnodes(grid, Face(), Face(), Center()))
    lon_b = Array(λnodes(grid, Face(), Face(), Center()))

    lat,   lon   = two_dimensionalize(lat,   lon)
    lat_b, lon_b = two_dimensionalize(lat_b, lon_b)

    return structured_coordinate_dataset(lat, lon, lat_b, lon_b)
end

function structured_coordinate_dataset(lat, lon, lat_b, lon_b)
    numpy  = get_package("numpy")
    xarray = get_package("xarray")

    lat   = Array(lat')
    lon   = Array(lon')
    lat_b = Array(lat_b')
    lon_b = Array(lon_b')

    lat = numpy.array(lat)
    lon = numpy.array(lon)

    lat_b = numpy.array(lat_b)
    lon_b = numpy.array(lon_b)

    ds_lat = xarray.DataArray(
        lat,
        dims=["y", "x"],
        coords=Dict(
            "lat" => (["y", "x"], lat),
            "lon" => (["y", "x"], lon)
        ),
        name="latitude"
    )
    
    ds_lon = xarray.DataArray(
        lon,
        dims=["y", "x"],
        coords=Dict(
            "lat" => (["y", "x"], lat),
            "lon" => (["y", "x"], lon)
        ),
        name="longitude"
    )
    
    ds_lat_b = xarray.DataArray(
        lat_b,
        dims=["y_b", "x_b"],
        coords=Dict(
            "lat_b" => (["y_b", "x_b"], lat_b),
            "lon_b" => (["y_b", "x_b"], lon_b)
        ),
    )

    ds_lon_b = xarray.DataArray(
        lon_b,
        dims=["y_b", "x_b"],
        coords=Dict(
            "lat_b" => (["y_b", "x_b"], lat_b),
            "lon_b" => (["y_b", "x_b"], lon_b)
        ),
    )

    return xarray.Dataset(
        Dict("lat"   => ds_lat, 
             "lon"   => ds_lon,
             "lat_b" => ds_lat_b,
             "lon_b" => ds_lon_b)
    )
end

function regridder_weights(dst::Grids, src::Grids; method::String="bilinear")

    src_ds = coordinate_dataset(src)
    dst_ds = coordinate_dataset(dst)

    regridder =  get_package("xesmf").Regridder(src_ds, dst_ds, method, periodic=true)

    # Move back to Julia
    # Convert the regridder weights to a Julia sparse matrix
    coords = regridder.weights.data
    shape  = pyconvert(Tuple{Int, Int}, coords.shape)
    vals   = pyconvert(Array{Float64}, coords.data)
    coords = pyconvert(Array{Float64}, coords.coords)
    rows = coords[1,:].+1
    cols = coords[2,:].+1

    W = sparse(rows, cols, vals, shape[1], shape[2])

    return W
end
