module JRA55

using Oceananigans
using NCDatasets

filenames = Dict(
    :shortwave_radiation => "RYF.rsds.1990_1991.nc",
)

variable_names = Dict(
    :shortwave_radiation => "rsds",
)

urls = Dict(
    :shortwave_radiation => "https://www.dropbox.com/scl/fi/z6fkvmd9oe3ycmaxta131/" *
                            "RYF.rsds.1990_1991.nc?rlkey=r7q6zcbj6a4fxsq0f8th7c4tc&dl=0"
)

"""
    jra55_field_time_series(name, architecture=CPU();
                            time_indices = :,    
                            url = urls[name],
                            filename = filenames[name],
                            variable_name = variable_names[name])

Return a FieldTimeSeries representing JRA55 data at the interface between the
atmosphere and ocean.
"""
function jra55_field_time_series(name, arch=CPU();
                                 time_indices = :,    
                                 url = urls[name],
                                 filename = filenames[name],
                                 variable_name = variable_names[name])

    isfile(filename) || download(url, filename)

    ds = Dataset(filename)
    data = ds[variable_name][:, :, time_indices]

    # Make source field
    λ = ds["lon_bnds"][1, :]
    φ = ds["lat_bnds"][1, :]
    times = ds["time"][time_indices]
    close(ds)

    push!(φ, 90)
    push!(λ, λ[1] + 360)

    Nxs = length(λ) - 1
    Nys = length(φ) - 1

    grid = LatitudeLongitudeGrid(arch,
                                 size = (Nxs, Nys);
                                 longitude = λ,
                                 latitude = φ,
                                 topology = (Periodic, Bounded, Flat))

    fts = FieldTimeSeries{Center, Center, Nothing}(grid, times)

    # JRA55 fields have 1 grid point in the z-direction
    interior(fts, :, :, 1, :) .= data[:, :, :]

    return fts
end

end # module

