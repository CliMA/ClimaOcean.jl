
function write_timeseries!(data, loc, grid, times, path, name)

    dims = length(size(data)) - 1
    spatial_indices = Tuple(Colon() for i in 1:dims)

    f_tmp = Field{loc...}(grid)
    fts_tmp = FieldTimeSeries(loc, grid, times; 
                              backend = OnDisk(),
                              path,
                              name)

    for t in eachindex(times)
        set!(f_tmp, data[spatial_indices..., t])                           
        set!(fts_tmp, f_tmp, t)
    end

    return nothing
end