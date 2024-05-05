using CFTime
using Dates

function generate_ECCO_restoring_dataset(variable_name; 
                                         architecture = CPU(),
                                         dataset_filename = "ciao.nc",
                                         start_date = DateTimeAllLeap(1992, 1, 2),
                                         end_date = DateTimeAllLeap(2017, 12, 31))

    Nt = 0
    date = start_date                                     
    while date <= end_date
        Nt += 1
        date = date + Day(1)
    end

    field = empty_ecco4_field(variable_name)

    # Open a new dataset
    ds = Dataset(dataset_filename, "c")

    # Define the dimension "lon", "lat", "z" and "time",
    # with the size ECCO_Nx, ECCO_Ny, ECCO_Nz, and Nt, respectively.
    defDim(ds, "lon",  ECCO_Nx)
    defDim(ds, "lat",  ECCO_Ny)
    defDim(ds, "z",    ECCO_Nz)
    defDim(ds, "time", Nt)

    # Define the variables with the attribute units
    data = defVar(ds, string(variable_name), Float32, ("lon", "lat", "z", "time"))
    time = defVar(ds, "time", Float32, ("time",))
    lat  = defVar(ds, "lat",  Float32, ("lat",))
    lon  = defVar(ds, "lon",  Float32, ("lon",))
    zvar = defVar(ds, "z",    Float32, ("z", ))

    # Fill in the spatial nodes
    λ, φ, z = nodes(field)
    lat  .= λ
    lon  .= φ
    zvar .= z

    # Start filling in time-specific data
    n = 1
    date = start_date                                     
    while date <= end_date
        metadata = ECCOMetadata(variable_name, date)

        @info "retrieving inpainted $(variable_name) at $(date)"
        field = inpainted_ecco4_field(metadata; architecture, filename)
        copyto!(data[:, :, :, n], interior(field))
        time[n] = date.instant.periods.value / 1000 # Converting milliseconds to seconds
    
        n += 1

        date = date + Day(1)
    end

    close(ds)

    return nothing
end

function ECCO_forcing(variable_name; 
                      time_indices = :,  
                      mask = nothing,
                      architecture = CPU(),
                      backend = NetCDFBackend(4))

    restoring_field = reanalysis_field_time_series(variable_name; time_indices, backend)

    return restoring_field
end