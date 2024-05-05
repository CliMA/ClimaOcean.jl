using CFTime
using Dates

import Base: reverse

@inline reverse(::Face)   = Center()
@inline reverse(::Center) = Face()

@inline instantiate(T::DataType) = T()
@inline instantiate(T) = T

function download_dataset!(metadata::ECCOMetadata)
    dates = metadata.dates
    output_filename = file_name(metadata)
    isfile(output_filename) ||  
        generate_ECCO_restoring_dataset!(metadata.name; dates, output_filename)

    return nothing
end

const all_ecco_4_dates = DateTimeAllLeap(1992, 1, 2) : Day(1) : DateTimeAllLeap(2017, 12, 31)

function generate_ECCO_restoring_dataset!(variable_name; 
                                          architecture = CPU(),
                                          dates = all_ecco_4_dates,
                                          output_filename = "ecco4_$(metadata.name).nc")

    Nt    = length(dates)
    field = empty_ecco4_field(variable_name)

    # Open a new dataset
    ds = Dataset(output_filename, "c")

    # Define the dimension "lon", "lat", "z" and "time",
    # with the size ECCO_Nx, ECCO_Ny, ECCO_Nz, and Nt, respectively.
    defDim(ds, "lon",  ECCO_Nx)
    defDim(ds, "lat",  ECCO_Ny)
    defDim(ds, "z",    ECCO_Nz)    
    defDim(ds, "lon_bnds",  ECCO_Nx)
    defDim(ds, "lat_bnds",  ECCO_Ny)
    defDim(ds, "z_bnds",    ECCO_Nz + 1)
    defDim(ds, "time", Nt)

    # Define the variables with the attribute units
    data = defVar(ds, string(variable_name), Float32, ("lon", "lat", "z", "time"))
    time = defVar(ds, "time", Float32, ("time",))
    lon  = defVar(ds, "lon",  Float32, ("lon",))
    lat  = defVar(ds, "lat",  Float32, ("lat",))
    zn   = defVar(ds, "z",    Float32, ("z", ))
    
    lon_bounds  = defVar(ds, "lon_bnds", Float32, ("lon_bnds",))
    lat_bounds  = defVar(ds, "lat_bnds", Float32, ("lat_bnds",))
    z_bounds    = defVar(ds, "z_bnds",   Float32, ("z_bnds", ))

    # Fill in the spatial nodes
    loc = instantiate.(location(field))
    λ,  φ,  z = nodes(field)
    λr, φr, zr = nodes(field.grid, reverse.(loc)...)
    lon  .= λ
    lat  .= φ
    zn   .= z
    lon_bounds .= λr
    lat_bounds .= φr[1:end-1]
    z_bounds   .= zr

    # Start filling in time-specific data
    for (n, date) in enumerate(dates)

        metadata = ECCOMetadata(variable_name, date)
        
        @info "retrieving inpainted $(variable_name) at $(date)"
        field = inpainted_ecco4_field(metadata; architecture)
        data[:, :, :, n] = Float32.(Array(interior(field)))
        time[n] = date.instant.periods.value / 1000 # Converting milliseconds to seconds
    end

    close(ds)

    return nothing
end