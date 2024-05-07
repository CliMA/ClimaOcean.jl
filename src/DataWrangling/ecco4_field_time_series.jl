using Oceananigans.Fields: interpolate!
using Oceananigans.OutputReaders: Cyclical, TotallyInMemory, AbstractInMemoryBackend, FlavorOfFTS, time_indices

using CUDA: @allowscalar

using NCDatasets
using JLD2 
using Dates

import Oceananigans.Fields: set!
import Oceananigans.OutputReaders: new_backend, update_field_time_series!

struct ECCO4NetCDFBackend <: AbstractInMemoryBackend{Int}
    start :: Int
    length :: Int
end

"""
    ECCO4NetCDFBackend(length)

Represents an ECCO4 FieldTimeSeries backed by ECCO4 native .nc files.
"""
ECCO4NetCDFBackend(length) = ECCO4NetCDFBackend(1, length)

Base.length(backend::ECCO4NetCDFBackend) = backend.length
Base.summary(backend::ECCO4NetCDFBackend) = string("ECCO4NetCDFBackend(", backend.start, ", ", backend.length, ")")

const ECCO4NetCDFFTS = FlavorOfFTS{<:Any, <:Any, <:Any, <:Any, <:ECCO4NetCDFBackend}

new_backend(::ECCO4NetCDFBackend, start, length) = ECCO4NetCDFBackend(start, length)

function set!(fts::ECCO4NetCDFFTS, path::ECCOMetadata=fts.path, name::String=fts.name) 

    backend = fts.backend
    start = backend.start

    for t in start:start+length(backend)-1
        
        # find the file associated with the time index
        metadata = @inbounds path[t] 

        arch = architecture(fts)
        f = inpainted_ecco4_field(metadata; architecture = arch, maxiter = 10)
        set!(fts[t], f)
    end

    fill_halo_regions!(fts)

    return nothing
end

"""
    download_dataset!(metadata::ECCOMetadata)

Download the dataset specified by the given metadata. If the metadata contains a single date, 
the dataset is downloaded directly. If the metadata contains multiple dates, the dataset is 
downloaded for each date individually.

# Arguments
- `metadata::ECCOMetadata`: The metadata specifying the dataset to be downloaded.
"""
function download_dataset!(metadata::ECCOMetadata)
    
    for data in metadata
        filename = file_name(data)

        if !isfile(filename) 
            cmd = `podaac-data-downloader -c $(remote_folder) -d ./ --start-date $(datestr)T00:00:00Z --end-date $(datestr)T00:00:00Z -e .nc`
            @info "downloading $(filename) from $(remote_folder)"
            try
                run(cmd)
            catch error
                @info "Note: to download ECCO4 data please install podaac-data-downloader using \\ 
                    `pip install podaac`. Provide a username and password to the python environment. \\
                    For details about the installation refer to "
                throw(ArgumentError("The error is $error"))
            end
        else
            @info "File $(filename) already exists"
        end
    end

    return nothing
end

"""
    ecco4_times(metadata; start_time = metadata.dates[1])

Extracts the time values from the given metadata and calculates the time difference
from the start time.

# Arguments
- `metadata`: The metadata containing the date information.
- `start_time`: The start time for calculating the time difference. Defaults to the first date in the metadata.

# Returns
An array of time differences in seconds.
"""
function ecco4_times(metadata; start_time = metadata.dates[1])
    times = []
    for data in metadata
        date = data.dates
        time = date - start_time
        time = Second(time).value
        push!(times, time)
    end

    return times
end

"""
    ECCO4_field_time_series(variable_name;
                            architecture = CPU(),
                            location = nothing,
                            url = nothing,
                            filename = nothing,
                            shortname = nothing,
                            backend = InMemory(),
                            preprocess_chunk_size = 10,
                            preprocess_architecture = CPU(),
                            time_indices = nothing)

Return a `FieldTimeSeries` containing oceanic reanalysis data for `variable_name`,
which describes one of the variables in the "ECCO4" dataset.
"""
function ECCO4_field_time_series(metadata::ECCOMetadata;
                                 architecture = CPU(),
                                 backend = ECCO4NetCDFBackend(2),
                                 time_indexing = Cyclical())

    # ECCO4 data is too chunky to allow other backends
    if !(backend isa ECCO4NetCDFBackend) 
        msg = string("We cannot load the ECCO4 dataset with an $(backend) backend, only ECCO4NetCDFBackend is allowed!")
        throw(ArgumentError(msg))
    end

    # Making sure all the required individual files are downloaded
    download_dataset!(metadata)

    location = field_location(metadata)
    ftmp = empty_ecco4_field(first(metadata); architecture)
    shortname = short_name(metadata)

    ECCO4_native_grid = ftmp.grid
    boundary_conditions = FieldBoundaryConditions(ECCO4_native_grid, location)
    times = ecco4_times(metadata)
    
    fts = FieldTimeSeries{location...}(ECCO4_native_grid, times;
                                       backend,
                                       time_indexing,
                                       boundary_conditions,
                                       path = metadata,
                                       name = shortname)
            
    # Let's set the data
    set!(fts)

    return fts
end