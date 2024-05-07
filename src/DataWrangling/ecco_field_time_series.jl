using Oceananigans.Units
using Oceananigans.Fields: interpolate!
using Oceananigans.OutputReaders: Cyclical, TotallyInMemory, AbstractInMemoryBackend, FlavorOfFTS, time_indices
using Oceananigans.Utils: Time

using CUDA: @allowscalar

using NCDatasets
using JLD2 
using Dates

using ClimaOcean: stateindex

import Oceananigans.Fields: set!
import Oceananigans.OutputReaders: new_backend, update_field_time_series!

struct ECCONetCDFBackend <: AbstractInMemoryBackend{Int}
    start :: Int
    length :: Int
end

"""
    ECCONetCDFBackend(length)

Represents an ECCO FieldTimeSeries backed by ECCO native .nc files.
"""
ECCONetCDFBackend(length) = ECCONetCDFBackend(1, length)

Base.length(backend::ECCONetCDFBackend) = backend.length
Base.summary(backend::ECCONetCDFBackend) = string("ECCONetCDFBackend(", backend.start, ", ", backend.length, ")")

const ECCONetCDFFTS = FlavorOfFTS{<:Any, <:Any, <:Any, <:Any, <:ECCONetCDFBackend}

new_backend(::ECCONetCDFBackend, start, length) = ECCONetCDFBackend(start, length)

function set!(fts::ECCONetCDFFTS, path::ECCOMetadata=fts.path, name::String=fts.name) 

    backend = fts.backend
    start = backend.start

    for t in start:start+length(backend)-1
        
        # find the file associated with the time index
        metadata = @inbounds path[t] 

        arch = architecture(fts)
        f = inpainted_ecco_field(metadata; architecture = arch, maxiter = 5)
        set!(fts[t], f)
    end

    fill_halo_regions!(fts)

    return nothing
end

"""
    ecco_times(metadata; start_time = metadata.dates[1])

Extracts the time values from the given metadata and calculates the time difference
from the start time.

# Arguments
- `metadata`: The metadata containing the date information.
- `start_time`: The start time for calculating the time difference. Defaults to the first date in the metadata.

# Returns
An array of time differences in seconds.
"""
function ecco_times(metadata; start_time = metadata.dates[1])
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
    ECCO_field_time_series(variable_name;
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
which describes one of the variables in the "ECCO" dataset.
"""
function ECCO_field_time_series(metadata::ECCOMetadata;
                                 architecture = CPU(),
                                 backend = ECCONetCDFBackend(2),
                                 time_indexing = Cyclical())

    # ECCO data is too chunky to allow other backends
    if !(backend isa ECCONetCDFBackend) 
        msg = string("We cannot load the ECCO dataset with an $(backend) backend, only ECCONetCDFBackend is allowed!")
        throw(ArgumentError(msg))
    end

    # Making sure all the required individual files are downloaded
    download_dataset!(metadata)

    location = field_location(metadata)
    ftmp = empty_ecco_field(first(metadata); architecture)
    shortname = short_name(metadata)

    ECCO_native_grid = ftmp.grid
    boundary_conditions = FieldBoundaryConditions(ECCO_native_grid, location)
    times = ecco_times(metadata)
    
    fts = FieldTimeSeries{location...}(ECCO_native_grid, times;
                                       backend,
                                       time_indexing,
                                       boundary_conditions,
                                       path = metadata,
                                       name = shortname)
            
    # Let's set the data
    set!(fts)

    return fts
end

@inline variable_name(i) = ifelse(i == 1, :T, ifelse(i == 2, :S, ifelse(i == 3, :u, :v)))

oceananigans_fieldindex = Dict(
    :temperature => 1,
    :salinity    => 2,
    :u_velocity  => 3,
    :v_velocity  => 4
)

"""
    struct ECCORestoring{FTS, I, M, N}

A struct representing the ECCO forcing function for a specific field.

## Fields
- `ecco_fts`: The ECCO field time series data.
- `field_idx`: The index of the field.
- `mask`: The mask value.
- `λ`: The timescale.

"""
struct ECCORestoring{FTS, I, M, N}
    ecco_fts  :: FTS
    field_idx :: I
    mask      :: M
    λ         :: N
end

@inline function (p::ECCORestoring)(i, j, k, grid, clock, fields)
    
    # Figure out all the inputs: variable name, time, location, and node
    var_name  = variable_name(p.field_idx)
    time      = Time(clock.time)
    loc       = location(p.ecco_fts)
    X         = node(i, j, k, grid, loc)

    # Extracting the ECCO field time series data
    ecco_grid          = p.ecco_fts.grid
    ecco_data          = p.ecco_fts.data
    ecco_backend       = p.ecco_fts.backend
    ecco_time_indexing = p.ecco_fts.time_indexing

    # Extracting the field value at the current node
    @inbounds var = getproperty(fields, var_name)[i, j, k]
    
    # Interpolating the ECCO field time series data
    ecco_var = interpolate(X, time, ecco_data, loc, grid, ecco_grid, ecco_backend, ecco_time_indexing)
    
    # Extracting the mask value at the current node
    mask = stateindex(p.mask, i, j, k, grid, clock.time, loc)

    return 1 / p.λ * mask * (ecco_var - var)
end

"""
    ECCO_restoring_forcing(metadata::ECCOMetadata;
                            architecture = CPU(), 
                            backend = ECCONetCDFBackend(2),
                            time_indexing = Cyclical(),
                            mask = 1,
                            timescale = 5days)

Create a restoring forcing term for ECCO field time series.

## Arguments
- `metadata`: The metadata for the ECCO field time series.
- `architecture`: The architecture.
- `backend`: The backend.
- `time_indexing`: The time indexing.
- `mask`: The mask value.
- `timescale`: The timescale.

## Returns
- The restoring forcing term.

"""
function ECCO_restoring_forcing(metadata::ECCOMetadata;
                                 architecture = CPU(), 
                                 backend = ECCONetCDFBackend(2),
                                 time_indexing = Cyclical(),
                                 mask = 1,
                                 timescale = 5days)

    ecco_fts = ECCO_field_time_series(metadata; architecture, backend, time_indexing)                  

    variable_name = metadata.name
    field_idx = oceananigans_fieldindex[variable_name]
    ecco_restoring = ECCORestoring(ecco_fts, field_idx, mask, timescale)

    restoring_forcing = Forcing(ecco_restoring; 
                                discrete_form=true, 
                                parameters=(; ecco_fts, field_idx, mask, λ=timescale))

    return restoring_forcing
end