using Oceananigans.Units
using Oceananigans.Grids: node
using Oceananigans.Fields: interpolate!, interpolate, location
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
Each time instance is stored in an individual file.
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
        f = inpainted_ecco_field(metadata; architecture = arch)
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
function ecco_times(metadata; start_time = first(metadata).dates)
    times = []
    for data in metadata
        date = data.dates
        time = date - start_time
        time = Second(time).value
        push!(times, time)
    end

    times = tuple(times...)

    return times
end

"""
        ECCO_field_time_series(metadata::ECCOMetadata;
                            architecture = CPU(),
                            time_indices_in_memory = 2,
                            time_indexing = Cyclical())

Create a field time series object for ECCO data.

Args:
- metadata: An ECCOMetadata object containing information about the ECCO dataset.
- architecture: The architecture to use for computations (default: CPU()).
- time_indices_in_memory: The number of time indices to keep in memory (default: 2).
- time_indexing: The time indexing scheme to use (default: Cyclical()).

Returns:
- fts: A FieldTimeSeries object representing the ECCO field time series.

Example:
```
metadata = ECCOMetadata(...)
fts = ECCO_field_time_series(metadata)
```
"""
function ECCO_field_time_series(metadata::ECCOMetadata;	
                                 architecture = CPU(),	
                                 time_indices_in_memory = 2,	
                                 time_indexing = Cyclical())	

    # ECCO data is too chunky to allow other backends	
    backend = ECCONetCDFBackend(time_indices_in_memory)

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

ECCO_field_time_series(variable_name::Symbol, version=ECCO4Monthly(); kw...) = 
    ECCO_field_time_series(ECCOMetadata(variable_name, all_ecco_dates(version), version); kw...)

oceananigans_fieldname = Dict(
    :temperature => :T, 
    :salinity    => :S, 
    :u_velocity  => :u, 
    :v_velocity  => :v)

"""
    struct ECCORestoring{FTS, I, M, N}

A struct representing the ECCO forcing function for a specific field.

## Fields
- `ecco_fts`: The ECCO field time series data.
- `field_idx`: The index of the field.
- `mask`: The mask value.
- `λ`: The timescale.

"""
struct ECCORestoring{FTS, M, N}
    ecco_fts :: FTS
    mask     :: M
    λ        :: N
end

struct ECCOGPURestoring{D, L, T, G, B, I, M, N}
    ecco_data           :: D
    location            :: L
    ecco_times          :: T
    ecco_grid           :: G
    ecco_backend        :: B
    ecco_time_indexing  :: I
    mask                :: M
    λ                   :: N
end

Adapt.adapt_structure(to, p::ECCORestoring) = 
    ECCOGPURestoring(Adapt.adapt(to, p.ecco_fts.data), 
                     Adapt.adapt(to, instantiated_location(p.ecco_fts)),
                     Adapt.adapt(to, p.ecco_fts.times),
                     Adapt.adapt(to, p.ecco_fts.grid),
                     Adapt.adapt(to, p.ecco_fts.backend),
                     Adapt.adapt(to, p.ecco_fts.time_indexing),
                     Adapt.adapt(to, p.mask),
                     Adapt.adapt(to, p.λ))

@inline function (p::ECCORestoring)(x, y, z, t, var)
    
    # Figure out all the inputs: time, location, and node
    time = Time(t)
    loc  = instantiated_location(p.ecco_fts)
    X    = (x, y, z)

    # Extracting the ECCO field time series data and parameters
    ecco_times         = p.ecco_fts.times
    ecco_grid          = p.ecco_fts.grid
    ecco_data          = p.ecco_fts.data
    ecco_backend       = p.ecco_fts.backend
    ecco_time_indexing = p.ecco_fts.time_indexing

    # Interpolating the ECCO field time series data ont the current node and time
    ecco_var = interpolate(X, time, ecco_data, loc, ecco_grid, ecco_times, ecco_backend, ecco_time_indexing)
    
    # Extracting the mask value at the current node
    # TODO: make this work with numbers, arrays, or functions
    # (something like the reverse of `stateindex`)
    mask = p.mask(X...)

    return 1 / p.λ * mask * (ecco_var - var)
end

@inline function (p::ECCOGPURestoring)(x, y, z, t, var)
    
    # Figure out all the inputs: time, location, and node
    time = Time(t)
    loc  = p.location
    X    = (x, y, z)

    # Extracting the ECCO field time series data and parameters
    ecco_times         = p.ecco_times
    ecco_grid          = p.ecco_grid
    ecco_data          = p.ecco_data
    ecco_backend       = p.ecco_backend
    ecco_time_indexing = p.ecco_time_indexing

    # Interpolating the ECCO field time series data ont the current node and time
    ecco_var = interpolate(X, time, ecco_data, loc, ecco_grid, ecco_times, ecco_backend, ecco_time_indexing)
    
    # Extracting the mask value at the current node
    # TODO: make this work with numbers, arrays, or functions
    # (something like the reverse of `stateindex`)
    mask = p.mask(X...)

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
function ECCO_restoring_forcing(variable_name::Symbol, version=ECCO4Monthly(); kw...) 
     metadata = ECCOMetadata(variable_name, all_ecco_dates(version), version)
    return ECCO_restoring_forcing(metadata; kw...)
end

function ECCO_restoring_forcing(metadata::ECCOMetadata;
                                architecture = CPU(), 
                                time_indices_in_memory = 2, # Not more than this if we want to use GPU!
                                time_indexing = Cyclical(),
                                mask = 1,
                                timescale = 5days)

    ecco_fts = ECCO_field_time_series(metadata; architecture, time_indices_in_memory, time_indexing)                  

    variable_name = metadata.name
    field_name = oceananigans_fieldname[variable_name]
    
    ecco_restoring = ECCORestoring(ecco_fts, mask, timescale)
    # Defining the forcing that depends on the restoring field.
    restoring_forcing = Forcing(ecco_restoring; field_dependencies = field_name)

    return restoring_forcing
end