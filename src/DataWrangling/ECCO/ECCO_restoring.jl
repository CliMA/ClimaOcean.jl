using Oceananigans.Grids: node, on_architecture
using Oceananigans.Fields: interpolate!, interpolate, location, instantiated_location
using Oceananigans.OutputReaders: Cyclical, TotallyInMemory, AbstractInMemoryBackend, FlavorOfFTS, time_indices
using Oceananigans.Utils: Time

using Base

using NCDatasets
using JLD2 

using Dates: Second
using ClimaOcean: stateindex

import Oceananigans.Fields: set!
import Oceananigans.OutputReaders: new_backend, update_field_time_series!

@inline instantiate(T::DataType) = T()
@inline instantiate(T) = T

struct ECCONetCDFBackend{N} <: AbstractInMemoryBackend{Int}
    start :: Int
    length :: Int

    ECCONetCDFBackend{N}(start::Int, length::Int) where N = new{N}(start, length)
end

"""
    ECCONetCDFBackend(length)

Represents an ECCO FieldTimeSeries backed by ECCO native .nc files.
Each time instance is stored in an individual file.
"""
ECCONetCDFBackend(length; on_native_grid = false) = ECCONetCDFBackend{on_native_grid}(1, length)

Base.length(backend::ECCONetCDFBackend)  = backend.length
Base.summary(backend::ECCONetCDFBackend) = string("ECCONetCDFBackend(", backend.start, ", ", backend.length, ")")

const ECCONetCDFFTS{N} = FlavorOfFTS{<:Any, <:Any, <:Any, <:Any, <:ECCONetCDFBackend{N}} where N

new_backend(::ECCONetCDFBackend{N}, start, length) where N = ECCONetCDFBackend{N}(start, length)
on_native_grid(::ECCONetCDFBackend{N}) where N = N

function set!(fts::ECCONetCDFFTS, path::ECCOMetadata=fts.path, name::String=fts.name) 

    backend = fts.backend
    start = backend.start

    for t in start:start+length(backend)-1
        
        # find the file associated with the time index
        metadata = @inbounds path[t] 
        set!(fts[t], metadata)

        # Make sure we clean up after ourselves!
        GC.gc()
    end

    fill_halo_regions!(fts)

    return nothing
end

"""
    ECCO_times(metadata; start_time = metadata.dates[1])

Extracts the time values from the given metadata and calculates the time difference
from the start time.

# Arguments
- `metadata`: The metadata containing the date information.
- `start_time`: The start time for calculating the time difference. Defaults to the first date in the metadata.

# Returns
An array of time differences in seconds.
"""
function ECCO_times(metadata; start_time = first(metadata).dates)
    times = zeros(length(metadata))
    for (t, data) in enumerate(metadata)
        date = data.dates
        time = date - start_time
        time = Second(time).value
        times[t] = time
    end

    return times
end

"""
    ECCO_field_time_series(metadata::ECCOMetadata;
                           architecture = CPU(),
                           time_indices_in_memory = 2,
                           time_indexing = Cyclical(),
                           grid = nothing)

Create a field time series object for ECCO data.

# Arguments:
- metadata: An ECCOMetadata object containing information about the ECCO dataset.

# Keyword Arguments:
- architecture: The architecture to use for computations (default: CPU()).
- time_indices_in_memory: The number of time indices to keep in memory (default: 2).
- time_indexing: The time indexing scheme to use (default: Cyclical()).
- grid: if not a `nothing`, the ECCO data is directly interpolated on the `grid`,
"""
function ECCO_field_time_series(metadata::ECCOMetadata;	
                                architecture = CPU(),	
                                time_indices_in_memory = 2,	
                                time_indexing = Cyclical(),
                                grid = nothing)	

    # ECCO data is too chunky to allow other backends	
    backend = ECCONetCDFBackend(time_indices_in_memory; 
                                on_native_grid = isnothing(grid))

    # Making sure all the required individual files are downloaded
    download_dataset!(metadata)

    location = field_location(metadata)
    ftmp = empty_ECCO_field(first(metadata); architecture)
    shortname = short_name(metadata)

    ECCO_native_grid = ftmp.grid
    boundary_conditions = FieldBoundaryConditions(ECCO_native_grid, location)
    times = ECCO_times(metadata)

    fts_grid = isnothing(grid) ? ECCO_native_grid : grid

    fts = FieldTimeSeries{location...}(fts_grid, times;	
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
    ECCO_field_time_series(ECCOMetadata(variable_name, all_ECCO_dates(version), version); kw...)

# Variable names for restoreable data
struct Temperature end
struct Salinity end
struct UVelocity end
struct VVelocity end

oceananigans_fieldname = Dict(
    :temperature => Temperature(), 
    :salinity    => Salinity(), 
    :u_velocity  => UVelocity(), 
    :v_velocity  => VVelocity())

@inline Base.getindex(fields, i, j, k, ::Temperature) = @inbounds fields.T[i, j, k]
@inline Base.getindex(fields, i, j, k, ::Salinity)    = @inbounds fields.S[i, j, k]
@inline Base.getindex(fields, i, j, k, ::UVelocity)   = @inbounds fields.u[i, j, k]
@inline Base.getindex(fields, i, j, k, ::VVelocity)   = @inbounds fields.v[i, j, k]

"""
    struct ECCORestoring{FTS, G, M, V, N} <: Function

A struct representing ECCO restoring.

# Fields
- `ECCO_fts`: The ECCO FTS on the native ECCO grid.
- `ECCO_grid`: The native ECCO grid to interpolate from.
- `mask`: A mask (could be a number, an array, a function or a field).
- `variable_name`: The variable name of the variable that needs restoring.
- `λ⁻¹`: The reciprocal of the restoring timescale.
"""
struct ECCORestoring{FTS, G, M, V, N} <: Function
    ECCO_fts      :: FTS
    ECCO_grid     :: G
    mask          :: M
    variable_name :: V
    λ⁻¹           :: N
end

Adapt.adapt_structure(to, p::ECCORestoring) = 
    ECCORestoring(Adapt.adapt(to, p.ECCO_fts), 
                  Adapt.adapt(to, p.ECCO_grid),
                  Adapt.adapt(to, p.mask),
                  Adapt.adapt(to, p.variable_name),
                  Adapt.adapt(to, p.λ⁻¹))

@inline function (p::ECCORestoring)(i, j, k, grid, clock, fields)
    
    # Figure out all the inputs: time, location, and node
    time = Time(clock.time)
    loc  = location(p.ECCO_fts)

    # Retrieve the variable to force
    @inbounds var = fields[i, j, k, p.variable_name]

    ECCO_backend = p.ECCO_fts.backend
    native_grid = on_native_grid(ECCO_backend) 

    ECCO_var = get_ECCO_variable(Val(native_grid), p.ECCO_fts, i, j, k, p.ECCO_grid, grid, time)

    # Extracting the mask value at the current node
    mask = stateindex(p.mask, i, j, k, grid, clock.time, loc)

    return p.λ⁻¹ * mask * (ECCO_var - var)
end

# Differentiating between restoring done with an ECCO FTS
# that lives on the native ECCO grid, that requires interpolation in space
# _inside_ the restoring function and restoring based on an ECCO 
# FTS defined on the model grid that requires only time interpolation
@inline function get_ECCO_variable(::Val{true}, ECCO_fts, i, j, k, ECCO_grid, grid, time)
    # Extracting the ECCO field time series data and parameters
    ECCO_times         = ECCO_fts.times
    ECCO_data          = ECCO_fts.data
    ECCO_time_indexing = ECCO_fts.time_indexing
    ECCO_backend       = ECCO_fts.backend
    ECCO_location      = instantiated_location(ECCO_fts)

    X = node(i, j, k, grid, ECCO_location...)

    # Interpolating the ECCO field time series data onto the current node and time
    return interpolate(X, time, ECCO_data, ECCO_location, ECCO_grid, ECCO_times, ECCO_backend, ECCO_time_indexing)
end    

@inline get_ECCO_variable(::Val{false}, ECCO_fts, i, j, k, ECCO_grid, grid, time) = @inbounds ECCO_fts[i, j, k, time]

"""
    ECCO_restoring_forcing(metadata::ECCOMetadata;
                            architecture = CPU(), 
                            backend = ECCONetCDFBackend(2),
                            time_indexing = Cyclical(),
                            mask = 1,
                            timescale = 5days)

Create a restoring forcing term that restores to values stored in an ECCO field time series.

# Arguments:
=============
- `metadata`: The metadata for the ECCO field time series.

# Keyword Arguments:
====================
- `architecture`: The architecture. Typically `CPU` or `GPU`
- `time_indices_in_memory`: The number of time indices to keep in memory. trade-off between performance
                            and memory footprint.    
- `time_indexing`: The time indexing scheme for the field time series, see [`FieldTimeSeries`](@ref)
- `mask`: The mask value. Can be a function of `(x, y, z, time)`, an array or a number
- `timescale`: The restoring timescale.
"""
function ECCO_restoring_forcing(variable_name::Symbol, version=ECCO4Monthly(); kw...) 
     metadata = ECCOMetadata(variable_name, all_ECCO_dates(version), version)
    return ECCO_restoring_forcing(metadata; kw...)
end

function ECCO_restoring_forcing(metadata::ECCOMetadata;
                                architecture = CPU(), 
                                time_indices_in_memory = 2, # Not more than this if we want to use GPU!
                                time_indexing = Cyclical(),
                                mask = 1,
                                timescale = 20days,
                                grid = nothing)

    ECCO_fts  = ECCO_field_time_series(metadata; grid, architecture, time_indices_in_memory, time_indexing)                  
    ECCO_grid = ECCO_fts.grid

    # Grab the correct Oceananigans field to restore
    variable_name = metadata.name
    field_name = oceananigans_fieldname[variable_name]
    ECCO_restoring = ECCORestoring(ECCO_fts, ECCO_grid, mask, field_name, 1 / timescale)
    
    # Defining the forcing that depends on the restoring field.
    restoring_forcing = Forcing(ECCO_restoring; discrete_form = true)

    return restoring_forcing
end

