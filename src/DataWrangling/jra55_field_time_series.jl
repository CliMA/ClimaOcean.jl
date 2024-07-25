using Oceananigans.Units
using Oceananigans.Grids: node, on_architecture
using Oceananigans.Fields: interpolate!, interpolate, location, instantiated_location
using Oceananigans.OutputReaders: Cyclical, TotallyInMemory, AbstractInMemoryBackend, FlavorOfFTS, time_indices
using Oceananigans.Utils: Time

using CUDA: @allowscalar
using Base

using NCDatasets
using JLD2 
using Dates

using ClimaOcean: stateindex
using ClimaOcean.DataWrangling: native_times

import Oceananigans.Fields: set!
import Oceananigans.OutputReaders: new_backend, update_field_time_series!

@inline instantiate(T::DataType) = T()
@inline instantiate(T) = T

struct JRA55NetCDFBackend{N} <: AbstractInMemoryBackend{Int}
    start :: Int
    length :: Int

    JRA55NetCDFBackend{N}(start::Int, length::Int) where N = new{N}(start, length)
end

"""
    JRA55NetCDFBackend(length)

Represents an JRA55 FieldTimeSeries backed by JRA55 native .nc files.
Each time instance is stored in an individual file.
"""
JRA55NetCDFBackend(length; on_native_grid = false) = JRA55NetCDFBackend{on_native_grid}(1, length)

Base.length(backend::JRA55NetCDFBackend)  = backend.length
Base.summary(backend::JRA55NetCDFBackend) = string("JRA55NetCDFBackend(", backend.start, ", ", backend.length, ")")

const JRA55NetCDFFTS{N} = FlavorOfFTS{<:Any, <:Any, <:Any, <:Any, <: JRA55NetCDFBackend{N}} where N

new_backend(::JRA55NetCDFBackend{N}, start, length) where N = JRA55ONetCDFBackend{N}(start, length)
on_native_grid(::JRA55NetCDFBackend{N}) where N = N

function set!(fts::JRA55NetCDFFTS, path::JRA55Metadata=fts.path, name::String=fts.name) 

    backend = fts.backend
    start = backend.start

    # Set the JRA55 dataset based on the backend!s

    for t in start:start+length(backend)-1
        
        # find the file associated with the time index
        metadata = @inbounds path[t] 

        arch = architecture(fts)

        # f = inpainted_ecco_field(metadata; architecture = arch)
        if on_native_grid(backend)
            set!(fts[t], f)
        else
            interpolate!(fts[t], f)
        end
    end

    fill_halo_regions!(fts)

    return nothing
end

"""
    JRA55_field_time_series(metadata::ECCOMetadata;
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
function JRA55_field_time_series(metadata::JRA55Metadata;	
                                 architecture = CPU(),	
                                 backend = JRA55NetCDFBackend(20),	
                                 time_indexing = Cyclical(),
                                 grid = nothing)	

    # ECCO data is too chunky to allow other backends	
    backend = if backend isa JRA55NetCDFBackend
        JRA55NetCDFBackend(backend.length; 
                           on_native_grid = isnothing(grid))
    else 
        backend
    end

    # Making sure all the required individual files are downloaded
    download_dataset!(metadata)

    location  = field_location(metadata)
    shortname = short_name(metadata)

    JRA55_native_grid = jra55_native_grid()
    boundary_conditions = FieldBoundaryConditions(JRA55_native_grid, location)
    times = native_times(metadata)

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
    ECCO_field_time_series(ECCOMetadata(variable_name, all_ecco_dates(version), version); kw...)

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
- `ecco_fts`: The ECCO FTS on the native ECCO grid.
- `ecco_grid`: The native ECCO grid to interpolate from.
- `mask`: A mask (could be a number, an array, a function or a field).
- `variable_name`: The variable name of the variable that needs restoring.
- `λ⁻¹`: The reciprocal of the restoring timescale.
"""
struct ECCORestoring{FTS, G, M, V, N} <: Function
    ecco_fts  :: FTS
    ecco_grid :: G
    mask      :: M
    variable_name   :: V
    λ⁻¹       :: N
end

Adapt.adapt_structure(to, p::ECCORestoring) = 
    ECCORestoring(Adapt.adapt(to, p.ecco_fts), 
                  Adapt.adapt(to, p.ecco_grid),
                  Adapt.adapt(to, p.mask),
                  Adapt.adapt(to, p.variable_name),
                  Adapt.adapt(to, p.λ⁻¹))

@inline function (p::ECCORestoring)(i, j, k, grid, clock, fields)
    
    # Figure out all the inputs: time, location, and node
    time = Time(clock.time)
    loc  = location(p.ecco_fts)

    # Retrieve the variable to force
    @inbounds var = fields[i, j, k, p.variable_name]

    ecco_backend = p.ecco_fts.backend
    native_grid = on_native_grid(ecco_backend) 

    ecco_var = get_ecco_variable(Val(native_grid), p.ecco_fts, i, j, k, p.ecco_grid, grid, time)

    # Extracting the mask value at the current node
    mask = stateindex(p.mask, i, j, k, grid, clock.time, loc)

    return p.λ⁻¹ * mask * (ecco_var - var)
end

# Differentiating between restoring done with an ECCO FTS
# that lives on the native ecco grid, that requires interpolation in space
# _inside_ the restoring function and restoring based on an ECCO 
# FTS defined on the model grid that requires only time interpolation
@inline function get_ecco_variable(::Val{true}, ecco_fts, i, j, k, ecco_grid, grid, time)
    # Extracting the ECCO field time series data and parameters
    ecco_times         = ecco_fts.times
    ecco_data          = ecco_fts.data
    ecco_time_indexing = ecco_fts.time_indexing
    ecco_backend       = ecco_fts.backend
    ecco_location      = instantiated_location(ecco_fts)

    X = node(i, j, k, grid, ecco_location...)

    # Interpolating the ECCO field time series data onto the current node and time
    return interpolate(X, time, ecco_data, ecco_location, ecco_grid, ecco_times, ecco_backend, ecco_time_indexing)
end    

@inline get_ecco_variable(::Val{false}, ecco_fts, i, j, k, ecco_grid, grid, time) = @inbounds ecco_fts[i, j, k, time]

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
     metadata = ECCOMetadata(variable_name, all_ecco_dates(version), version)
    return ECCO_restoring_forcing(metadata; kw...)
end

function ECCO_restoring_forcing(metadata::ECCOMetadata;
                                architecture = CPU(), 
                                time_indices_in_memory = 2, # Not more than this if we want to use GPU!
                                time_indexing = Cyclical(),
                                mask = 1,
                                timescale = 20days,
                                grid = nothing)

    ecco_fts  = ECCO_field_time_series(metadata; grid, architecture, time_indices_in_memory, time_indexing)                  
    ecco_grid = ecco_fts.grid

    # Grab the correct Oceananigans field to restore
    variable_name = metadata.name
    field_name = oceananigans_fieldname[variable_name]
    
    ecco_restoring = ECCORestoring(ecco_fts, ecco_grid, mask, field_name, 1 / timescale)
    
    # Defining the forcing that depends on the restoring field.
    restoring_forcing = Forcing(ecco_restoring; discrete_form = true)

    return restoring_forcing
end