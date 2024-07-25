module JRA55

using Oceananigans
using Oceananigans.Units
 
using Oceananigans.Architectures: arch_array
using Oceananigans.DistributedComputations
using Oceananigans.DistributedComputations: child_architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: λnodes, φnodes, on_architecture
using Oceananigans.Fields: interpolate!
using Oceananigans.OutputReaders: Cyclical, TotallyInMemory, AbstractInMemoryBackend, FlavorOfFTS, time_indices

using ClimaOcean.OceanSeaIceModels:
    PrescribedAtmosphere,
    TwoBandDownwellingRadiation

using CUDA: @allowscalar

using NCDatasets
using JLD2 
using Dates
using Scratch

import Oceananigans.Fields: set!
import Oceananigans.OutputReaders: new_backend, update_field_time_series!
using Downloads: download

include("jra55_metadata.jl")

download_jra55_cache::String = ""
function __init__()
    global download_jra55_cache = @get_scratch!("JRA55")
end

compute_bounding_nodes(::Nothing, ::Nothing, LH, hnodes) = nothing
compute_bounding_nodes(bounds, ::Nothing, LH, hnodes) = bounds

# TODO: remove the allowscalar
function compute_bounding_nodes(::Nothing, grid, LH, hnodes)
    hg = hnodes(grid, LH())
    h₁ = @allowscalar minimum(hg)
    h₂ = @allowscalar maximum(hg)
    return h₁, h₂
end

function compute_bounding_indices(::Nothing, hc)
    Nh = length(hc)
    return 1, Nh
end

function compute_bounding_indices(bounds, hc)
    h₁, h₂ = bounds
    Nh = length(hc)

    # The following should work. If ᵒ are the extrema of nodes we want to
    # interpolate to, and the following is a sketch of the JRA55 native grid,
    #
    #      1         2         3         4         5
    # |         |         |         |         |         |
    # |    x  ᵒ |    x    |    x    |    x  ᵒ |    x    |
    # |         |         |         |         |         |
    # 1         2         3         4         5         6 
    #
    # then for example, we should find that (iᵢ, i₂) = (1, 5).
    # So we want to reduce the first index by one, and limit them
    # both by the available data. There could be some mismatch due
    # to the use of different coordinate systems (ie whether λ ∈ (0, 360)
    # which we may also need to handle separately.
    i₁ = searchsortedfirst(hc, h₁)
    i₂ = searchsortedfirst(hc, h₂)
    i₁ = max(1, i₁ - 1)
    i₂ = min(Nh, i₂)

    return i₁, i₂
end

infer_longitudinal_topology(::Nothing) = Periodic

function infer_longitudinal_topology(λbounds)
    λ₁, λ₂ = λbounds
    TX = λ₂ - λ₁ ≈ 360 ? Periodic : Bounded
    return TX
end

function compute_bounding_indices(longitude, latitude, grid, LX, LY, λc, φc)
    λbounds = compute_bounding_nodes(longitude, grid, LX, λnodes)
    φbounds = compute_bounding_nodes(latitude, grid, LY, φnodes)

    i₁, i₂ = compute_bounding_indices(λbounds, λc)
    j₁, j₂ = compute_bounding_indices(φbounds, φc)
    TX = infer_longitudinal_topology(λbounds)

    return i₁, i₂, j₁, j₂, TX
end

# Convert dates to range until Oceananigans supports dates natively
function jra55_times(native_times, start_time=native_times[1])

    times = []
    for native_time in native_times
        time = native_time - start_time
        time = Second(time).value
        push!(times, time)
    end

    return times
end

struct JRA55NetCDFBackend <: AbstractInMemoryBackend{Int}
    start :: Int
    length :: Int
end

"""
    JRA55NetCDFBackend(length)

Represents a JRA55 FieldTimeSeries backed by JRA55 native .nc files.
"""
JRA55NetCDFBackend(length) = JRA55NetCDFBackend(1, length)

Base.length(backend::JRA55NetCDFBackend) = backend.length
Base.summary(backend::JRA55NetCDFBackend) = string("JRA55NetCDFBackend(", backend.start, ", ", backend.length, ")")

const JRA55NetCDFFTS = FlavorOfFTS{<:Any, <:Any, <:Any, <:Any, <:JRA55NetCDFBackend}

function set!(fts::JRA55NetCDFFTS, path::String=fts.path, name::String=fts.name) 

    ds = Dataset(path)

    # Note that each file should have the variables
    #   - ds["time"]:     time coordinate 
    #   - ds["lon"]:      longitude at the location of the variable
    #   - ds["lat"]:      latitude at the location of the variable
    #   - ds["lon_bnds"]: bounding longitudes between which variables are averaged
    #   - ds["lat_bnds"]: bounding latitudes between which variables are averaged
    #   - ds[shortname]:  the variable data

    # Nodes at the variable location
    λc = ds["lon"][:]
    φc = ds["lat"][:]
    LX, LY, LZ = location(fts)
    i₁, i₂, j₁, j₂, TX = compute_bounding_indices(nothing, nothing, fts.grid, LX, LY, λc, φc)

    ti = time_indices(fts)
    ti = collect(ti)
    data = ds[name][i₁:i₂, j₁:j₂, ti]
    close(ds)

    copyto!(interior(fts, :, :, 1, :), data)
    fill_halo_regions!(fts)

    return nothing
end

new_backend(::JRA55NetCDFBackend, start, length) = JRA55NetCDFBackend(start, length)

const AA = Oceananigans.Architectures.AbstractArchitecture

JRA55_prescribed_atmosphere(time_indices=Colon(); kw...) =
    JRA55_prescribed_atmosphere(CPU(), time_indices; kw...)

JRA55_prescribed_atmosphere(arch::Distributed, time_indices=Colon(); kw...) =
    JRA55_prescribed_atmosphere(child_architecture(arch), time_indices; kw...)

# TODO: allow the user to pass dates
function JRA55_prescribed_atmosphere(architecture::AA, time_indices=Colon();
                                     backend = nothing,
                                     time_indexing = Cyclical(),
                                     reference_height = 10,  # meters
                                     include_rivers_and_icebergs = true, # rivers and icebergs are not needed in single column simulations
                                     other_kw...)

    if isnothing(backend) # apply a default
        Ni = try
            length(time_indices)
        catch
            Inf
        end

        # Manufacture a default for the number of fields to keep InMemory
        Nf = min(24, Ni)
        backend = JRA55NetCDFBackend(Nf)
    end

    kw = (; time_indices, time_indexing, backend, architecture)
    kw = merge(kw, other_kw) 

    ua  = JRA55_field_time_series(:eastward_velocity;               kw...)
    va  = JRA55_field_time_series(:northward_velocity;              kw...)
    Ta  = JRA55_field_time_series(:temperature;                     kw...)
    qa  = JRA55_field_time_series(:specific_humidity;               kw...)
    ra  = JRA55_field_time_series(:relative_humidity;               kw...)
    pa  = JRA55_field_time_series(:sea_level_pressure;              kw...)
    Fra = JRA55_field_time_series(:rain_freshwater_flux;            kw...)
    Fsn = JRA55_field_time_series(:snow_freshwater_flux;            kw...)
    Ql  = JRA55_field_time_series(:downwelling_longwave_radiation;  kw...)
    Qs  = JRA55_field_time_series(:downwelling_shortwave_radiation; kw...)

    freshwater_flux = (rain = Fra,
                       snow = Fsn)

    # Remember that rivers and icebergs are on a different grid and have
    # a different frequency than the rest of the JRA55 data
    if include_rivers_and_icebergs
        Fri = JRA55_field_time_series(:river_freshwater_flux;   kw...)
        Fic = JRA55_field_time_series(:iceberg_freshwater_flux; kw...)
        runoff_flux = (rivers   = Fri,
                       icebergs = Fic)
    else
        runoff_flux = nothing
    end

    times = ua.times

    velocities = (u = ua,
                  v = va)

    tracers = (T = Ta,
               q = qa,
               r = ra)

                       
    pressure = pa

    downwelling_radiation = TwoBandDownwellingRadiation(shortwave=Qs, longwave=Ql)

    FT = eltype(ua)
    reference_height = convert(FT, reference_height)

    atmosphere = PrescribedAtmosphere(times, FT;
                                      velocities,
                                      freshwater_flux,
                                      runoff_flux,
                                      tracers,
                                      downwelling_radiation,
                                      reference_height,
                                      pressure)

    return atmosphere
end

end # module

