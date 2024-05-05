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
    TwoStreamDownwellingRadiation

using CUDA: @allowscalar

using NCDatasets
using JLD2 
using Dates

import Oceananigans.Fields: set!
import Oceananigans.OutputReaders: new_backend, update_field_time_series!
using Downloads: download

const AA = Oceananigans.Architectures.AbstractArchitecture

JRA55_prescribed_atmosphere(time_indices=Colon(); kw...) =
    JRA55_prescribed_atmosphere(CPU(), time_indices; kw...)

JRA55_prescribed_atmosphere(arch::Distributed, time_indices=Colon(); kw...) =
    JRA55_prescribed_atmosphere(child_architecture(arch), time_indices; kw...)

# TODO: allow the user to pass dates
function JRA55_prescribed_atmosphere(architecture::AA, time_indices=Colon();
                                     backend = nothing,
                                     time_indexing = Cyclical(),
                                     measurement_height = 10,  # meters
                                     other_kw...)

    if isnothing(backend) # apply a default
        Ni = try
            length(time_indices)
        catch
            Inf
        end

        # Manufacture a default for the number of fields to keep InMemory
        Nf = min(24, Ni)
        backend = NetCDFBackend(Nf)
    end

    kw = (; time_indices, time_indexing, backend, architecture)
    kw = merge(kw, other_kw) 

    ua  = field_time_series_from_metadata(JRA55Metadata(:eastward_velocity, time_indices);               kw...)
    va  = field_time_series_from_metadata(JRA55Metadata(:northward_velocity, time_indices);              kw...)
    Ta  = field_time_series_from_metadata(JRA55Metadata(:temperature, time_indices);                     kw...)
    qa  = field_time_series_from_metadata(JRA55Metadata(:specific_humidity, time_indices);               kw...)
    ra  = field_time_series_from_metadata(JRA55Metadata(:relative_humidity, time_indices);               kw...)
    pa  = field_time_series_from_metadata(JRA55Metadata(:sea_level_pressure, time_indices);              kw...)
    Fra = field_time_series_from_metadata(JRA55Metadata(:rain_freshwater_flux, time_indices);            kw...)
    Fsn = field_time_series_from_metadata(JRA55Metadata(:snow_freshwater_flux, time_indices);            kw...)
    Fri = field_time_series_from_metadata(JRA55Metadata(:river_freshwater_flux, time_indices);           kw...)
    Fic = field_time_series_from_metadata(JRA55Metadata(:iceberg_freshwater_flux, time_indices);         kw...)
    Ql  = field_time_series_from_metadata(JRA55Metadata(:downwelling_longwave_radiation, time_indices);  kw...)
    Qs  = field_time_series_from_metadata(JRA55Metadata(:downwelling_shortwave_radiation, time_indices); kw...)

    times = ua.times

    velocities = (u = ua,
                  v = va)

    tracers = (T = Ta,
               q = qa,
               r = ra)

    freshwater_flux = (rain = Fra,
                       snow = Fsn,
                       rivers = Fri,
                       icebergs = Fic)
                       
    pressure = pa

    downwelling_radiation = TwoStreamDownwellingRadiation(shortwave=Qs, longwave=Ql)

    FT = eltype(ua)
    measurement_height = convert(FT, measurement_height)

    atmosphere = PrescribedAtmosphere(times, FT;
                                      velocities,
                                      freshwater_flux,
                                      tracers,
                                      downwelling_radiation,
                                      measurement_height,
                                      pressure)

    return atmosphere
end