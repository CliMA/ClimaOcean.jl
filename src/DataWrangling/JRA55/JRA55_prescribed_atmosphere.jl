
const AA = Oceananigans.Architectures.AbstractArchitecture

JRA55PrescribedAtmosphere(time_indices=Colon(); kw...) =
    JRA55PrescribedAtmosphere(CPU(), time_indices; kw...)

JRA55PrescribedAtmosphere(arch::Distributed, time_indices=Colon(); kw...) =
    JRA55PrescribedAtmosphere(child_architecture(arch), time_indices; kw...)

# TODO: allow the user to pass dates
"""
    JRA55PrescribedAtmosphere(architecture::AA, time_indices=Colon();
                              backend = nothing,
                              time_indexing = Cyclical(),
                              dates = all_dates(metadata), 
                              reference_height = 10,  # meters
                              include_rivers_and_icebergs = false,
                              other_kw...)

Return a `PrescribedAtmosphere` representing JRA55 reanalysis data.
"""
function JRA55PrescribedAtmosphere(metadata::JRA55Metadata, arch_or_grid = CPU();
                                   backend = nothing,
                                   time_indexing = Cyclical(),
                                   reference_height = 10,  # meters
                                   include_rivers_and_icebergs = false,
                                   other_kw...)

    time_indices = JRA55_time_indices(metadata.version, metadata.dates)

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

    ua  = JRA55FieldTimeSeries(:eastward_velocity;               kw...)
    va  = JRA55FieldTimeSeries(:northward_velocity;              kw...)
    Ta  = JRA55FieldTimeSeries(:temperature;                     kw...)
    qa  = JRA55FieldTimeSeries(:specific_humidity;               kw...)
    pa  = JRA55FieldTimeSeries(:sea_level_pressure;              kw...)
    Fra = JRA55FieldTimeSeries(:rain_freshwater_flux;            kw...)
    Fsn = JRA55FieldTimeSeries(:snow_freshwater_flux;            kw...)
    Ql  = JRA55FieldTimeSeries(:downwelling_longwave_radiation;  kw...)
    Qs  = JRA55FieldTimeSeries(:downwelling_shortwave_radiation; kw...)

    freshwater_flux = (rain = Fra,
                       snow = Fsn)

    # Remember that rivers and icebergs are on a different grid and have
    # a different frequency than the rest of the JRA55 data. We use `PrescribedAtmospheres`
    # "auxiliary_freshwater_flux" feature to represent them.
    if include_rivers_and_icebergs
        Fri = JRA55FieldTimeSeries(:river_freshwater_flux;   kw...)
        Fic = JRA55FieldTimeSeries(:iceberg_freshwater_flux; kw...)
        auxiliary_freshwater_flux = (rivers = Fri, icebergs = Fic)
    else
        auxiliary_freshwater_flux = nothing
    end

    times = ua.times

    velocities = (u = ua,
                  v = va)

    tracers = (T = Ta,
               q = qa)

    pressure = pa

    downwelling_radiation = TwoBandDownwellingRadiation(shortwave=Qs, longwave=Ql)

    FT = eltype(ua)
    reference_height = convert(FT, reference_height)

    atmosphere = PrescribedAtmosphere(times, FT;
                                      velocities,
                                      freshwater_flux,
                                      auxiliary_freshwater_flux,
                                      tracers,
                                      downwelling_radiation,
                                      reference_height,
                                      pressure)

    return atmosphere
end
