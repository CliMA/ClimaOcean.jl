const AA = Oceananigans.Architectures.AbstractArchitecture

JRA55PrescribedAtmosphere(arch::Distributed; kw...) =
    JRA55PrescribedAtmosphere(child_architecture(arch); kw...)

"""
    JRA55PrescribedAtmosphere(architecture::AA;
                              version = JRA55RepeatYear(),
                              dates = all_dates(version),
                              backend = nothing,
                              time_indexing = Cyclical(),
                              reference_height = 10,  # meters
                              include_rivers_and_icebergs = false,
                              other_kw...)

Return a `PrescribedAtmosphere` representing JRA55 reanalysis data.
"""
function JRA55PrescribedAtmosphere(architecture;
                                   version = JRA55RepeatYear(),
                                   dates = all_dates(version),
                                   backend = JRA55NetCDFBackend(10),
                                   time_indexing = Cyclical(),
                                   reference_height = 10,  # meters
                                   include_rivers_and_icebergs = false,
                                   other_kw...)

    kw = (; time_indexing, backend, dates, version)
    kw = merge(kw, other_kw) 

    ua  = JRA55FieldTimeSeries(:eastward_velocity, architecture;               kw...)
    va  = JRA55FieldTimeSeries(:northward_velocity, architecture;              kw...)
    Ta  = JRA55FieldTimeSeries(:temperature, architecture;                     kw...)
    qa  = JRA55FieldTimeSeries(:specific_humidity, architecture;               kw...)
    pa  = JRA55FieldTimeSeries(:sea_level_pressure, architecture;              kw...)
    Fra = JRA55FieldTimeSeries(:rain_freshwater_flux, architecture;            kw...)
    Fsn = JRA55FieldTimeSeries(:snow_freshwater_flux, architecture;            kw...)
    Ql  = JRA55FieldTimeSeries(:downwelling_longwave_radiation, architecture;  kw...)
    Qs  = JRA55FieldTimeSeries(:downwelling_shortwave_radiation, architecture; kw...)

    freshwater_flux = (rain = Fra,
                       snow = Fsn)

    # Remember that rivers and icebergs are on a different grid and have
    # a different frequency than the rest of the JRA55 data. We use `PrescribedAtmospheres`
    # "auxiliary_freshwater_flux" feature to represent them.
    if include_rivers_and_icebergs
        Fri = JRA55FieldTimeSeries(:river_freshwater_flux, architecture;   kw...)
        Fic = JRA55FieldTimeSeries(:iceberg_freshwater_flux, architecture; kw...)
        auxiliary_freshwater_flux = (rivers = Fri, icebergs = Fic)
    else
        auxiliary_freshwater_flux = nothing
    end

    times = ua.times
    grid  = ua.grid

    velocities = (u = ua,
                  v = va)

    tracers = (T = Ta,
               q = qa)

    pressure = pa

    downwelling_radiation = TwoBandDownwellingRadiation(shortwave=Qs, longwave=Ql)

    FT = eltype(ua)
    reference_height = convert(FT, reference_height)

    atmosphere = PrescribedAtmosphere(grid, times;
                                      velocities,
                                      freshwater_flux,
                                      auxiliary_freshwater_flux,
                                      tracers,
                                      downwelling_radiation,
                                      reference_height,
                                      pressure)

    return atmosphere
end
