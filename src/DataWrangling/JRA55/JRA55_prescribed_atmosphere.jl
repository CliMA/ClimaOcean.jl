const AA = Oceananigans.Architectures.AbstractArchitecture

JRA55PrescribedAtmosphere(arch::Distributed, FT = Float32; kw...) =
    JRA55PrescribedAtmosphere(child_architecture(arch); kw...)


"""
    JRA55PrescribedAtmosphere([architecture = CPU(), FT = Float32];
                              dataset = JRA55RepeatYear(),
                              start_date = first_date(dataset, :temperature),
                              end_date = last_date(dataset, :temperature),
                              backend = JRA55NetCDFBackend(10),
                              time_indexing = Cyclical(),
                              surface_layer_height = 10,  # meters
                              include_rivers_and_icebergs = false,
                              other_kw...)

Return a [`PrescribedAtmosphere`](@ref) representing JRA55 reanalysis data.
The atmospheric data will be held in `JRA55FieldTimeSeries` objects containing.
For a detailed description of the keyword arguments, see the [`JRA55FieldTimeSeries`](@ref) constructor.
"""
function JRA55PrescribedAtmosphere(architecture = CPU(), FT = Float32;
                                   dataset = JRA55RepeatYear(),
                                   start_date = first_date(dataset, :temperature),
                                   end_date = last_date(dataset, :temperature),
                                   backend = InMemory(),
                                   time_indexing = Cyclical(),
                                   surface_layer_height = 10,  # meters
                                   include_rivers_and_icebergs = false,
                                   other_kw...)

    kw = (; time_indexing, backend, start_date, end_date, dataset)
    kw = merge(kw, other_kw) 

    ua  = JRA55FieldTimeSeries(:eastward_velocity, architecture, FT;               kw...)
    va  = JRA55FieldTimeSeries(:northward_velocity, architecture, FT;              kw...)
    Ta  = JRA55FieldTimeSeries(:temperature, architecture, FT;                     kw...)
    qa  = JRA55FieldTimeSeries(:specific_humidity, architecture, FT;               kw...)
    pa  = JRA55FieldTimeSeries(:sea_level_pressure, architecture, FT;              kw...)
    Fra = JRA55FieldTimeSeries(:rain_freshwater_flux, architecture, FT;            kw...)
    Fsn = JRA55FieldTimeSeries(:snow_freshwater_flux, architecture, FT;            kw...)
    Ql  = JRA55FieldTimeSeries(:downwelling_longwave_radiation, architecture, FT;  kw...)
    Qs  = JRA55FieldTimeSeries(:downwelling_shortwave_radiation, architecture, FT; kw...)

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
    surface_layer_height = convert(FT, surface_layer_height)

    atmosphere = PrescribedAtmosphere(grid, times;
                                      velocities,
                                      freshwater_flux,
                                      auxiliary_freshwater_flux,
                                      tracers,
                                      downwelling_radiation,
                                      surface_layer_height,
                                      pressure)

    return atmosphere
end
