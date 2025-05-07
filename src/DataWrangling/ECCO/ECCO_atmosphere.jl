using ClimaOcean.DataWrangling: DatasetBackend
using Oceananigans.OutputReaders
using ClimaOcean.OceanSeaIceModels: TwoBandDownwellingRadiation

"""
    ECCOPrescribedAtmosphere([architecture = CPU(), FT = Float32];
                              dataset = ECCO4Montly(),
                              start_date = first_date(dataset, :temperature),
                              end_date = last_date(dataset, :temperature),
                              backend = DatasetBackend(10),
                              time_indexing = Cyclical(),
                              surface_layer_height = 10,  # meters
                              include_rivers_and_icebergs = false,
                              other_kw...)

Return a [`PrescribedAtmosphere`](@ref) representing JRA55 reanalysis data.
The atmospheric data will be held in `JRA55FieldTimeSeries` objects containing.
For a detailed description of the keyword arguments, see the [`JRA55FieldTimeSeries`](@ref) constructor.
"""
function ECCOPrescribedAtmosphere(architecture = CPU(), FT = Float32;
                                  dataset = ECCO4Monthly(),
                                  start_date = first_date(dataset, :air_temperature),
                                  end_date = last_date(dataset, :air_temperature),
                                  time_indexing = Cyclical(),
                                  time_indices_in_memory = 10,
                                  surface_layer_height = 2,  # meters
                                  other_kw...)

    ua_meta = Metadata(:eastward_wind;         dataset, start_date, end_date)    
    va_meta = Metadata(:northward_wind;        dataset, start_date, end_date)    
    Ta_meta = Metadata(:air_temperature;       dataset, start_date, end_date)    
    qa_meta = Metadata(:air_specific_humidity; dataset, start_date, end_date)    
    pa_meta = Metadata(:sea_level_pressure;    dataset, start_date, end_date)    
    Ql_meta = Metadata(:downwelling_longwave;  dataset, start_date, end_date)
    Qs_meta = Metadata(:downwelling_shortwave; dataset, start_date, end_date)
    Fr_meta = Metadata(:rain_freshwater_flux;  dataset, start_date, end_date)

    kw = (; time_indices_in_memory, time_indexing)
    kw = merge(kw, other_kw)

    ua = FieldTimeSeries(ua_meta, architecture; kw...)
    va = FieldTimeSeries(va_meta, architecture; kw...)
    Ta = FieldTimeSeries(Ta_meta, architecture; kw...)
    qa = FieldTimeSeries(qa_meta, architecture; kw...)
    pa = FieldTimeSeries(pa_meta, architecture; kw...)
    Ql = FieldTimeSeries(Ql_meta, architecture; kw...)
    Qs = FieldTimeSeries(Qs_meta, architecture; kw...)
    Fr = FieldTimeSeries(Fr_meta, architecture; kw...)
    
    auxiliary_freshwater_flux = nothing
    freshwater_flux = (; rain = Fr)

    times = ua.times
    grid  = ua.grid

    velocities = (u = ua, v = va)
    tracers = (T = Ta, q = qa)
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
